"""
Compare the NSIDC and GLORYS12 sea-ice readings on a few cohorts, before the
NSIDC classification is run on the whole record.

Usage
-----
    python3 validate_nsidc.py --nsidc $KRICO_NSIDC --glorys $KRICO_GLORYS12 \\
        [--stride 10] [--csv validate_nsidc.csv] cohort1.nc cohort2.nc ...

--stride N keeps every N-th particle. A cohort is ~545k particles x 201
obs, and this script holds three along-trajectory readings at once, which
exceeds what a login node allows; stride 10 is plenty to characterise the
comparison and runs in a few hundred MB.

For each trajectory file this reports, on stdout and as one CSV row:

  timing      whether obs k of the trajectory falls on the calendar day
              release_date + k, the convention the sampler relies on
              (checked on the trajectory's own time variable; must be 0
              mismatches; the time of day, noon in the KRICO output, is
              reported but irrelevant to a daily field)
  release     where the release positions fall on the NSIDC grid: ocean,
              near-coast (ocean, 1-3 cells from land), coast, land, lake,
              off-grid. Coast/land/lake read NaN in the CDR.
  track       the same for every valid trajectory position
  values      on positions where both readings are valid: bias (NSIDC minus
              GLORYS12), mean absolute difference, correlation, and the
              agreement of the two readings on the 15% (advance) and 40% (M4)
              thresholds
  M1          share killed at spawning (>= 80%) under each reading
  advance     share of particles with a detected advance, and the median
              difference in advance day, under each reading

The NSIDC columns are given twice: as read (coast/land -> NaN) and with
nearest-ocean redirection, so that the effect of the coast cells can be
judged before deciding on --nsidc-nearest-ocean for the full run.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import filters, nsidc_ice, sea_ice, trajectory  # noqa: E402
from krico_recruitment import io as krico_io                        # noqa: E402

SURFACE = {50: "ocean", 75: "lake", 200: "coast", 250: "land", 0: "offgrid"}


def _surface_shares(surface, adj, mask):
    """Share of positions (within mask) per surface class, near-coast split out."""
    n = int(mask.sum())
    out = {}
    for code, name in SURFACE.items():
        sel = mask & (surface == code)
        if code == 50:
            near = sel & (adj >= 1) & (adj <= 3)
            out["ocean"] = (sel & ~near).sum() / n
            out["nearcoast"] = near.sum() / n
        else:
            out[name] = sel.sum() / n
    return out


def validate(path, nsidc_dir, glorys_dir, stride=1):
    ds = krico_io.read_trajectory_file(path)
    sel = slice(None, None, stride)
    lon = ds["lon"][sel].values
    lat = ds["lat"][sel].values
    depth = ds["z"][sel].values
    sic_g = ds["sea_ice_area_fraction"][sel].values
    time_ns = ds["time"][sel].values
    n_particles, n_obs = lon.shape
    release_date = pd.Timestamp(time_ns[0, 0]).normalize()
    row = {"cohort": Path(path).name, "release_date": str(release_date.date()),
           "n_particles": n_particles, "stride": stride}

    # Timing convention: obs k falls on the calendar day release + k.
    day_ns = 86_400_000_000_000
    t = time_ns.astype("datetime64[ns]").astype("int64")
    finite_t = np.isfinite(lon) & np.isfinite(lat)
    day_index = np.floor_divide(t - release_date.value, day_ns)
    off = (day_index != np.arange(n_obs)[None, :]) & finite_t & (t > 0)
    row["timing_mismatches"] = int(off.sum())
    tod = np.mod(t - release_date.value, day_ns)[finite_t & (t > 0)]
    row["timing_time_of_day_h"] = float(np.median(tod) / 3.6e12) if tod.size else np.nan

    # Positions the NSIDC grid cannot place: NaN at release, and finite
    # positions outside the grid rectangle (north of ~50-55S depending on
    # longitude). The latter have no ice in either product.
    row["release_nan"] = float((~finite_t[:, 0]).mean())
    r0, c0 = nsidc_ice.cell_index(lon[:, 0], lat[:, 0])
    off0 = finite_t[:, 0] & (r0 < 0)
    row["release_offgrid_finite"] = float(off0.mean())
    row["release_offgrid_lat_max"] = float(lat[:, 0][off0].max()) if off0.any() else np.nan
    row["release_offgrid_lat_min"] = float(lat[:, 0][off0].min()) if off0.any() else np.nan
    row["release_offgrid_lon_min"] = float(lon[:, 0][off0].min()) if off0.any() else np.nan
    row["release_offgrid_lon_max"] = float(lon[:, 0][off0].max()) if off0.any() else np.nan

    # Where positions fall on the NSIDC grid.
    surface, adj = nsidc_ice.cell_surface_type(lon, lat, nsidc_dir)
    for k, v in _surface_shares(surface[:, 0], adj[:, 0],
                                np.ones(n_particles, bool)).items():
        row[f"release_{k}"] = v
    for k, v in _surface_shares(surface, adj, finite_t).items():
        row[f"track_{k}"] = v

    # Along-trajectory readings.
    sic_n = nsidc_ice.along_trajectory_sic(lon, lat, release_date, nsidc_dir)
    sic_no = nsidc_ice.along_trajectory_sic(lon, lat, release_date, nsidc_dir,
                                            nearest_ocean=True)
    valid_g = np.isfinite(sic_g) & finite_t
    row["track_glorys_nan"] = 1 - valid_g.sum() / finite_t.sum()
    row["track_nsidc_nan"] = 1 - (np.isfinite(sic_n) & finite_t).sum() / finite_t.sum()
    row["track_nsidc_no_nan"] = 1 - (np.isfinite(sic_no) & finite_t).sum() / finite_t.sum()

    both = valid_g & np.isfinite(sic_n)
    d = sic_n[both] - sic_g[both]
    row["bias"] = float(d.mean())
    row["mad"] = float(np.abs(d).mean())
    row["corr"] = float(np.corrcoef(sic_n[both], sic_g[both])[0, 1])
    for thr, name in ((sea_ice.ADVANCE_SIC_THRESHOLD, "15"), (0.40, "40")):
        row[f"agree_{name}pct"] = float(((sic_n[both] > thr) == (sic_g[both] > thr)).mean())

    # M1 under each reading.
    m1_g = filters.evaluate_M1(sea_ice.spawning_sic(lon[:, 0], lat[:, 0],
                                                    release_date, glorys_dir))
    m1_n = filters.evaluate_M1(nsidc_ice.spawning_sic(lon[:, 0], lat[:, 0],
                                                      release_date, nsidc_dir))
    m1_no = filters.evaluate_M1(nsidc_ice.spawning_sic(lon[:, 0], lat[:, 0],
                                                       release_date, nsidc_dir,
                                                       nearest_ocean=True))
    row["m1_glorys"] = float(m1_g.mean())
    row["m1_nsidc"] = float(m1_n.mean())
    row["m1_nsidc_no"] = float(m1_no.mean())
    row["m1_agree"] = float((m1_g == m1_n).mean())

    # Advance under each reading.
    april1 = sea_ice.find_april1_day_index(release_date, n_obs)
    adv_g = sea_ice.detect_sea_ice_advance(sic_g, april1)
    adv_n = sea_ice.detect_sea_ice_advance(sic_n, april1)
    adv_no = sea_ice.detect_sea_ice_advance(sic_no, april1)
    row["adv_glorys"] = float((adv_g >= 0).mean())
    row["adv_nsidc"] = float((adv_n >= 0).mean())
    row["adv_nsidc_no"] = float((adv_no >= 0).mean())
    bothadv = (adv_g >= 0) & (adv_n >= 0)
    row["adv_both"] = float(bothadv.mean())
    row["adv_day_median_diff"] = float(np.median(adv_n[bothadv] - adv_g[bothadv])) if bothadv.any() else np.nan
    row["adv_day_mad"] = float(np.abs(adv_n[bothadv] - adv_g[bothadv]).mean()) if bothadv.any() else np.nan

    # Not strictly needed, but cheap: particles alive at the end.
    last_valid = trajectory.last_valid_day(lon, lat, depth)
    row["alive_at_end"] = float((last_valid >= n_obs - 2).mean())
    return row


def _print(row):
    pct = lambda v: f"{100 * v:6.2f}%"  # noqa: E731
    print(f"\n{row['cohort']}  released {row['release_date']}  "
          f"{row['n_particles']} particles" + (f" (stride {row['stride']})" if row['stride'] > 1 else ""))
    print(f"  timing   : {row['timing_mismatches']} obs off the daily convention "
          f"(time of day {row['timing_time_of_day_h']:.1f} h)")
    print(f"  release  : NaN position {pct(row['release_nan'])}; finite but off the NSIDC grid "
          f"{pct(row['release_offgrid_finite'])}"
          + (f", lat {row['release_offgrid_lat_min']:.1f} to {row['release_offgrid_lat_max']:.1f}, "
             f"lon {row['release_offgrid_lon_min']:.1f} to {row['release_offgrid_lon_max']:.1f}"
             if row['release_offgrid_finite'] > 0 else ""))
    for scope in ("release", "track"):
        print(f"  {scope:8s} : " + "  ".join(
            f"{k} {pct(row[f'{scope}_{k}'])}"
            for k in ("ocean", "nearcoast", "coast", "land", "lake", "offgrid")))
    print(f"  track NaN: GLORYS {pct(row['track_glorys_nan'])}  NSIDC {pct(row['track_nsidc_nan'])}"
          f"  NSIDC/nearest-ocean {pct(row['track_nsidc_no_nan'])}")
    print(f"  values   : bias {row['bias']:+.3f}  MAD {row['mad']:.3f}  r {row['corr']:.3f}  "
          f"agree>15% {pct(row['agree_15pct'])}  agree>40% {pct(row['agree_40pct'])}")
    print(f"  M1 killed: GLORYS {pct(row['m1_glorys'])}  NSIDC {pct(row['m1_nsidc'])}  "
          f"NSIDC/nearest-ocean {pct(row['m1_nsidc_no'])}  per-particle agreement {pct(row['m1_agree'])}")
    print(f"  advance  : GLORYS {pct(row['adv_glorys'])}  NSIDC {pct(row['adv_nsidc'])}  "
          f"NSIDC/nearest-ocean {pct(row['adv_nsidc_no'])}  both {pct(row['adv_both'])}  "
          f"day diff median {row['adv_day_median_diff']:+.0f}  MAD {row['adv_day_mad']:.1f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("cohorts", nargs="+", help="trajectory files (YYYY_MM_DD.nc)")
    parser.add_argument("--nsidc", required=True)
    parser.add_argument("--glorys", required=True)
    parser.add_argument("--csv", default=None, help="write one row per cohort")
    parser.add_argument("--stride", type=int, default=1,
                        help="keep every N-th particle (default 1 = all)")
    args = parser.parse_args()

    rows = []
    for path in args.cohorts:
        row = validate(path, args.nsidc, args.glorys, args.stride)
        _print(row)
        rows.append(row)
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"\nWritten {args.csv}")


if __name__ == "__main__":
    main()
