"""
Tests of pass 1 (process_cohort.py) on real cohorts.

Usage
-----
    # 1. Regression against v2.0.0, before the tag and the full run.
    #    Classifies each cohort in memory and writes nothing.
    python3 validate_pass1.py regression \\
        $KRICO_RUNS/KRICO_0001/1993_11_15.nc $KRICO_RUNS/KRICO_0003/1994_01_15.nc

    # 2. Report on files written by process_cohort.py: the one-cohort test.
    python3 validate_pass1.py report ../data_v3.0_test/1993_11_15.nc

Regression
----------
With M3 switched off, the new classification must return the v2.0.0 one:

  - outcome, through the code table outcome.V2_TO_V3, for every particle;
  - the eight other variables of v2.0.0, value for value;
  - the particle identifiers.

M1 is now read from the release-cell file instead of the monthly sea-ice
files, M4 comes from run lengths instead of a loop, and the outcome is
derived from the stored diagnostics instead of the arrays in memory: the
test shows that none of this changed a single particle. It then switches
M3 on and reports what M3 removes and from which v2.0.0 outcome.

The v2.0.0 files are read from recruitment/data/ (--v2 DIR for another
folder).

Report
------
For each file written by process_cohort.py:

  - the outcome recomputed from the diagnostics read back from the file
    equals the stored one (classify.outcome);
  - cross-checks between diagnostics that the classification computes by
    different routes: the advance recomputed from the daily sea-ice bits,
    the shelf bit on the advance day against the sea floor at advance, the
    day furcilia IV is reached against the development at advance;
  - the bytes each variable takes on disk, and per particle (needs h5py).

Where KRICO_GLORYS12 is set, the regression also compares the sea ice at
spawning of every particle with the reading of v2.0.0 from the monthly
sea-ice files, value for value.

Environment
-----------
    KRICO_GLORYS12_RELEASE_CELLS   required for the regression
    KRICO_GLORYS12                 optional; monthly sea-ice files
"""

import argparse
import dataclasses
import os
import sys
import time
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import classify, development, sea_ice, io as krico_io  # noqa: E402
from krico_recruitment.outcome import (OUTCOME_CODES, OUTCOME_MEANINGS,  # noqa: E402
                                       V2_MEANINGS, V2_TO_V3)
import process_cohort  # noqa: E402

V2_VARIABLES = ("release_lon", "release_lat", "release_depth", "final_lon",
                "final_lat", "final_depth", "travel_time",
                "trajectory_path_length")
DEFAULT_V2 = Path(__file__).resolve().parents[1] / "data"


# --------------------------------------------------------------- regression

def regression(trajectory_path, release_cells_dir, v2_dir):
    started = time.time()
    trajectory_path = Path(trajectory_path)
    v2_path = Path(v2_dir) / trajectory_path.name
    print(f"== {trajectory_path.name}: regression against {v2_path}")
    if not v2_path.is_file():
        print("  verdict: FAIL, the v2.0.0 file is missing")
        return False
    no_m3 = dataclasses.replace(classify.BASELINE, m3_temperature=None)
    particle_id, new, info = process_cohort.classify_cohort(
        trajectory_path, release_cells_dir, parameters=no_m3,
        log=lambda *a: None)
    n = particle_id.size
    with xr.open_dataset(v2_path) as ds:
        old = {name: ds[name].values for name in V2_VARIABLES + ("outcome", "particle")}

    reasons = []
    if old["particle"].shape != particle_id.shape:
        print(f"  verdict: FAIL, {old['particle'].size} particles in v2.0.0, "
              f"{n} in the trajectory file")
        return False
    if not np.array_equal(old["particle"], particle_id):
        reasons.append("particle identifiers differ")

    expected = V2_TO_V3[old["outcome"].astype(np.int64)]
    differ = expected != new["outcome"]
    print(f"  outcome, M3 off: {int(differ.sum())} of {n} particles differ "
          "from v2.0.0")
    if differ.any():
        reasons.append(f"outcome differs for {int(differ.sum())} particles")
        pairs, counts = np.unique(
            np.stack([old["outcome"][differ], new["outcome"][differ]]),
            axis=1, return_counts=True)
        for (a, b), c in sorted(zip(pairs.T, counts), key=lambda x: -x[1])[:10]:
            print(f"    v2.0.0 {V2_MEANINGS[a]} -> {OUTCOME_MEANINGS[b]}: {c}")
    for name in V2_VARIABLES:
        a, b = old[name], new[name]
        if a.dtype != b.dtype:
            reasons.append(f"{name}: type {b.dtype}, v2.0.0 has {a.dtype}")
        equal = (a == b) | (np.isnan(a) & np.isnan(b)) if a.dtype.kind == "f" else a == b
        bad = int((~equal).sum())
        print(f"  {name:24s} {bad} of {n} values differ")
        if bad:
            reasons.append(f"{name} differs for {bad} particles")

    # The outcome written in a file is the baseline, with M3.
    with_m3 = classify.outcome(new)
    removed = with_m3 == OUTCOME_CODES["M3_cold_sensitivity"]
    print(f"  M3 on: removes {int(removed.sum())} particles "
          f"({100 * removed.mean():.2f}%), which were in v2.0.0:")
    for code, label in enumerate(V2_MEANINGS):
        c = int((old["outcome"][removed] == code).sum())
        if c:
            print(f"    {label:24s} {c:8d}")
    untouched = ~removed & (with_m3 != new["outcome"])
    if untouched.any():
        reasons.append(f"M3 changed {int(untouched.sum())} particles it does "
                       "not remove")
    recruits = with_m3 <= OUTCOME_CODES["recruit_off_shelf_slope"]
    recruits_v2 = np.isin(old["outcome"], (0, 5))
    print(f"  recruits (both classes): {100 * recruits_v2.mean():.2f}% without "
          f"M3, {100 * recruits.mean():.2f}% with M3")
    for name, s in info["day0"].items():
        print(f"  day 0, {name}: {100 * s['exact_share']:.4f}% exactly equal, "
              f"largest difference {s['max_abs_difference']:.3e}")

    # Sea ice at spawning, value for value, against the reading of v2.0.0
    # from the monthly files (sea_ice.spawning_sic): the direct test that
    # the release-cell file and the particle-to-cell link reproduce M1.
    glorys_dir = os.environ.get("KRICO_GLORYS12")
    if glorys_dir:
        reference = sea_ice.spawning_sic(
            new["release_lon"], new["release_lat"], info["release_date"],
            glorys_dir).astype(np.float32)
        ours = new["sic_spawning"]
        bad = int((~((reference == ours)
                     | (np.isnan(reference) & np.isnan(ours)))).sum())
        print(f"  sic_spawning against the monthly files: {bad} of {n} "
              "values differ")
        if bad:
            reasons.append(f"sic_spawning differs from the monthly files for "
                           f"{bad} particles")
    else:
        print("  KRICO_GLORYS12 not set: sic_spawning not compared with the "
              "monthly files")
    ok = not reasons
    print("  verdict: " + ("PASS" if ok else "FAIL, " + "; ".join(reasons))
          + f"  ({time.time() - started:.0f} s)")
    return ok


# ------------------------------------------------------------------- report

def _storage(path):
    """Bytes on disk per variable, or None without h5py."""
    try:
        import h5py
    except ImportError:
        return None
    with h5py.File(path, "r") as f:
        return {name: f[name].id.get_storage_size() for name in f.keys()}


def _unpack(packed, n_days):
    return np.unpackbits(packed, axis=1)[:, :n_days].astype(bool)


def report(path):
    path = Path(path)
    print(f"== {path}")
    d, coordinates = krico_io.read_diagnostics(path)
    with netCDF4.Dataset(path) as nc:
        attrs = {k: nc.getncattr(k) for k in nc.ncattrs()}
        n_days = int(nc["ice15_daily"].n_days)
    n = d["outcome"].size
    reasons = []
    print(f"  {n} particles, release {attrs.get('release_date')}, "
          f"code {attrs.get('code_version')} ({str(attrs.get('code_commit'))[:7]}), "
          f"dataset {attrs.get('dataset_version')}")

    expected = krico_io.variable_names(n_days)
    names = [v for v in d if v != "particle"]
    if names != expected:
        reasons.append("variables differ from the schema: "
                       f"{sorted(set(names) ^ set(expected))}")

    # Outcome from the diagnostics, as any reader would recompute it.
    again = classify.outcome(d, **coordinates)
    bad = int((again != d["outcome"]).sum())
    print(f"  outcome recomputed from the stored diagnostics: {bad} of {n} differ")
    if bad:
        reasons.append(f"outcome is not reproduced for {bad} particles")
    for value, label in enumerate(OUTCOME_MEANINGS):
        c = int((d["outcome"] == value).sum())
        print(f"    {value} {label:26s} {c:8d}  ({100.0 * c / n:5.2f}%)")

    # Cross-checks between diagnostics.
    ice15 = _unpack(d["ice15_daily"], n_days)
    shelf = _unpack(d["shelf2000_daily"], n_days)
    release = np.datetime64(str(attrs["release_date"])[:10])
    opening = sea_ice.find_april1_day_index(str(release), n_days)
    advance = sea_ice.detect_sea_ice_advance(
        np.where(ice15, 1.0, 0.0).astype(np.float32), opening)
    advance = np.where(advance == -1, classify.DAY_NONE, advance)
    bad = int((advance != d["advance_day"]).sum())
    print(f"  advance recomputed from ice15_daily (opening on day {opening}): "
          f"{bad} differ")
    if bad:
        reasons.append(f"advance_day and ice15_daily disagree for {bad} particles")
    has = d["advance_day"] != classify.DAY_NONE
    rows = np.flatnonzero(has)
    with np.errstate(invalid="ignore"):
        bad = int((shelf[rows, d["advance_day"][rows]]
                   != (d["advance_bathymetry"][rows] < classify.SHELF_DEPTH)).sum())
    print(f"  shelf2000_daily on the advance day against advance_bathymetry: "
          f"{bad} differ")
    if bad:
        reasons.append(f"shelf2000_daily and advance_bathymetry disagree for {bad}")
    known = has & (d["advance_dev"] != development.DEV_MISSING)
    at_fiv = d["advance_dev"] >= development.FIV_INDEX * development.DEV_SCALE
    by_day = (d["fiv_day"] != classify.DAY_NONE) & (d["fiv_day"] <= d["advance_day"])
    bad = int((at_fiv != by_day)[known].sum())
    print(f"  furcilia IV at advance, advance_dev against fiv_day: {bad} differ")
    if bad:
        reasons.append(f"advance_dev and fiv_day disagree for {bad} particles")

    for name in ("sea_ice", "temperature"):
        print(f"  day 0, {name}: "
              f"{100 * float(attrs[f'day0_{name}_exact_share']):.4f}% exactly "
              f"equal, largest difference "
              f"{float(attrs[f'day0_{name}_max_abs_difference']):.3e}")

    sizes = _storage(path)
    total = path.stat().st_size
    if sizes is None:
        print(f"  h5py not available: {total} bytes in all, "
              f"{total / n:.2f} per particle")
    else:
        print("  bytes on disk per variable:")
        for name in ["particle"] + names:
            raw = d[name].nbytes
            print(f"    {name:32s} {sizes[name]:11d}  {sizes[name] / n:6.2f} "
                  f"per particle  (uncompressed {raw / n:5.1f})")
        v2 = sum(sizes[v] for v in V2_VARIABLES + ("outcome", "particle"))
        print(f"    {'the variables of v2.0.0':32s} {v2:11d}  {v2 / n:6.2f} per particle")
        print(f"    {'diagnostics':32s} {sum(sizes.values()) - v2:11d}  "
              f"{(sum(sizes.values()) - v2) / n:6.2f} per particle")
        print(f"    {'file':32s} {total:11d}  {total / n:6.2f} per particle")
    ok = not reasons
    print("  verdict: " + ("PASS" if ok else "FAIL, " + "; ".join(reasons)))
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("test", choices=("regression", "report"))
    parser.add_argument("files", nargs="+", type=Path,
                        help="trajectory files (regression) or pass 1 files (report)")
    parser.add_argument("--v2", type=Path, default=DEFAULT_V2, metavar="DIR",
                        help="folder of the v2.0.0 files (default: recruitment/data)")
    args = parser.parse_args()
    ok = True
    try:
        if args.test == "regression":
            directory = os.environ.get("KRICO_GLORYS12_RELEASE_CELLS")
            if not directory:
                sys.exit("ERROR: KRICO_GLORYS12_RELEASE_CELLS environment "
                         "variable not set.")
            for path in args.files:
                ok &= regression(path, directory, args.v2)
        else:
            for path in args.files:
                ok &= report(path)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")
    print("ALL PASS" if ok else "AT LEAST ONE FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
