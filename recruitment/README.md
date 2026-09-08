# Recruitment

Per-particle recruitment post-processing for the KRICO Lagrangian trajectory dataset.

Reads one Parcels trajectory NetCDF per release-day cohort and produces a per-particle recruitment outcome NetCDF following the Thorpe et al. (2019) framework, classifying each particle into one of eight mutually exclusive outcome states.

## Layout

```
recruitment/
├── README.md
├── krico_recruitment/             # Python package
│   ├── __init__.py
│   ├── outcome.py                 # 8-state outcome codes (CF flag)
│   ├── development.py             # T-dependent development (Thorpe Eq. 2, Table 1)
│   ├── filters.py                 # M1 and M4 evaluation
│   ├── sea_ice.py                 # Sea-ice at spawning, advance detection, censoring
│   ├── trajectory.py              # Path length, pick-at-day utilities
│   └── io.py                      # NetCDF read / write
├── scripts/
│   ├── process_cohort.py          # Single-cohort orchestration script
│   ├── run_year.sh                # SLURM array driver (1 task per spawning year)
│   ├── outcome_distribution.py    # Dataset-wide outcome tally
│   ├── outcome_distribution.sh    # SLURM driver for the tally (see Validation)
│   ├── archive_by_year.sh         # Group cohort files into spawning-year tar.gz archives
│   ├── extract_archives.sh        # Extract spawning-year archives back into individual files
│   └── download_from_zenodo.py    # Download archives from Zenodo
├── tests/
├── data/                          # Recruitment outputs (gitignored)
└── archives/                      # Spawning-year tar.gz archives (gitignored)
```

## Two workflow paths

There are two ways to obtain the recruitment data:

### Path A — Full pipeline from raw trajectories

For users with access to the raw KRICO simulation outputs. Reproduces the recruitment classification from scratch.

**Prerequisites:**
- HPC environment (SLURM cluster with Python 3.13.13)
- `KRICO_RUNS` environment variable set to the directory containing raw trajectory simulations (see main [README](../README.md))
- `KRICO_GLORYS12` environment variable set to the GLORYS12 preprocessing output directory containing the monthly sea-ice files (`glorys12_ice_YYYY_MM.nc`). Required because M1 is evaluated at the spawning date, which precedes the trajectory — see [M1 at spawning](#m1-at-spawning) below.

  These files must start at **October 1993**, one month earlier than the trajectories. The first cohort is released on 15 November 1993 and M1 reads back 24 days from there, to 22 October 1993, so `glorys12_ice_1993_10.nc` is required even though no particle is tracked in that month.

**Workflow:**
1. Submit the SLURM array job (one task per spawning year, 32 tasks total):
   ```bash
   cd recruitment/scripts
   sbatch run_year.sh
   ```
   Or for a subset:
   ```bash
   sbatch --array=1 run_year.sh         # First spawning year only (1994)
   sbatch --array=1-5 run_year.sh       # First five spawning years (1994-1998)
   ```
2. Outputs go to `recruitment/data/` (one NetCDF per cohort, ~3,848 files for the full hindcast).

Cohorts whose output already exists are skipped, so re-runs are idempotent. This also means the job does nothing against a populated output directory: when the classification itself changes, move the previous `recruitment/data/` aside first.

For testing on a single cohort:
```bash
python scripts/process_cohort.py \
    $KRICO_RUNS/KRICO_0001/1993_11_15.nc \
    data/1993_11_15.nc \
    --glorys $KRICO_GLORYS12
```

### Path B — Pre-processed data from Zenodo

For users who want to analyze the recruitment outcomes without re-running the pipeline.

**Prerequisites:**
- Local machine (no HPC needed)
- ~100 GB peak disk space (42 GB compressed download + 56 GB extracted; archives can be deleted after extraction to reclaim 42 GB)

**Workflow:**
1. Download the dataset from Zenodo ([10.5281/zenodo.22548763](https://doi.org/10.5281/zenodo.22548763)) into `recruitment/archives/`:
   ```bash
   cd recruitment/scripts
   python download_from_zenodo.py
   ```
2. Extract the archives into `recruitment/data/`:
   ```bash
   ./extract_archives.sh
   ```

After extraction, `recruitment/data/` contains 3,848 NetCDF files (one per release date).

## Dataset versions

| Version | M1 evaluated at | DOI |
|---|---|---|
| v1.0.0 | Release date, from the trajectory at day 0 | [10.5281/zenodo.20101159](https://doi.org/10.5281/zenodo.20101159) |
| **v2.0.0** | Spawning date, 24 days earlier, from the GLORYS12 field | [**10.5281/zenodo.22548763**](https://doi.org/10.5281/zenodo.22548763) |

v2.0.0 is what the current code produces and what Path B downloads by default, so the two paths agree. v1.0.0 remains available for anyone who needs the release-date classification; pass `--doi` to `download_from_zenodo.py` to retrieve it.

The two are not interchangeable. Evaluating M1 at release rather than at spawning underestimates early-season spawning suppression: M1 is 5.19 percentage points lower under v1.0.0, and every outcome downstream of it shifts accordingly. Files from the two versions must not be combined in a single analysis.

Output files from v2.0.0 carry the evaluation provenance as global attributes (`m1_spawning_offset_days`, `m1_spawning_date`, `m1_sic_source`); v1.0.0 files have none, which distinguishes them.

## M1 at spawning

M1 represents a constraint acting on the spawning adult: spawning does not occur where sea-ice concentration is at or above 80%. The 23-26 day descent-ascent cycle from spawning to calyptopis I is not simulated, so particles enter the model as calyptopis I and their release date postdates spawning. Reading sea-ice concentration from the trajectory at day 0 therefore evaluates the constraint at the wrong time, and because concentration falls through the early season it systematically underestimates spawning suppression.

M1 is instead evaluated at the release position on the spawning date, sampling the GLORYS12 sea-ice field directly. The release position stands as the spawning position, consistent with the model's neglect of transport during the descent-ascent interval.

The offset is a constant, `SPAWNING_OFFSET_DAYS = 24` in `krico_recruitment/sea_ice.py`, taken from the 23-26 day range given by Thorpe et al. (2019). Varying it across that range changes the domain-wide M1 fraction by 0.72 percentage points, from 28.11% to 28.83%; see `S1_m1_offset_sensitivity/` in the [krico-paper1](https://github.com/ogourgue/krico-paper1) repo.

Sampling is nearest-neighbour on the reanalysis grid, whereas Parcels interpolates when sampling along trajectories. Compared against the trajectory value at day 0 across the full 32-year record, the two agree to a mean absolute difference of 1.2 × 10⁻⁵ in sea-ice concentration, and 1.8 × 10⁻⁵ of particles are classified differently. The largest disagreements are in coastal cells where the field has a sharp gradient.

M1 is the only filter affected. M4 acts from calyptopis I onward, which is the release, and M5 and M6 are triggered by calendar events during tracking.

## Archiving (project maintainer)

To re-create the spawning-year tar.gz archives from `recruitment/data/` (e.g., before uploading to Zenodo):

```bash
cd recruitment/scripts
./archive_by_year.sh
```

This produces 32 archives in `recruitment/archives/` (one per spawning year, 1994–2025), each containing 120 cohort files — 121 in the eight leap years, which is how 32 spawning years give 3,848 cohorts.

## Output schema

One row per particle. See the Paper 1 methodology document for full column definitions.

| Variable | Type | Description |
|---|---|---|
| `particle` (coord) | int64 | Particle id from upstream trajectory file |
| `release_lon`, `release_lat`, `release_depth` | float32 | Release position |
| `outcome` | int8 (flag) | 8-state outcome (see below) |
| `final_lon`, `final_lat`, `final_depth` | float32 | Position at moment fate is determined |
| `travel_time` | int32 | Days from release to moment fate is determined |
| `trajectory_path_length` | float32 | Along-track path length (m) up to moment fate is determined |

Outcome flag codes:

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | censored |
| 2 | killed_M1 |
| 3 | killed_M4 |
| 4 | killed_M5_no_FIV |
| 5 | killed_M5_not_on_shelf |
| 6 | killed_M6_no_advance |
| 7 | exited_domain |

`censored` is a subcategory of `killed_M6_no_advance`: both refer to particles that survived the calyptopis window and were alive at end of tracking with no detected sea-ice advance event. The `censored` label additionally marks particles whose classification is provisional because SIC was rising in the last 30 days of tracking, suggesting an advance event might have followed shortly after the cutoff. The two are kept as separate flag codes here so downstream analyses can distinguish them; figures in the Paper 1 repo fold censored into killed_M6 by default.

## End-of-tracking tolerance

A particle is considered alive at end of tracking if its last valid trajectory day is within `END_OF_TRACKING_TOLERANCE = 1` day of the final tracking day. Parcels can produce a NaN position on the very last output step due to time-interpolation precision at the simulation end, even for particles that were never deleted. Without this tolerance, such particles would be misclassified as `exited_domain` with their fate position recorded one day before the end of tracking — which for ACC-drifting particles can place them deep in the domain interior rather than at a model boundary.

For end-of-tracking outcomes (`censored`, `killed_M6_no_advance`), `fate_day` is set to the particle's last valid trajectory day (typically `n_obs - 1`, but `n_obs - 2` for particles affected by the precision artifact). Position lookup at `fate_day` therefore returns valid coordinates in both cases.

The constant lives in `krico_recruitment/sea_ice.py`.

## Dependencies

Python 3.13.13, `numpy`, `xarray`, `pandas`, `netCDF4`. Parcels is not required.

## Validation

Dataset-wide outcome distribution over the full 32-year run (1994–2025, 3 848 cohorts, ≈ 2.1 × 10⁹ particles). Counts and fractions are for v2.0.0, the current classification with M1 at spawning; the last column gives the v1.0.0 fraction, with M1 at release, for comparison.

| outcome | count (v2.0.0) | fraction (v2.0.0) | fraction (v1.0.0) |
|---|---:|---:|---:|
| success | 110,417,165 | 5.26% | 5.83% |
| censored | 27,221,228 | 1.30% | 1.42% |
| killed_M1 | 595,485,067 | 28.35% | 23.16% |
| killed_M4 | 134,504,220 | 6.40% | 8.46% |
| killed_M5_no_FIV | 131,306,367 | 6.25% | 7.32% |
| killed_M5_not_on_shelf | 282,050,117 | 13.43% | 14.37% |
| killed_M6_no_advance | 803,211,965 | 38.24% | 38.67% |
| exited_domain | 16,069,207 | 0.77% | 0.78% |

Total: 2,100,265,336 particles. `killed_M6_no_advance` and `censored` together account for 39.54%.

Evaluating M1 at spawning rather than at release raises it by 5.19 percentage points, drawn from the categories downstream of it roughly in proportion to their size. `exited_domain` is unchanged, as it must be: domain exit is determined by the trajectory, not by M1.

Regenerate with:
```bash
cd recruitment/scripts
sbatch outcome_distribution.sh          # writes to recruitment/data/
```

Per-cohort qualitative behaviors match Thorpe (2019): success rate peaks in late January; M1 declines as ice retreats through the season; M4 peaks in late summer; M5_no_FIV rises monotonically into March as the time available to reach FIV shrinks. See `F2_phenology_curve` and `F3_outcome_composition` in the [krico-paper1](https://github.com/ogourgue/krico-paper1) repo for the visual breakdown.
