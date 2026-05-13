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
│   ├── sea_ice.py                 # Sea-ice advance detection, censoring
│   ├── trajectory.py              # Path length, pick-at-day utilities
│   └── io.py                      # NetCDF read / write
├── scripts/
│   ├── process_cohort.py          # Single-cohort orchestration script
│   ├── run_year.sh                # SLURM array driver (1 task per spawning year)
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
- HPC environment (SLURM cluster with Python 3.11+)
- `KRICO_RUNS` environment variable set to the directory containing raw trajectory simulations (see main [README](../README.md))

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

For testing on a single cohort:
```bash
python scripts/process_cohort.py \
    $KRICO_RUNS/KRICO_0001/1993_11_15.nc \
    data/1993_11_15.nc
```

### Path B — Pre-processed data from Zenodo

For users who want to analyze the recruitment outcomes without re-running the pipeline.

**Prerequisites:**
- Local machine (no HPC needed)
- ~100 GB peak disk space (42 GB compressed download + 56 GB extracted; archives can be deleted after extraction to reclaim 42 GB)

**Workflow:**
1. Download the dataset from Zenodo (10.5281/zenodo.20101159) into `recruitment/archives/`:
   ```bash
   cd recruitment/scripts
   python download_from_zenodo.py
   ```
2. Extract the archives into `recruitment/data/`:
   ```bash
   ./extract_archives.sh
   ```

After extraction, `recruitment/data/` contains 3,848 NetCDF files (one per release date), identical to what Path A would produce.

## Archiving (project maintainer)

To re-create the spawning-year tar.gz archives from `recruitment/data/` (e.g., before uploading to Zenodo):

```bash
cd recruitment/scripts
./archive_by_year.sh
```

This produces 32 archives in `recruitment/archives/` (one per spawning year, 1994–2025), each containing ~121 cohort files.

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

`censored` is a subcategory of `killed_M6_no_advance`: both refer to particles that survived the calyptope window and were alive at end of tracking with no detected sea-ice advance event. The `censored` label additionally marks particles whose classification is provisional because SIC was rising in the last 30 days of tracking, suggesting an advance event might have followed shortly after the cutoff. The two are kept as separate flag codes here so downstream analyses can distinguish them; figures in the Paper 1 repo fold censored into killed_M6 by default.

## End-of-tracking tolerance

A particle is considered alive at end of tracking if its last valid trajectory day is within `END_OF_TRACKING_TOLERANCE = 1` day of the final tracking day. Parcels can produce a NaN position on the very last output step due to time-interpolation precision at the simulation end, even for particles that were never deleted. Without this tolerance, such particles would be misclassified as `exited_domain` with their fate position recorded one day before the end of tracking — which for ACC-drifting particles can place them deep in the domain interior rather than at a model boundary.

For end-of-tracking outcomes (`censored`, `killed_M6_no_advance`), `fate_day` is set to the particle's last valid trajectory day (typically `n_obs - 1`, but `n_obs - 2` for particles affected by the precision artifact). Position lookup at `fate_day` therefore returns valid coordinates in both cases.

The constant lives in `krico_recruitment/sea_ice.py`.

## Dependencies

Python ≥ 3.11, `numpy`, `xarray`, `pandas`, `netCDF4`.

## Validation

Dataset-wide outcome distribution over the full 32-year run (1994–2025, 3 848 cohorts, ≈ 2.1 × 10⁹ particles):

| outcome | fraction |
|---|---|
| success | 5.83% |
| censored | 1.42% |
| killed_M1 | 23.16% |
| killed_M4 | 8.46% |
| killed_M5_no_FIV | 7.32% |
| killed_M5_not_on_shelf | 14.37% |
| killed_M6_no_advance | 38.67% |
| exited_domain | 0.78% |

Per-cohort qualitative behaviors match Thorpe (2019): success rate peaks in mid-January; M1 declines as ice retreats through the season; M4 peaks in late summer; M5_no_FIV rises monotonically into March as the time available to reach FIV shrinks. See `F2_phenology_curve` and `F3_outcome_composition` in the [krico-paper1](https://github.com/ogourgue/krico-paper1) repo for the visual breakdown.