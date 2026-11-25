# Recruitment

Per-particle recruitment post-processing for the KRICO Lagrangian trajectory dataset.

Reads one Parcels trajectory NetCDF per release-day cohort and produces a per-particle NetCDF following the Thorpe et al. (2019) framework: the outcome of each particle, one of nine mutually exclusive states, and the diagnostics the outcome is derived from.

> **State of this branch.** The code is that of dataset **v3.0.0**, in production: [pass 0](#release-cell-extraction-glorys12-sea-ice-and-temperature) (inputs at the release cells) and [pass 1](#dataset-v300-pass-1) (the classification) are written; pass 2 (chlorophyll) is not. The published dataset is still v2.0.0, which the tag `v2.0.1` produces; sections that describe it say so.

## Layout

```
recruitment/
├── README.md
├── krico_recruitment/             # Python package
│   ├── __init__.py
│   ├── outcome.py                 # 9-state outcome codes (CF flag), v2.0.0 correspondence
│   ├── development.py             # T-dependent development (Thorpe Eq. 2, Table 1), rate change
│   ├── filters.py                 # M1, M3 and M4 thresholds, runs under sea ice
│   ├── pre_release.py             # Before release: sea ice at spawning, M3 temperatures, from the release cells
│   ├── classify.py                # Diagnostics along the trajectory, outcome from the diagnostics
│   ├── inputs.py                  # Checks of the inputs of one cohort
│   ├── sea_ice.py                 # Sea-ice at spawning, advance detection, censoring
│   ├── trajectory.py              # Path length, pick-at-day utilities
│   ├── io.py                      # NetCDF read / write, schema of the v3.0.0 files
│   ├── nsidc_grid.py              # NSIDC 25 km polar-stereographic grid geometry
│   ├── nsidc_ice.py               # NSIDC sea ice at spawning and along trajectories
│   ├── release_cells.py           # Release cells on the GLORYS12 grid, particle-to-cell link
│   └── provenance.py              # Commit and version written into the output files
├── scripts/
│   ├── process_cohort.py          # Pass 1: classification of one cohort
│   ├── run_year.sh                # SLURM array driver (1 task per spawning year)
│   ├── check_inputs.py / .sh      # Inputs of pass 1, checked for every cohort before submission
│   ├── validate_pass1.py / .sh    # Regression against v2.0.0, report on written files
│   ├── run_year_nsidc.sh          # v2 only, not ported: NSIDC sea ice (see Ice-source test)
│   ├── run_year_glorys12_25km.sh  # v2 only, not ported: GLORYS12 ice on the NSIDC grid
│   ├── download_nsidc.py / .sh    # Daily NSIDC G02202 v6 files (Python + SLURM array)
│   ├── nsidc_daily_provenance.txt # DOI, URL, access date, sha256 of every daily file
│   ├── download_pisces.py / .sh   # Daily 3D chlorophyll of the PISCES hindcast (Copernicus Marine toolbox)
│   ├── download_globcolour.py / .sh  # Daily chlorophyll climatology of GlobColour (Copernicus Marine toolbox)
│   ├── preprocess_nsidc.py / .sh  # Daily files -> monthly stacks on the native grid
│   ├── preprocess_glorys12_25km.py / .sh  # GLORYS12 ice binned onto the NSIDC grid
│   ├── preprocess_glorys12_release_cells.py / .sh  # GLORYS12 ice and temperature at the release cells
│   ├── glorys12_release_cells.csv # The 35,303 release cells (grid indices, position, depth)
│   ├── validate_glorys12_release_cells.py / .sh  # Acceptance tests of that extraction
│   ├── validate_nsidc.py          # NSIDC vs GLORYS12 ice readings on given cohorts
│   ├── validate_nsidc.csv         # The four-cohort validation behind run_year_nsidc.sh
│   ├── outcome_distribution.py    # Dataset-wide outcome tally
│   ├── outcome_distribution.sh    # SLURM driver for the tally (see Validation)
│   ├── archive_by_year.sh         # Group cohort files into spawning-year tar.gz archives
│   ├── extract_archives.sh        # Extract spawning-year archives back into individual files
│   └── download_from_zenodo.py    # Download archives from Zenodo
├── data/                          # Recruitment outputs, v2.0.0 (gitignored)
├── data_v3.0_glorys12/            # Pass 1 outputs, v3.0.0 (gitignored)
├── data_nsidc/                    # Recruitment outputs with NSIDC sea ice (gitignored)
├── data_glorys12_25km/            # Same, with GLORYS12 ice on the NSIDC grid (gitignored)
└── archives/                      # Spawning-year tar.gz archives (gitignored)
```

## Two workflow paths

There are two ways to obtain the recruitment data:

### Path A — Full pipeline from raw trajectories

For users with access to the raw KRICO simulation outputs. Reproduces the classification of v3.0.0 from scratch.

**Prerequisites:**
- HPC environment (SLURM cluster with Python 3.13.13)
- `KRICO_RUNS` set to the directory containing raw trajectory simulations (see main [README](../README.md))
- `KRICO_GLORYS12_RELEASE_CELLS` set to the directory of the release-cell files, `glorys12_release_cells_YYYY.nc`. These hold the GLORYS12 sea ice and temperature at the release cells and are produced by pass 0 from the GLORYS12 forcing (`KRICO_GLORYS12`), see [Release-cell extraction](#release-cell-extraction-glorys12-sea-ice-and-temperature). Pass 1 itself reads no forcing file.

**Workflow** (from `recruitment/scripts`; details in [Dataset v3.0.0: pass 1](#dataset-v300-pass-1)):
```bash
mkdir -p logs ../data_v3.0_glorys12
sbatch check_inputs.sh                               # every cohort through the checks of the job
sbatch --array=1 run_year.sh data_v3.0_glorys12      # first spawning year (1994)
sbatch --array=2-32 run_year.sh data_v3.0_glorys12   # the rest
```
Outputs go to the folder named on the command line (one NetCDF per cohort, 3,848 files for the full hindcast). There is no default folder. Cohorts whose output already exists are skipped, so a failed task is simply resubmitted.

To classify a single cohort:
```bash
python3 scripts/process_cohort.py \
    $KRICO_RUNS/KRICO_0001/1993_11_15.nc \
    data_v3.0_glorys12/1993_11_15.nc
```

### Path B — Pre-processed data from Zenodo

For users who want to analyze the recruitment outcomes without re-running the pipeline. This path still serves **v2.0.0**, the published dataset; v3.0.0 is not on Zenodo yet.

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

| v3.0.0 (in production) | As v2.0.0; adds M3, recruitment in two classes, the diagnostics | not yet published |

v2.0.0 is what the tag `v2.0.1` produces and what Path B downloads by default. The code of this branch produces v3.0.0, whose outcome codes differ (see [Dataset v3.0.0: pass 1](#dataset-v300-pass-1)). v1.0.0 remains available for anyone who needs the release-date classification; pass `--doi` to `download_from_zenodo.py` to retrieve it.

The two are not interchangeable. Evaluating M1 at release rather than at spawning underestimates early-season spawning suppression: M1 is 5.19 percentage points lower under v1.0.0, and every outcome downstream of it shifts accordingly. Files from the two versions must not be combined in a single analysis.

Output files from v2.0.0 carry the evaluation provenance as global attributes (`m1_spawning_offset_days`, `m1_spawning_date`, `m1_sic_source`); v1.0.0 files have none, which distinguishes them.

## M1 at spawning

M1 represents a constraint acting on the spawning adult: spawning does not occur where sea-ice concentration is at or above 80%. The 23-26 day descent-ascent cycle from spawning to calyptopis I is not simulated, so particles enter the model as calyptopis I and their release date postdates spawning. Reading sea-ice concentration from the trajectory at day 0 therefore evaluates the constraint at the wrong time, and because concentration falls through the early season it systematically underestimates spawning suppression.

M1 is instead evaluated at the release position on the spawning date, sampling the GLORYS12 sea-ice field directly. The release position stands as the spawning position, consistent with the model's neglect of transport during the descent-ascent interval.

The offset is a constant, `SPAWNING_OFFSET_DAYS = 24` in `krico_recruitment/sea_ice.py`, taken from the 23-26 day range given by Thorpe et al. (2019). Varying it across that range changes the domain-wide M1 fraction by 0.72 percentage points, from 28.11% to 28.83%; see `S1_m1_offset_sensitivity/` in the [krico-paper1](https://github.com/ogourgue/krico-paper1) repo.

Sampling is nearest-neighbour on the reanalysis grid. Parcels does not interpolate a tracer in space either: on the C-grid a particle reads the value of the cell that holds it, and for particles within float32 rounding of a cell edge that cell can be the neighbour of the nearest grid point (see [Release-cell extraction](#release-cell-extraction-glorys12-sea-ice-and-temperature)). Compared against the trajectory value at day 0 across the full 32-year record, the two readings agree to a mean absolute difference of 1.2 × 10⁻⁵ in sea-ice concentration, and 1.8 × 10⁻⁵ of particles are classified differently. The largest disagreements are where the field has a sharp gradient.

M1 is the only filter affected. M4 acts from calyptopis I onward, which is the release, and M5 and M6 are triggered by calendar events during tracking.

## Ice-source test (NSIDC sea ice)

GLORYS12 is forced by ERA-Interim up to 31 December 2018 and by ERA5 from 1 January 2019. To test whether that switch leaves a mark in the ice-dependent outcomes, the classification can be rerun on the same trajectories with every sea-ice concentration it reads, M1 at spawning and the along-trajectory series behind M4, the advance detection and censoring, taken from the NOAA/NSIDC Climate Data Record of Passive Microwave Sea Ice Concentration, version 6 (G02202, doi:10.7265/b18j-z797) instead of GLORYS12. The trajectories, and so the circulation behind M5b, are unchanged.

Inputs are built once, into `$KRICO_NSIDC` (e.g. `/scratch/cvan/KRICO/Pre/NSIDC`):

```bash
cd scripts && mkdir -p logs
sbatch download_nsidc.sh                 # daily files, 1993-10-01 to 2025-10-31, ~4 GB
python3 download_nsidc.py --provenance   # when all tasks are done; commit the provenance file
python3 preprocess_nsidc.py --check-grid # grid geometry against the ancillary lat/lon
sbatch preprocess_nsidc.sh               # monthly stacks nsidc_ice_YYYY_MM.nc, 385 files
```

Then `sbatch run_year_nsidc.sh` runs the classification into `data_nsidc/`, with `--ice-source nsidc` passed to `process_cohort.py`. Output files carry `sic_source = "NSIDC G02202 v6 (...)"` and `nsidc_nearest_ocean`; the production files carry `sic_source = "GLORYS12"` (v2.0.0 files predate the attribute and have only `m1_sic_source`). The two sets must not be combined.

A third classification separates resolution from product in the comparison: `preprocess_glorys12_25km.py` bins the GLORYS12 monthly ice onto the NSIDC grid (plain mean of the ORCA12 cells inside each 25 km cell, NSIDC land/coast/lake as fill) into `$KRICO_GLORYS12_25KM`, and `run_year_glorys12_25km.sh` runs the classification on it through the same sampler (`--ice-source glorys12_25km`) into `data_glorys12_25km/`, with `sic_source = "GLORYS12v1 ice binned onto the NSIDC G02202 25 km grid"`. Native GLORYS12 minus this run is the resolution and sampling effect; this run minus NSIDC is the product difference at equal resolution. On the full record the resolution effect is at or below 0.1 percentage points for every outcome domain-wide, so the whole GLORYS12-versus-NSIDC difference is in the ice values.

Sampling reads the 25 km cell containing each position, with no interpolation (`krico_recruitment/nsidc_ice.py`). Coast and land cells are fill in the CDR and read NaN, which the filters treat as below threshold; a nearest-ocean redirection exists (`--nsidc-nearest-ocean`) but is not used, because on four cohorts spanning the record coast cells hold 0.3-0.4% of positions and move M1 by at most 0.4 points. On those cohorts the two readings agree along trajectories to r = 0.97-0.98 and a mean absolute difference of 0.035-0.044 in concentration, M1 differs by 2-4 points in either direction, and advance detection by under 1 point (`scripts/validate_nsidc.csv`). The 2.2% of release positions north of the NSIDC grid edge (50-51°S, 8°W-7°E, in 48.6N) read NaN and get the same outcome in both products, since neither has ice there.

## Release-cell extraction (GLORYS12 sea ice and temperature)

The rules that act before release, in the 24 days that are not simulated, depend only on the date and on the grid cell of the release position. Their inputs are extracted once, without reading any trajectory: `preprocess_glorys12_release_cells.py` copies the daily GLORYS12 sea-ice concentration and the temperature column at the release cells into one file per spawning year, `$KRICO_GLORYS12_RELEASE_CELLS/glorys12_release_cells_YYYY.nc` (e.g. `/scratch/cvan/KRICO/Pre/GLORYS12_RELEASE_CELLS`). The files are an intermediate product, like the trajectories and the forcing: they are not part of the published dataset. The extraction evaluates no rule.

- **Cells:** the 35,303 cells in which `run.py` of krico-templates releases particles (`scripts/glorys12_release_cells.csv`: sea floor of 1000 to 2000 m, centre inside CCAMLR subareas 48.1 to 48.6 and 88.3) and their eight neighbours, 55,034 cells. A particle is linked to the grid point nearest to its release position, the sampling rule of M1, and that point can be a neighbour of the cell the particle was generated in (`krico_recruitment/release_cells.py`). The table was built in [krico-paper2](https://github.com/ogourgue/krico-paper2) (`local/subarea/raw/spawning_cells.csv`, sha256 `6c4f0e01f954…`) and reproduces the 546,199 particles generated per release day.
- **Period:** 22 October (first spawning date) to 30 September (last tracked day) of each spawning year, daily.
- **Levels:** surface to the sea floor of the deepest release cell, 40 levels; the level bounds (`depthw`) are stored with them.
- **Values:** copied unchanged as float32, with lossless compression; 2.1 GB on disk and 15 minutes for spawning year 1994.
- **Time stamps:** every sea-ice record of the source from 19 December 2018 is stamped 11:52:30 instead of 12:00. It is taken as the mean of its day, and each file summarises such records in `source_stamps_not_at_noon`.

```bash
cd scripts && mkdir -p logs
python3 preprocess_glorys12_release_cells.py --check             # inputs of the 32 years, no output
sbatch validate_glorys12_release_cells.sh coverage ../data       # acceptance test 1, coverage
sbatch --array=1 preprocess_glorys12_release_cells.sh            # spawning year 1994, as a test
sbatch validate_glorys12_release_cells.sh day0 \
    $KRICO_RUNS/KRICO_0001/1993_11_15.nc \
    $KRICO_RUNS/KRICO_0003/1994_01_15.nc \
    $KRICO_RUNS/KRICO_0004/1994_03_14.nc                         # acceptance test 2, day 0
sbatch --array=2-32 preprocess_glorys12_release_cells.sh         # the other years
```

`validate_glorys12_release_cells.py` runs two tests, whose reports are the job outputs `logs/cells_val_<job id>.out`. Coverage (`--releases`): the nearest grid point of every release position of the dataset is among the extracted cells. Day 0 (trajectory files): the sea ice and temperature a trajectory recorded at its release equal the extraction on the release day, at the cell and level of the particle; Parcels samples tracers on the C-grid without spatial interpolation and particles are released at 12:00, the time stamp of the daily means, so the comparison is exact, and comparisons shifted by a day or a level must be worse.

Both tests passed on 25 November 2026. Coverage, over the 2,100,265,336 release positions of the v2.0.0 dataset: none has its nearest grid point outside the extracted cells; 242,288 have it in a cell south of a release cell, 1,955 in a cell west of one and 4 in another neighbour. Day 0, on the cohorts of 15 November 1993, 15 January 1994 and 14 March 1994: at the nearest grid point, temperature is exactly equal for 99.975 to 99.982% of particles and sea ice for 99.987 to 99.997%; one day or one level away, at most 0.007% of temperatures are equal. The particles not reproduced exactly lie within float32 rounding of a cell edge or a level bound, where the stored position cannot tell which cell Parcels used. For cohorts released from 19 December 2018, sea ice at day 0 is close to the extraction but not exactly equal, because Parcels interpolates in time between the shifted stamps.

Every file records the commit that produced it (`code_commit`, `code_version`); the extraction stops if the code is modified, untracked or not pushed.

## Dataset v3.0.0: pass 1

Pass 1 (`process_cohort.py`) reads one trajectory file and the release-cell file of its spawning year, and writes one file per cohort into `recruitment/data_v3.0_glorys12/`.

### What changed from v2.0.0

- **M3, cold sensitivity, is applied**, after M1 and before M4. Like M1 it acts during the 24 days between spawning and release, which are not simulated, and is evaluated at the release cell.
- **Recruitment has two levels.** A larva that passes M1, M3 and M4 and is at furcilia IV when sea ice advances over it is a recruit, as in Thorpe et al. (2019). Recruits are then split by where they are at that moment: over the shelf or slope (sea floor shallower than 2000 m) or beyond it. The former `killed_M5_not_on_shelf` is the second class of recruits, not a mortality.
- **The outcome codes are renumbered** (table below; `krico_recruitment.outcome.V2_TO_V3` gives the correspondence).
- **Diagnostics are stored for every particle**, including those removed by M1 or M3, and the outcome is computed from them. Any threshold can therefore be evaluated from the files, without the trajectories, with the function that produced the dataset:
  ```python
  from krico_recruitment import classify, io
  d, coordinates = io.read_diagnostics("data_v3.0_glorys12/1993_11_15.nc")
  code = classify.outcome(d, classify.Parameters(m1_sic=0.70), **coordinates)
  ```
- With M3 switched off, the classification is that of v2.0.0, particle for particle (`validate_pass1.py regression`).

| Code | Outcome | v2.0.0 code |
|---|---|---|
| 0 | `recruit_on_shelf_slope` | 0 `success` |
| 1 | `recruit_off_shelf_slope` | 5 `killed_M5_not_on_shelf` |
| 2 | `M1_spawning_under_ice` | 2 |
| 3 | `M3_cold_sensitivity` | new |
| 4 | `M4_calyptopis_starvation` | 3 |
| 5 | `M5_overwintering_stage` | 4 `killed_M5_no_FIV` |
| 6 | `M6_overwintering_habitat` | 6 |
| 7 | `censored` | 1 |
| 8 | `exited_domain` | 7 |

M5 and M6 are evaluated at each particle's own sea-ice advance, detected from 1 April, not at a fixed date. A particle with no advance is M6 (or censored) without its stage being tested.

### M3

Thorpe et al. (2019) remove embryos and larvae whose mean temperature is −1 °C or colder during either the descent or the ascent. Days are counted from the spawning day, day 0: the descent takes the daily means of days 0 to 4, the ascent those of days 5 to 23, and release is day 24. Their embryo sinks at a density-dependent rate and hatches at 530 to 1000 m; the KRICO forcing has no salinity, so:

- the hatching depth is bracketed by 1000 m (baseline) and 530 m, both stored;
- sinking is steady, so the descent mean is the thickness-weighted mean of the water column from the surface to the hatching depth, over the descent days. The ascent mean is the same column mean over the ascent days.

The column is cut at the sea floor of the cell; a cell with no value removes nothing. M2 (hatching on the sea floor) holds by construction: embryos hatch at 1000 m or shallower and the release cells are 1000 m or deeper. Details in `krico_recruitment/pre_release.py`.

### Variables

All on the `particle` dimension. Days are day indices along the trajectory (0 = release), `uint8`, 255 for none. Development is in thousandths of a stage from the start of calyptopis I, rounded down, `uint16`, 65535 for none (3000 = furcilia I, 6000 = furcilia IV).

| Variable | Type | Serves |
|---|---|---|
| `release_lon`, `release_lat`, `release_depth`, `outcome`, `final_lon`, `final_lat`, `final_depth`, `travel_time`, `trajectory_path_length` | as in v2.0.0 | A particle removed by M3 is recorded like one removed by M1: release position, day 0 |
| `release_cell` | int32 | Index of the particle's cell in the release-cell files |
| `sic_spawning` | float32 | M1 at any threshold |
| `m3_descent_temperature`, `m3_ascent_temperature` (particle, hatching_depth) | float32 | M3 at any threshold, for hatching at 530 or 1000 m |
| `spawning_temperature_surface`, `_mean_100m`, `_mean_200m` | float32 | Warm limit at spawning (a sensitivity, not in `outcome`) |
| `m4_longest_run` (particle, ice_threshold) | uint8 | Longest run under ice above 0.20, 0.30, 0.40, 0.50, 0.60 during the calyptopis stages: M4 at those thresholds, any duration |
| `m4_first_run_day`, `m4_first_run_dev` | day, development | First run of more than 10 days above 0.40 over the whole track: M4 under another development rate |
| `calyptopis_end_day`, `fiv_day` | day | First day at furcilia I, at furcilia IV |
| `advance_day`, `advance_dev`, `advance_bathymetry`, `advance_lon`, `advance_lat` | day, development, float32 | The state at the advance: any stage, any shelf/slope depth, any development rate |
| `alive_at_end`, `sic_rising_at_end` | int8 | M6, censoring, domain exit; not masked by the earlier rules |
| `release_bathymetry` | float32 | Sea floor at release |
| `days_below_fit` | uint8 | Days below −1 °C, the lower limit of the fitted range of the development function |
| `ice15_daily`, `shelf2000_daily` (particle, packed_day) | uint8 | One bit per day, `numpy.unpackbits(x, axis=1)[:, :201]`: sea ice above 0.15, sea floor shallower than 2000 m. Another definition of the advance; history on and off the shelf/slope |

Each diagnostic carries as attributes the parameters baked into it; `outcome` carries the baseline parameters. Global attributes hold the provenance (`code_commit`, `code_version`, `dataset_version`, `trajectory_source`, `release_cells_source`, `history`) and the result of the day-0 check: the share of particles whose sea ice and temperature at day 0 equal the release-cell file exactly, the share within 0.01, and the mean and largest difference.

A development rate other than 1 is evaluated by rescaling the development clock (`development.apply_rate`), which is exact up to the daily step of the integration: against a full re-integration on synthetic cohorts, about 0.05% of particles change outcome. It cannot be combined with another M4 threshold or duration.

### Tests

```bash
cd recruitment/scripts
sbatch check_inputs.sh                                     # before anything else
sbatch validate_pass1.sh regression $KRICO_RUNS/KRICO_0001/1993_11_15.nc
sbatch validate_pass1.sh cohort $KRICO_RUNS/KRICO_0001/1993_11_15.nc data_v3.0_test
```

- `check_inputs` runs on every cohort the checks the job itself applies, and lists every problem.
- `regression` classifies a cohort in memory with M3 off and requires the v2.0.0 file in `recruitment/data/`: the outcome through the code table, the eight other variables value for value. With `KRICO_GLORYS12` set it also compares the sea ice at spawning with the monthly files. It then reports what M3 removes.
- `cohort` writes one file into a test folder and reports on it: outcome recomputed from the stored diagnostics, cross-checks between diagnostics, bytes on disk per variable, run time and peak memory.

### Not ported

`run_year_nsidc.sh` and `run_year_glorys12_25km.sh` stop with a message: `process_cohort.py` reads GLORYS12 sea ice only. The classifications of the [ice-source test](#ice-source-test-nsidc-sea-ice) were produced with the v2 code.

## Archiving (project maintainer)

To re-create the spawning-year tar.gz archives from `recruitment/data/` (e.g., before uploading to Zenodo):

```bash
cd recruitment/scripts
./archive_by_year.sh
```

This produces 32 archives in `recruitment/archives/` (one per spawning year, 1994–2025), each containing 120 cohort files — 121 in the eight leap years, which is how 32 spawning years give 3,848 cohorts.

## Output schema (v2.0.0)

The schema of v3.0.0 is in [Dataset v3.0.0: pass 1](#dataset-v300-pass-1). One row per particle. See the Paper 1 methodology document for full column definitions.

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

Python 3.13.13, `numpy`, `xarray`, `pandas`, `netCDF4`; `h5py` for the bytes per variable in `validate_pass1.py report`. Parcels is not required.

## Validation

The numbers of this section are those of v2.0.0 and v1.0.0; those of v3.0.0 follow its full run.

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
