# recruitment

Per-particle recruitment post-processing for the KRICO Lagrangian trajectory dataset.

Reads one Parcels trajectory NetCDF per release-day cohort and produces a per-particle recruitment outcome NetCDF following the schema defined in the Paper 1 methodology document.

## Layout

```
recruitment/
├── README.md
├── krico_recruitment/          # python package
│   ├── __init__.py
│   ├── outcome.py              # 7-state outcome codes (CF flag)
│   ├── development.py          # T-dependent development (Thorpe Eq. 2, Table 1)
│   ├── filters.py              # M1 and M4 evaluation
│   ├── sea_ice.py              # sea-ice advance detection, censoring
│   ├── trajectory.py           # path length, pick-at-day utilities
│   └── io.py                   # NetCDF read / write
├── scripts/
│   └── process_cohort.py       # single-cohort orchestration script
└── data/                       # recruitment outputs (one .nc per cohort)
```

## Single-cohort usage

```bash
python scripts/process_cohort.py <input_trajectory.nc> <output_recruitment.nc>
```

Example:
```bash
python scripts/process_cohort.py \
    /scratch/cvan/KRICO/Runs/KRICO_0001/1993_11_15.nc \
    data/1993_11_15.nc
```

## Output schema

One row per particle. See the Paper 1 methodology document for full column definitions.

| Variable | Type | Description |
|---|---|---|
| `particle` (coord) | int64 | Particle id from upstream trajectory file |
| `release_lon`, `release_lat`, `release_depth` | float32 | Release position |
| `outcome` | int8 (flag) | 7-state outcome (see below) |
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

## Dependencies

Python ≥ 3.11, `numpy`, `xarray`, `pandas`, `netCDF4`.

## Validation

First quantitative validation: 5 cohorts processed across the 1993-94 spawning year (Nov 15 → Mar 14). Results:
- Success-rate bell shape centered on Jan 15: 2.2% → 6.9% → **8.3%** → 7.8% → 4.1%, matching Thorpe (2019).
- All filter-specific behaviors (M1 declining with ice retreat, M4 peaking in Feb, M5_no_FIV monotonically rising into Mar) match Thorpe's qualitative predictions.
- Genuine `killed_M6_no_advance` is essentially zero (~0.01%) — virtually all particles that survive M1, M4, and stay in the domain encounter sea-ice advance.

Quantitative equivalence with the PoC is not expected because:
- M4 is now T-dependent (PoC used fixed 33-day window).
- M5 is applied as a developmental threshold (PoC did not).
- Winter onset is per-particle sea-ice advance (PoC had no equivalent in its unfiltered matrix, and used fixed 15 May equivalent in filtered curves).
- `exited_domain` is now a separate outcome (PoC implicitly mixed it with other categories).