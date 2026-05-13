# KRICO: Post-processing

Post-processing and analysis pipeline for the KRICO project (Antarctic Krill Connectivity). Based on a 32-year Lagrangian particle tracking hindcast (1994–2025) of Antarctic krill larval dispersal across CCAMLR Areas 48 and 88, forced by GLORYS12v1 ocean reanalysis.

Author: Olivier Gourgue (RBINS)

Related repositories:

* __[krico-templates](https://github.com/ogourgue/krico-templates)__ — Simulation templates (Parcels + GLORYS12v1)
* __[krico-paper1](https://github.com/ogourgue/krico-paper1)__ — Paper 1 figure reproduction

## Repository structure

```
krico-post-production/
├── README.md
├── recruitment/                   # Per-particle recruitment classification
│   ├── README.md
│   ├── krico_recruitment/         # Python package
│   ├── scripts/                   # Processing and data management scripts
│   ├── data/                      # Recruitment NetCDF files (gitignored)
│   ├── archives/                  # Tar.gz archives (gitignored)
│   └── tests/
└── ...                            # Additional analysis subfolders (to be added)
```

Each analysis subfolder is self-contained, with its own scripts, data, and documentation. Data and archives folders are gitignored and created automatically by the scripts.

## Setup

Clone the repository wherever you like:

```bash
git clone https://github.com/ogourgue/krico-post-production.git
```

Some analyses require access to raw trajectory data, in which case the `KRICO_RUNS` environment variable must be set to point to the simulation root:

```bash
export KRICO_RUNS="/path/to/raw/trajectories"      # e.g., /scratch/cvan/KRICO/Runs
```

Add this line to your shell configuration file (e.g. `~/.bash_profile`, `~/.bashrc`, or `~/.zshrc`) to make it persistent across sessions. See each analysis subfolder's README for whether `KRICO_RUNS` is required.

The raw trajectory data are not currently publicly archived due to storage constraints. They are available on demand for the duration of the KRICO project (until end of 2026), after which public archiving will be revisited.

## Analyses

* __[recruitment/](recruitment/)__ — Per-particle classification into 8 outcome states (success, censored, killed by mortality filters M1/M4/M5/M6, exited domain). Outputs available on Zenodo (DOI: [10.5281/zenodo.20101159](https://doi.org/10.5281/zenodo.20101159)).

Additional analysis subfolders will be added as the project progresses.