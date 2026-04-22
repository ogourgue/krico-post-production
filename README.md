# KRICO/Post/Production

Post-processing and analysis pipeline for the KRICO project (Antarctic Krill Connectivity). Based on 30-year Lagrangian particle tracking simulations (1994–2025) over CCAMLR Areas 48 and 88, forced by GLORYS12v1 ocean reanalysis.

**Author:** Olivier Gourgue (RBINS)

---

## Repository structure

```
krico-post-production/
├── config.py                  # Paths and shared constants
└── ...                        # Analysis scripts (to be added)
```

Each analysis subfolder contains scripts and a `data/` subdirectory (gitignored) for outputs.

---

## Setup

### 1. Clone the repository

```bash
cd /path/to/KRICO/Post/
git clone https://github.com/ogourgue/krico-post-production.git Production
```

### 2. Set the environment variable

Point the scripts to your local copy of the KRICO data:

```bash
export KRICO_ROOT="/path/to/KRICO"
```

Add this line to your shell configuration file (e.g. `~/.bashrc` or `~/.zshrc`) to make it persistent across sessions.

---

## Simulations

128 simulation folders (`KRICO_0001` to `KRICO_0128`) covering 32 spawning years (1994–2025), with 4 release month folders per year (November, December, January, February). Each folder contains daily NetCDF output files (~1.6–1.7 GB each) with ~500,000 particle trajectories tracked over 200 days.
