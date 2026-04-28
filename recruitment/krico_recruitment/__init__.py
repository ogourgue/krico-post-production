"""
krico_recruitment
=================

Post-processing pipeline for KRICO Lagrangian trajectory outputs.
Produces per-particle recruitment outcome tables in NetCDF format.

One row per particle, with outcome classified into one of eight states:
success, censored, five killed_by categories (M1, M4, M5_no_FIV,
M5_not_on_shelf, M6_no_advance), or exited_domain. See the methodology
document for full definitions.
"""

from krico_recruitment.outcome import OUTCOME_CODES, OUTCOME_MEANINGS