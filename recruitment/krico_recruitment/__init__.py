"""
krico_recruitment
=================

Post-processing pipeline for KRICO Lagrangian trajectory outputs.
Produces per-particle recruitment outcome tables in NetCDF format.

One row per particle: the outcome, one of nine states (two classes of
recruits, M1, M3, M4, M5, M6, censored, exited_domain; see
krico_recruitment.outcome), and the diagnostics it is derived from (see
krico_recruitment.classify and krico_recruitment.io).
"""

from krico_recruitment.outcome import OUTCOME_CODES, OUTCOME_MEANINGS