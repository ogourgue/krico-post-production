"""
krico_recruitment
=================

Post-processing pipeline for KRICO Lagrangian trajectory outputs.
Produces per-particle recruitment outcome tables in NetCDF format.

One row per particle, with outcome classified into one of seven states
(success, censored, or five killed_by categories). See the methodology
document for full definitions.
"""

from krico_recruitment.outcome import OUTCOME_CODES, OUTCOME_MEANINGS
