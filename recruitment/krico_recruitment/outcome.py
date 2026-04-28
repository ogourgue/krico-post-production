"""
Outcome categories and integer flag encoding.

Each particle ends in exactly one of eight mutually exclusive states.
Stored as a CF-compliant integer flag in NetCDF output.
"""

# CF-compliant flag_values and flag_meanings for NetCDF encoding.
# The order defines the integer codes (0..7).
OUTCOME_MEANINGS = (
    "success",
    "censored",
    "killed_M1",
    "killed_M4",
    "killed_M5_no_FIV",
    "killed_M5_not_on_shelf",
    "killed_M6_no_advance",
    "exited_domain",
)

# Mapping from string label to integer code (for readable code).
OUTCOME_CODES = {label: code for code, label in enumerate(OUTCOME_MEANINGS)}

# Reverse mapping for inspection / debugging.
OUTCOME_LABELS = {code: label for code, label in enumerate(OUTCOME_MEANINGS)}