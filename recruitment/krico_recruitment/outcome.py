"""
Outcome categories and integer flag encoding (dataset v3.0.0).

Each particle ends in exactly one of nine mutually exclusive states, stored
as a CF flag. Recruits come first, then the rules in the order they act.
The rule labels follow Table 2 of Thorpe et al. (2019).

Recruitment has two levels. A larva that passes M1, M3 and M4 and has
reached furcilia IV when sea ice advances over it is a recruit (codes 0
and 1), as in Thorpe et al. (2019). The two codes then say where the
recruit is at that moment: over the shelf or slope (sea floor shallower
than 2000 m) or beyond it. That distinction is ours and is not a mortality.

The codes were renumbered from v2.0.0, whose eight labels are kept in
V2_MEANINGS with the correspondence V2_TO_V3. Code 3 (M3) is new.
"""

import numpy as np

# CF flag_values / flag_meanings. The order defines the integer codes (0..8).
OUTCOME_MEANINGS = (
    "recruit_on_shelf_slope",
    "recruit_off_shelf_slope",
    "M1_spawning_under_ice",
    "M3_cold_sensitivity",
    "M4_calyptopis_starvation",
    "M5_overwintering_stage",
    "M6_overwintering_habitat",
    "censored",
    "exited_domain",
)

# Mapping from label to integer code, and back.
OUTCOME_CODES = {label: code for code, label in enumerate(OUTCOME_MEANINGS)}
OUTCOME_LABELS = {code: label for code, label in enumerate(OUTCOME_MEANINGS)}

# Labels of dataset v2.0.0, in the order of their codes (0..7).
V2_MEANINGS = (
    "success",
    "censored",
    "killed_M1",
    "killed_M4",
    "killed_M5_no_FIV",
    "killed_M5_not_on_shelf",
    "killed_M6_no_advance",
    "exited_domain",
)

# v3.0.0 code of each v2.0.0 code: V2_TO_V3[v2_code]. With M3 switched off,
# a v2.0.0 outcome and a v3.0.0 outcome are the same state under this map.
V2_TO_V3 = np.array([
    OUTCOME_CODES["recruit_on_shelf_slope"],    # 0 success
    OUTCOME_CODES["censored"],                  # 1 censored
    OUTCOME_CODES["M1_spawning_under_ice"],     # 2 killed_M1
    OUTCOME_CODES["M4_calyptopis_starvation"],  # 3 killed_M4
    OUTCOME_CODES["M5_overwintering_stage"],    # 4 killed_M5_no_FIV
    OUTCOME_CODES["recruit_off_shelf_slope"],   # 5 killed_M5_not_on_shelf
    OUTCOME_CODES["M6_overwintering_habitat"],  # 6 killed_M6_no_advance
    OUTCOME_CODES["exited_domain"],             # 7 exited_domain
], dtype=np.int8)
