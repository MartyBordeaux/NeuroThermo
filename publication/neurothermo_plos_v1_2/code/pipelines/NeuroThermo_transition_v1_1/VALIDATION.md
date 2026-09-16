# Validation requirements

The package validates the following before producing results:

1. the input contains the completed transition v1.0 files;
2. `transition_paths.csv` has 121,524 rows;
3. there are 988 support-state scenarios and 72 biological WT×SCA3 pairs;
4. the frozen q=.75 reference contains 18 primary cells and 12 `core_q75_secure` cells;
5. the v1.0 ISI projection is reproduced numerically to machine precision;
6. `support_rate_hz` is copied exactly to `active_support_rate_hz`, not recomputed;
7. the best endpoint active-window rate agrees with the independently frozen dynamic-v2.1 model active rate within the expected small numerical discrepancy;
8. the corrected experimental active-rate WT/SCA3 reference clouds do not overlap under the stringent endpoint rule.

This stage performs no HR integration and no refitting.
