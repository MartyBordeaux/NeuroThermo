# Animal-ID recovery provenance

Animal identifiers were recovered by matching original and renamed ABF recordings using exact SHA-256 file identity together with the retained experiment-day code. Within each experimental group, recordings carrying the same recovered day code were assigned to the same animal and different day codes to different animals.

The frozen primary multi-sweep cohort contains all six SCA3 cells with resolved animal identity: four from `SCA3_AN01` and two from `SCA3_AN02`. In WT, eight of twelve primary cells have resolved identity: four from `WT_AN01` and four from `WT_AN02`; `WT_03`, `WT_04`, `WT_06`, and `WT_07` remain unresolved from the currently available archive.

Machine-readable assignments are in `data/animal_id_recovery/accepted_cohort.csv` and `data/animal_id_recovery/animal_groups.csv`.
