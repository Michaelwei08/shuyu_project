# Association Metadata Requirements

## Required For Any Clinical Association

- `sample_id`
- disease group
- sample source
- sequencing depth or usable read-pair count
- cohort or study identifier
- batch or sequencing center if multiple batches are combined

## Strongly Preferred

- EBV status by pathology, EBER, or prior annotation
- stage or disease severity
- treatment status
- survival or progression outcome
- transplant status for PTLD
- immunosuppression status for PTLD and immune-state interpretation

## Feasibility Labels

- `association_ready`: viral calls, disease labels, depth, and core covariates are present.
- `descriptive_only`: viral calls are present but clinical metadata are incomplete.
- `blocked`: data access or required metadata are missing.
