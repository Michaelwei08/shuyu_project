# Direct-Download Smoke Test

## Bottom Line

The directly downloadable datasets support method smoke testing and viral-positive controls, but they do not fully support the original clinical-cohort hypothesis. The missing piece is still a directly downloadable large patient-level HL/DLBCL/PTLD WGS cohort with usable clinical metadata.

## Hypothesis Matrix

| Hypothesis | Status | Evidence |
|---|---:|---|
| Direct raw WGS exists for a target HL/DLBCL/PTLD disease context | partial | PRJNA809536 supports HL/EBV but is a 3-run case-report-scale dataset; PRJNA172563 supports DLBCL/NHL but mainly cell lines. |
| Direct raw WGS is large enough for detector smoke testing | pass | At least one retained dataset has direct WGS runs >=40 GB; PRJNA172563 has 12 and PRJNA809536 has 2. |
| Direct data can validate EBV/herpesvirus positive controls | pass | HL EBV context, Raji/EB-3 EBV-positive cell lines, BC-3 HHV8-positive PEL, and EBV viral genomes are directly downloadable. |
| Direct data can support anellovirus detection in lymphoma WGS | weak | No direct dataset has known anellovirus-positive lymphoma WGS annotation; anellovirus can only be explored opportunistically in host WGS. |
| Direct data can support clinical association replication | fail | Direct datasets are too small, cell-line-heavy, or viral-only; controlled datasets are still needed for clinical association replication. |
| Direct data can support HTLV/HIV/retroviral discrimination controls | partial | PRJDB11215 provides direct HTLV-1-infected cell-line WGS; no direct HIV cohort WGS with healthy controls was retained. |

## Top Dataset Verdicts

- `PRJNA809536`: supports_smoke_test_not_full_hypothesis; 3 direct WGS runs, 2 >=40 GB, 142.1 GB total. Primary use: `HL_EBV_positive_control`.
- `PRJNA172563`: partial_support; 13 direct WGS runs, 12 >=40 GB, 627.5 GB total. Primary use: `DLBCL_WGS_method_smoke_test`.
- `PRJNA772081`: partial_support; 9 direct WGS runs, 0 >=40 GB, 221.5 GB total. Primary use: `immunodeficiency_lymphoma_bridge`.
- `PRJNA998843`: partial_support; 2 direct WGS runs, 1 >=40 GB, 74.8 GB total. Primary use: `herpesvirus_mapping_control`.
- `PRJNA1212064`: partial_support; 1 direct WGS runs, 1 >=40 GB, 63.6 GB total. Primary use: `EBV_mapping_positive_control`.
- `PRJNA528941`: partial_support; 1 direct WGS runs, 1 >=40 GB, 51.9 GB total. Primary use: `EBV_mapping_positive_control`.
- `PRJDB18873`: usable_for_EBV_reference_only; 295 direct WGS runs, 0 >=40 GB, 143.7 GB total. Primary use: `EBV_reference_validation`.
- `PRJDB11215`: weak_support; 3 direct WGS runs, 0 >=40 GB, 88.1 GB total. Primary use: `HTLV_discrimination_control`.
- `PRJNA551423`: usable_for_PTLD_EBV_variant_only; 0 direct WGS runs, 0 >=40 GB, 0 GB total. Primary use: `PTLD_EBV_variant_reference`.
