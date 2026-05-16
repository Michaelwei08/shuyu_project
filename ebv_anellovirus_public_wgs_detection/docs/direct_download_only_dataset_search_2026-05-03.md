# Direct-Download-Only Dataset Search

## Scope

This pass keeps only datasets that can be downloaded directly without EGA, dbGaP, GDC controlled-access tokens, UK Biobank approval, All of Us workspace access, or JGA/HGVD approval. I interpreted "directly downloadable" as open SRA/DDBJ/ENA-style runs with public download paths or SRA Toolkit access.

Controlled raw WGS datasets remain scientifically better for HL/DLBCL/PTLD, but they are out of scope for this direct-download-only pass.

## Best Direct Choices

1. `PRJNA809536`: EBV-positive classic Hodgkin lymphoma in XMEN twins. This is the best direct HL/EBV human WGS hit. It has 3 direct WGS runs, 2 above the 40 GB compressed-size heuristic, and 142.1 GB total direct WGS data.

2. `PRJNA172563`: CGCI NHL/DLBCL cell-line WGS public subset. This is the best direct DLBCL/NHL WGS-scale dataset. It has 13 direct WGS runs, 12 above 40 GB, and 627.5 GB total direct WGS data. Important caveat: most WGS runs in the project are `@dbgap@` and are excluded here.

3. `PRJNA772081`: inborn-errors-of-immunity lymphoma data. It has 9 direct human WGS runs and 221.5 GB direct WGS data. The individual WGS runs are below the 40 GB compressed-size heuristic, but the biology is relevant to immunodeficiency-associated lymphoma.

4. `PRJDB11215`: HTLV-1-infected cell-line WGS. This is directly useful for P1 HTLV discrimination controls. It has 3 direct WGS runs and 88.1 GB total direct WGS data.

5. EBV-positive lymphoma cell-line controls: `PRJNA1212064` EB-3, `PRJNA998843` BC-3, and `PRJNA528941` Raji. These are not disease cohorts, but they are directly downloadable herpesvirus-positive lymphoma WGS controls.

6. EBV-only resources: `PRJDB18873` and `PRJNA551423`. These should be used for EBV reference/variant validation only, not for host-WGS or anellovirus claims.

## Explicit Exclusions

- `CGCI-HTMCP-DLBCL`, `TCGA-DLBC`, `CGCI-BLGSP`: raw WGS is controlled in GDC/dbGaP. Excluded despite being the best scientific fits.
- EGA HL/DLBCL datasets: controlled access. Excluded.
- UK Biobank, All of Us, SPARK, JGA/HGVD: qualified or controlled access. Excluded.
- GDC open VCF/TXT/TSV/BEDPE files: directly downloadable, but not raw WGS reads and not suitable for EBV/anellovirus read detection.

## Download Commands

Use SRA Toolkit for reproducible download:

```powershell
prefetch SRR18110684
fasterq-dump SRR18110684 --split-files --threads 8
```

For direct SRA file URLs, use the `download_path` column in the downloaded RunInfo manifests under `data/external_manifests/`. I verified representative direct URLs with HTTP `HEAD` responses returning `200`.

## Files Added

- `data/cohorts/direct_download_only_wgs_candidates_2026-05-03.csv`
- `data/external_manifests/ncbi_open_sra_search_direct_download_candidates_2026-05-03.csv`
- `data/external_manifests/ncbi_runinfo_PRJNA809536_direct_download_2026-05-03.csv`
- `data/external_manifests/ncbi_runinfo_PRJNA172563_direct_download_2026-05-03.csv`
- `data/external_manifests/ncbi_runinfo_PRJNA772081_direct_download_2026-05-03.csv`
- `data/external_manifests/ncbi_runinfo_PRJDB11215_direct_download_2026-05-03.csv`
- Additional direct-control RunInfo files for EBV/HHV8 lymphoma cell-line and EBV-only datasets.
