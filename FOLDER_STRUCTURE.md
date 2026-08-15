# Project folder structure

This document describes the local project layout used for the renal MRI Bachelor thesis. It is intentionally generic: it contains no participant identifiers, absolute paths, acquisition files, or study data.

```text
Everything related to data and analysis/
+-- code or notebooks/
|   +-- Automate/                 # Adapted quantitative-map processing workflow
|   +-- Registration/             # Thesis-specific rigid registration workflow
|   +-- Analysis/
|   |   +-- Separate/             # Final analysis: healthy volunteers and CKD patients
|   |   +-- Together/             # Earlier pooled analysis; not used for final results
|   +-- T1/ and T2/               # Earlier development and exploration code
|   +-- __pycache__/              # Local Python cache; do not publish
|   +-- .ipynb_checkpoints/       # Local notebook cache; do not publish
+-- Original data/                # Restricted source DICOM data; do not publish
+-- Masks/                        # Restricted kidney segmentation masks; do not publish
+-- Processed - Results/          # Subject-level generated maps, QC, registration, and ROI outputs; do not publish
+-- Results/
    +-- Separate/                 # Final aggregate tables and figures used in the thesis
```

## Processing flow

1. `Original data/` provides the DICOM input series.
2. `Automate/` creates quantitative T1, T2, MOLLI, B0, and B1 outputs in `Processed - Results/`.
3. `Registration/` registers the quantitative maps to the anatomical reference image.
4. Segmentation masks from `Masks/` are used to derive kidney ROI measurements.
5. `Analysis/Separate/` uses the final healthy-volunteer and CKD input tables to generate distributions, CoV, ICC, and Bland--Altman outputs.
6. Final thesis figures and tables are collected under `Results/Separate/`.

## Public-repository boundary

Only final thesis-specific scripts and documentation should be shared publicly. Raw data, masks, generated subject-level outputs, participant identifiers, cached files, preliminary notebooks, and previous pooled-analysis outputs remain private.
