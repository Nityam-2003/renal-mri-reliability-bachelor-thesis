# Renal MRI Reliability: Bachelor Thesis Code

This repository contains the thesis-specific processing, registration, and statistical-analysis scripts developed for a Bachelor thesis evaluating the reliability of quantitative renal MRI biomarkers.

The analysis focuses on respiratory-triggered T1 mapping, MOLLI T1 mapping, and respiratory-triggered T2 mapping in healthy volunteers and patients with chronic kidney disease (CKD). The final statistical analyses treat these two study populations independently.

## Scope

The repository documents the code used to:

1. process quantitative renal MRI data and generate T1, T2, MOLLI, B0, and B1 outputs;
2. register quantitative maps to an anatomical reference image;
3. generate distribution plots;
4. calculate coefficient of variation (CoV), intraclass correlation coefficient (ICC), and Bland--Altman agreement measures for repeatability and reproducibility analyses.

This is a thesis-specific workflow. It contains study-specific sequence naming rules and folder assumptions, and is not intended to be a general-purpose clinical processing pipeline.

## Repository contents

| File | Purpose |
|---|---|
| `functions2.py` | Adapted quantitative-map processing workflow for respiratory-triggered T1, T2 StimFit, MOLLI, B0, B1, B1 correction, and quality-control outputs. |
| `Rigid_registration.py` | Thesis-specific rigid-registration workflow for T1, T2, and MOLLI maps. |
| `stat_T1T2_by_group.R` | Final reliability analysis for healthy volunteers and CKD patients: CoV, ICC, and Bland--Altman analysis. |
| `distribution_plots_separate.R` | Final distribution plots for T1, MOLLI, and T2 by imaging centre and study population. |

## Processing and registration workflow

```text
Restricted DICOM data
        |
        v
Quantitative processing (functions2.py)
        |
        +-- T1, T2, MOLLI, B0, B1 maps and QC outputs
        |
        v
Rigid registration (Rigid_registration.py)
        |
        v
Kidney ROI measurements and final analysis input tables
        |
        v
Statistical analysis and figures (R scripts)
```

## External code and software provenance

This repository intentionally does **not** redistribute restricted upstream modules or third-party source code.

### RESPECT Processing Module

The quantitative processing workflow was adapted for this thesis from the RESPECT Processing Module. The upstream module provides helper functionality for loading and organising DICOM data. Access to the full module may be restricted.

- [RESPECT Processing Module](https://github.com/Computer-Assisted-Clinical-Medicine/RESPECT_Processing_Module)

### RESPECT Co-Registration Module

The rigid-registration workflow uses ITK/Elastix and relies on helper functionality associated with the RESPECT Co-Registration Module. The module itself, including its helper scripts and parameter maps, is not redistributed here. Access to the full module may be restricted.

- [RESPECT Co-Registration Module](https://github.com/Computer-Assisted-Clinical-Medicine/RESPECT_Co-Registration_Module)

### UKAT

Quantitative mapping uses the UKAT toolkit, including T1 mapping, T2 StimFit, B0 mapping, and utility functions. The UKAT source code is not redistributed in this repository.

- [UKAT: UKRIN-MAPS](https://github.com/UKRIN-MAPS/ukat)

### Thesis-developed code

The following scripts were developed and adapted for the requirements of this thesis:

- `functions2.py`: adaptation of the processing workflow for the renal MRI sequences, including the 3-parameter respiratory-triggered T1 fit, MOLLI handling, T2 StimFit workflow, B1 correction, quality-control outputs, and study-specific sequence selection.
- `Rigid_registration.py`: study-specific rigid registration of the final quantitative maps to the anatomical reference image, including centre-specific handling of available acquisition layouts.
- `stat_T1T2_by_group.R` and `distribution_plots_separate.R`: thesis-developed analysis scripts for the final healthy-volunteer and CKD-patient analyses. These scripts use standard R packages; they do not redistribute external analysis scripts.

## Software requirements

### Python

The processing and registration scripts require Python and the following packages:

```text
numpy
nibabel
pydicom
matplotlib
SimpleITK
itk
```

They also require authorised access to the RESPECT Processing and Co-Registration Modules and access to UKAT.

### R

The statistical-analysis scripts require R and the following packages:

```r
install.packages(c("tidyverse", "irr"))
```

## Local folder layout

The scripts were run within the following local layout:

```text
Everything related to data and analysis/
+-- Original data/                 # Restricted DICOM input data
+-- Masks/                         # Restricted kidney segmentation masks
+-- Processed - Results/           # Generated subject-level outputs
+-- Results/Separate/              # Final aggregate figures and tables
+-- code or notebooks/
    +-- Automate/
    |   +-- functions2.py
    +-- Registration/
    |   +-- Rigid_registration.py
    +-- Analysis/Separate/
        +-- stat_T1T2_by_group.R
        +-- distribution_plots_separate.R
```

The relative paths in the Python scripts reflect this project layout. Users with a different layout must adapt the paths or provide an equivalent configuration.

## Use

### Quantitative-map processing

Run from the `Automate` folder in a Python notebook or interactive Python session:

```python
from functions2 import process_patient

patient_folder = "../../Original data/<patient_folder_name>"
process_patient(patient_folder)
```

To process all available patient folders:

```python
import os
from functions2 import process_patient

original_data = "../../Original data"

for patient in sorted(os.listdir(original_data)):
    patient_folder = os.path.join(original_data, patient)
    if os.path.isdir(patient_folder):
        process_patient(patient_folder)
```

### Rigid registration

Run from the `Registration` folder:

```bash
python Rigid_registration.py
```

The script checks each available subject folder and registers each available T1, T2, and MOLLI map independently. Existing registered maps are skipped.

### Statistical analysis

Run the R scripts from the directory containing the four final, private input tables:

```r
source("stat_T1T2_by_group.R")
source("distribution_plots_separate.R")
```

The scripts expect four non-public long-format input CSV files: repeatability and reproducibility input tables for healthy volunteers and CKD patients. Each row represents a whole-kidney ROI measurement together with the study ID, imaging centre, vendor, visit, parameter, and measurement value.

## Data availability and privacy

Raw DICOM data, segmentation masks, participant identifiers, subject-level quantitative maps, input CSV tables, and generated participant-level outputs are not included. These materials are restricted because they are governed by study-data access requirements and participant privacy.

The upstream RESPECT modules are also not redistributed here. Researchers seeking to reproduce the complete workflow require authorised access to the appropriate data and upstream modules.

## Citation and acknowledgement

Please acknowledge the RESPECT Processing Module, RESPECT Co-Registration Module, UKAT, ITK/Elastix, and the R packages used when reusing or adapting this workflow. Consult the relevant upstream repositories for their citation and licence information.
