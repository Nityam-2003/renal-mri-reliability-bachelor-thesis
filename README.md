# Renal MRI Reliability: Bachelor Thesis Code

This repository contains the thesis-specific processing, registration, descriptive-analysis, and statistical-analysis scripts used for a Bachelor thesis evaluating the reliability of quantitative renal MRI biomarkers.

> **Code provenance:** This is not a from-scratch standalone software package. The workflow integrates and adapts existing RESPECT processing and co-registration components, the UKAT quantitative-mapping toolkit, and established reliability-analysis methods. In particular, the CoV implementation follows a published formulation previously applied to DWI data and was adapted here for the quantitative renal MRI data structure. This repository contains the thesis-specific adaptations and orchestration code only; upstream and restricted code is not redistributed.

The analysis focuses on respiratory-triggered T1 mapping, MOLLI T1 mapping, and respiratory-triggered T2 mapping in healthy volunteers and patients with chronic kidney disease (CKD). The final statistical analyses treat these two study populations independently.

## Scope

The repository documents the code used to:

1. process quantitative renal MRI data and generate T1, T2, MOLLI, B0, and B1 outputs;
2. register quantitative maps to an anatomical reference image;
3. apply whole-kidney masks and generate the final reliability-analysis input tables;
4. generate descriptive whole-kidney biomarker-value tables from all available processed examinations;
5. generate distribution plots;
6. calculate coefficient of variation (CoV), intraclass correlation coefficient (ICC), and Bland--Altman agreement measures for repeatability and reproducibility analyses.

This is a thesis-specific workflow. It contains study-specific sequence naming rules and folder assumptions, and is not intended to be a general-purpose clinical processing pipeline.

## Repository contents

| File | Purpose |
|---|---|
| `functions2.py` | Adapted quantitative-map processing workflow for respiratory-triggered T1, T2 StimFit, MOLLI, B0, B1, B1 correction, and quality-control outputs. |
| `Rigid_registration.py` | Thesis-specific rigid-registration workflow for T1, T2, and MOLLI maps. |
| `Analysis.py` | Applies kidney masks to registered maps and regenerates the four final group-specific repeatability and reproducibility input tables. |
| `Descriptive_roi_values.py` | Generates descriptive whole-kidney biomarker tables from every available masked ROI map, independently of visit pairing and the reliability-analysis input tables. |
| `stat_T1T2_by_group.R` | Final reliability analysis for healthy volunteers and CKD patients: CoV, ICC, and Bland--Altman analysis. |
| `distribution_plots_separate.R` | Final distribution plots for T1, MOLLI, and T2 by imaging centre and study population. |

## Processing and analysis workflow

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
Kidney ROI generation / measurements
        |
        +-------------------------------+
        |                               |
        v                               v
Reliability input tables           Descriptive ROI tables
(Analysis.py)                      (Descriptive_roi_values.py)
        |                               |
        v                               +-- all available processed examinations
CoV, ICC, Bland--Altman                 +-- centre-specific summaries
(stat_T1T2_by_group.R)                  +-- all-included-centres summaries
        |
        v
Distribution and reliability figures
(distribution_plots_separate.R and R outputs)
```

The descriptive tables and the repeatability/reproducibility analyses serve different purposes. The descriptive workflow retains all available processed examinations, including repeated visits, whereas the reliability workflow uses the study-specific visit-pair definitions required for repeatability and reproducibility analysis.

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

### Thesis-specific adaptations and integrations

The following scripts contain thesis-specific integration, adaptation, and analysis logic. They should not be interpreted as independent from-scratch implementations of the complete underlying processing, registration, or statistical methods.

- `functions2.py`: adaptation of the processing workflow for the renal MRI sequences, including the 3-parameter respiratory-triggered T1 fit, MOLLI handling, T2 StimFit workflow, B1 correction, quality-control outputs, and study-specific sequence selection.
- `Rigid_registration.py`: study-specific rigid registration of the final quantitative maps to the anatomical reference image, including centre-specific handling of available acquisition layouts.
- `Analysis.py`: study-specific application of kidney masks and generation of the four final healthy-volunteer and CKD-patient repeatability and reproducibility input tables.
- `Descriptive_roi_values.py`: thesis-specific descriptive analysis of the already masked whole-kidney ROI maps. It scans every available processed examination and does not use the repeatability/reproducibility visit-pair definitions. Repeated visits are intentionally retained because these tables describe all available processed measurements rather than one value per participant.
- `stat_T1T2_by_group.R` and `distribution_plots_separate.R`: thesis-specific R implementations of established distribution and reliability analyses for the renal T1, MOLLI, and T2 data structure, centre-specific visit conventions, and separate healthy-volunteer and CKD-patient analyses. The CoV implementation follows the published formulation described by de Boer et al., which was previously applied to DWI data and adapted here for this study.

### Descriptive whole-kidney analysis

`Descriptive_roi_values.py` is independent of the paired reliability analyses. It reads the existing masked whole-kidney ROI maps from `Processed - Results` and calculates one examination-level mean for each available parameter after applying the same validity ranges used when rebuilding the reliability-analysis tables:

```text
T1:     500-3000 ms
MOLLI:  500-3000 ms
T2:     15-150 ms
```

Only finite voxels within these ranges are included. Because the ROI maps are already masked, zero-valued background voxels are excluded by these ranges.

The script classifies IDs beginning with 5 as healthy volunteers and all other IDs as CKD patients. Centres 01 and 03 are assigned to GE, while Centres 02 and 04 are assigned to Siemens.

For each population, the script writes two private CSV files:

```text
Descriptive_HealthyVolunteers_master.csv
Descriptive_CKDPatients_master.csv
Descriptive_HealthyVolunteers_centre_summary.csv
Descriptive_CKDPatients_centre_summary.csv
```

The master tables retain every available examination and include the participant ID, population, centre, vendor, raw visit label, parameter, examination-level mean ROI value, valid-voxel count, and inclusion status.

The centre-summary tables report, for each centre and parameter, the number of available measurements, number of unique participants, mean, standard deviation, median, minimum, and maximum. An additional `All included centres` row is generated for each parameter.

Repeated examinations are intentionally retained. These descriptive tables therefore summarize all available processed measurements and should not be interpreted as population reference values based on a single examination per participant.

### R analysis workflow

As stated in the thesis, statistical analyses were performed in R using RStudio. The scripts implement distribution analysis, CoV, ICC, and Bland--Altman analysis for the final quantitative renal MRI input tables. The CoV calculation follows the formulation described by de Boer et al. for repeated-measurement reliability analysis, which was applied previously to DWI data and adapted for this study. ICC and Bland--Altman analyses use established statistical methods implemented with standard R packages.

## Software requirements

### Python

The processing, registration, ROI-table, and descriptive-analysis scripts require Python and the following packages:

```text
numpy
nibabel
pydicom
matplotlib
SimpleITK
itk
pandas
```

They also require authorised access to the RESPECT Processing and Co-Registration Modules and access to UKAT where applicable.

`Descriptive_roi_values.py` specifically requires:

```text
numpy
pandas
SimpleITK
```

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
+-- Processed - Results/           # Generated subject-level outputs and masked ROI maps
+-- Results/Separate/              # Final aggregate figures and tables
+-- code or notebooks/
    +-- Automate/
    |   +-- functions2.py
    +-- Registration/
    |   +-- Rigid_registration.py
    +-- Analysis/Separate/
        +-- Analysis.py
        +-- Descriptive_roi_values.py
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

### Reliability-analysis input tables

Run `Analysis.py` from `Analysis/Separate/`. It applies the available kidney masks to the registered maps and regenerates the four private input tables used by the R scripts:

```bash
python Analysis.py
```

The four non-public long-format input CSV files contain repeatability and reproducibility measurements for healthy volunteers and CKD patients. Each row represents a whole-kidney ROI measurement together with the study ID, imaging centre, vendor, visit, parameter, and measurement value.

### Descriptive whole-kidney tables

After the masked ROI maps have been generated, run:

```bash
python Descriptive_roi_values.py
```

The script scans every available examination folder in `Processed - Results`. It does not use the four reliability-analysis input tables and does not apply repeatability or reproducibility visit-pair selection.

It generates the following private files in `Analysis/Separate/`:

```text
Descriptive_HealthyVolunteers_master.csv
Descriptive_CKDPatients_master.csv
Descriptive_HealthyVolunteers_centre_summary.csv
Descriptive_CKDPatients_centre_summary.csv
```

The descriptive outputs retain repeated visits intentionally and are used to summarize the available whole-kidney T1, MOLLI, and T2 measurements by centre and across all included centres.

### Statistical analysis

Run the R scripts from `Analysis/Separate/` after generating the reliability input tables:

```r
source("stat_T1T2_by_group.R")
source("distribution_plots_separate.R")
```

The R scripts operate on the paired repeatability and reproducibility input tables. The descriptive CSV files generated by `Descriptive_roi_values.py` are separate outputs and are not used as replacements for those paired reliability tables.

## Data availability and privacy

Raw DICOM data, segmentation masks, participant identifiers, subject-level quantitative maps, masked ROI maps, reliability-analysis input CSV tables, descriptive master tables, descriptive centre-summary tables, and generated participant-level outputs are not included. These materials are restricted because they are governed by study-data access requirements and participant privacy.

The upstream RESPECT modules are also not redistributed here. Researchers seeking to reproduce the complete workflow require authorised access to the appropriate data and upstream modules.

`Descriptive_roi_values.py` reads restricted thesis data. If the local project policy requires this script itself to remain private, it should not be redistributed in a public repository; the description above documents its role in the thesis workflow without implying that restricted data are included.

## Citation and acknowledgement

Please acknowledge the RESPECT Processing Module, RESPECT Co-Registration Module, UKAT, ITK/Elastix, the published de Boer et al. CoV formulation, and the R packages used when reusing or adapting this workflow. Consult the relevant upstream repositories and the thesis references for citation and licence information.
