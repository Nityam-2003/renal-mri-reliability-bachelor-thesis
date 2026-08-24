# Renal MRI Reliability: Bachelor Thesis Code

This repository contains the thesis-specific processing, registration, ROI-extraction, descriptive-analysis, and statistical-analysis code used to evaluate quantitative renal MRI biomarkers in healthy volunteers and patients with chronic kidney disease (CKD).

The analysis covers respiratory-triggered T1 mapping, MOLLI T1 mapping, and respiratory-triggered T2 mapping. Repeatability and reproducibility were evaluated separately for healthy volunteers and CKD patients using coefficient of variation (CoV), intraclass correlation coefficient (ICC), and Bland--Altman analysis.

> **Scope and provenance:** This is a study-specific research workflow, not a standalone clinical application or a from-scratch implementation of every underlying method. It adapts and connects components from the RESPECT processing and co-registration workflows, UKAT quantitative mapping, 3D Slicer BRAINSFit, and established statistical methods. Restricted upstream modules, study data, and participant-level outputs are not distributed here.

## Final workflow

```text
Restricted DICOM examinations
        |
        v
Quantitative-map generation and QC
(Automate/functions2.py)
        |
        +-- respiratory-triggered T1: three-parameter fit + B1 correction
        +-- MOLLI T1: three-parameter fit + B1 correction
        +-- respiratory-triggered T2: UKAT StimFit
        +-- B0 and B1 maps
        |
        v
Prepared anatomical T1-weighted fixed images
        |
        v
Final map registration with 3D Slicer BRAINSFit
(Registration/brainsfit_batch_registration.py)
        |
        v
Whole-kidney mask application and ROI measurement
(Analysis/BRAINSFit_Registration/Analysis.py)
        |
        +-------------------------------+
        |                               |
        v                               v
Reliability input tables           Descriptive master tables
(paired visit definitions)         (all available examinations)
        |                               |
        v                               v
CoV, ICC, Bland--Altman            Centre-wise distribution plots
(stat_T1T2_by_group.R)             (distribution_plots_separate.R)
```

An earlier ITK/Elastix registration branch is retained in `Registration/Rigid_registration.py` because it formed part of the thesis workflow and technical evaluation. The final reported ROI measurements were generated from the BRAINSFit-registered maps. The repository is therefore not intended as a general comparison of registration software.

## Main repository contents

| Path | Purpose |
|---|---|
| `Automate/functions2.py` | Main quantitative-processing workflow for T1, MOLLI T1, T2 StimFit, B0, B1, B1 correction, and quality-control outputs. |
| `Registration/Rigid_registration.py` | Initial thesis-specific ITK/Elastix registration workflow, retained for provenance and comparison. |
| `Registration/brainsfit_batch_registration.py` | Final batch registration of unmasked quantitative maps directly to prepared anatomical T1-weighted fixed images using 3D Slicer BRAINSFit. |
| `Registration/run_brainsfit_registration.bat` | Windows launcher that runs the BRAINSFit batch script through 3D Slicer. |
| `Registration/Registration.ipynb` | Working notebook used during registration development and preparation or verification of fixed images. It is not the preferred batch entry point. |
| `Analysis/BRAINSFit_Registration/Analysis.py` | Applies kidney masks to the final BRAINSFit maps, saves distinct BRAINSFit ROI files, and builds the four repeatability/reproducibility input tables. |
| `Analysis/BRAINSFit_Registration/Descriptive_roi_values.py` | Builds examination-level and centre-summary descriptive tables from all available BRAINSFit ROI maps. |
| `Analysis/BRAINSFit_Registration/stat_T1T2_by_group.R` | Final CoV, ICC, and Bland--Altman analyses for healthy volunteers and CKD patients. |
| `Analysis/BRAINSFit_Registration/distribution_plots_separate.R` | Final centre-wise distribution plots based on the descriptive master tables. |
| `Analysis/Separate/` | Earlier Elastix-based analysis branch, retained for provenance. |
| `Analysis/T1_2_parameter/` | Superseded exploratory two-parameter T1 analysis; not used for the final thesis results. |
| `Analysis/Together/` | Earlier combined-population exploratory analysis; not used for the final separated healthy-volunteer and CKD results. |

Other notebooks and scripts record development, testing, data inspection, or earlier processing attempts. The files listed above are the principal entry points for reproducing the final workflow.

## Quantitative processing

`Automate/functions2.py` identifies the available sequences using study-specific series names and DICOM metadata. Missing or unprocessable sequences are skipped and recorded in the processing log.

The final processing choices were:

- **Respiratory-triggered T1:** three-parameter inversion-recovery fitting, followed by B1 correction.
- **MOLLI T1:** three-parameter fitting, followed by B1 correction.
- **Respiratory-triggered T2:** UKAT StimFit with vendor-specific models. StimFit accounts for stimulated-echo and B1-related refocusing effects, so no separate post-fit B1 correction was applied to T2.
- **B0:** maps were generated and stored but were not used in the final reliability analysis.
- **Quality control:** parameter maps, fit-quality outputs, numerical summaries, and processing logs were generated where applicable.

The code contains study-specific sequence names, centre conventions, and folder assumptions. These must be adapted for another dataset.

## Registration

### Initial ITK/Elastix branch

`Registration/Rigid_registration.py` uses ITK/Elastix to estimate rigid transformations from sequence-specific source images and applies them to the corresponding quantitative maps. It depends on helper functionality associated with the RESPECT Co-Registration Module. This branch is retained because it was used during the study and provided the initial registered outputs.

Run it from the `Registration` directory:

```bash
python Rigid_registration.py
```

### Final 3D Slicer BRAINSFit branch

The final branch registers each unmasked quantitative map directly to an already prepared T1-weighted fixed image for the same examination. The batch script uses the settings applied through 3D Slicer's **General Registration (BRAINS)** module:

- geometry-based initialisation;
- sampling percentage: `0.002`;
- rigid stage enabled;
- global-scale stage enabled;
- affine and B-spline stages disabled;
- linear interpolation;
- background fill value: `0`.

This is a linear rigid-plus-scale registration, not deformable registration.

Expected inputs within each examination folder are:

```text
Registration/T1/fixed.nii.gz
Registration/T2/fixed.nii.gz
Registration/MOLLI/fixed.nii.gz

T1_RespTrig/nifti/t1_map_b1corr.nii.gz
T2_RespTrig/nifti/stimfit_t2_map.nii.gz
MOLLI/nifti/molli_t1_map_b1corr.nii.gz
```

The registered maps and transforms are written without overwriting the Elastix outputs:

```text
Registration_BRAINS/T1/t1_map_registered_brains.nii.gz
Registration_BRAINS/T1/brainsfit_transform.h5

Registration_BRAINS/T2/stimfit_t2_map_registered_brains.nii.gz
Registration_BRAINS/T2/brainsfit_transform.h5

Registration_BRAINS/MOLLI/molli_t1_map_registered_brains.nii.gz
Registration_BRAINS/MOLLI/brainsfit_transform.h5
```

Before running the batch, edit these local settings if required:

1. `PROCESSED_RESULTS` in `brainsfit_batch_registration.py`;
2. `SLICER_EXE` in `run_brainsfit_registration.bat`.

Then run the Windows launcher:

```bat
run_brainsfit_registration.bat
```

The Python file must be executed by 3D Slicer rather than a standard Python or Jupyter environment. Existing outputs are protected by default, missing inputs are skipped, and each parameter/examination result is appended to:

```text
Processed - Results/BRAINSFit_batch_registration_log.csv
```

Software completion does not guarantee anatomically valid alignment. Registered maps and transferred kidney masks should therefore be visually checked across all retained slices, especially the outer slices.

## ROI extraction and analysis tables

Run the final analysis branch from `Analysis/BRAINSFit_Registration`:

```bash
python Analysis.py
```

The script:

1. reads the BRAINSFit-registered T1, MOLLI T1, and T2 maps;
2. locates the corresponding whole-kidney mask;
3. resamples the mask to the registered-map geometry using nearest-neighbour interpolation;
4. saves the masked maps as `t1_roi_brainsfit.nii.gz`, `molli_roi_brainsfit.nii.gz`, and `t2_roi_brainsfit.nii.gz` inside each examination's `Analysis` folder;
5. calculates the whole-kidney mean from finite voxels within the final validity ranges;
6. writes four long-format input tables for the separate population and reliability analyses.

Final validity ranges:

```text
T1:        500-2500 ms
MOLLI T1:  500-2500 ms
T2:         15-150 ms
```

The four private reliability input tables are:

```text
Repeatability_HealthyVolunteers_input.csv
Reproducibility_HealthyVolunteers_input.csv
Repeatability_Patients_input.csv
Reproducibility_Patients_input.csv
```

Centre-specific raw visit labels are mapped to the common `v2a` and `v3` columns expected by the R analysis. Rows from an available single visit may remain in an input table, but the R script forms complete pairs and excludes unmatched visits before calculating reliability statistics.

## Descriptive whole-kidney tables

After the BRAINSFit ROI maps have been generated, run:

```bash
python Descriptive_roi_values.py
```

Unlike the reliability analysis, this script does not apply visit-pair definitions. It retains every available processed examination and writes:

```text
Descriptive_HealthyVolunteers_master.csv
Descriptive_CKDPatients_master.csv
Descriptive_HealthyVolunteers_centre_summary.csv
Descriptive_CKDPatients_centre_summary.csv
```

The master tables contain participant, population, centre, vendor, raw visit, parameter, examination-level mean, valid-voxel count, and inclusion status. The summary tables report the number of measurements, number of unique participants, mean, standard deviation, median, minimum, and maximum for each centre and parameter, together with an `All included centres` row.

Repeated examinations are intentionally retained. The resulting means describe all available processed measurements and are not single-visit population reference values.

## Statistical analysis and distribution plots

Run the R scripts after generating their respective input tables:

```r
source("stat_T1T2_by_group.R")
source("distribution_plots_separate.R")
```

`stat_T1T2_by_group.R` reads the four repeatability/reproducibility tables, constructs complete visit pairs, and calculates:

- overall, centre-specific, and vendor-specific CoV;
- overall, centre-specific, and vendor-specific ICC with confidence intervals;
- Bland--Altman mean relative difference and 95% limits of agreement;
- corresponding tables and figures for healthy volunteers and CKD patients separately.

The CoV implementation follows the repeated-measurement formulation described by de Boer et al. and was adapted to the final renal MRI data structure.

`distribution_plots_separate.R` does **not** read or combine the four paired-analysis tables. It reads the two descriptive master tables so that each available processed examination is included once in the appropriate population-specific distribution dataset.

## Software requirements

### Python processing and Elastix registration

The local workflow used Python with the following principal packages:

```text
numpy
pandas
nibabel
pydicom
matplotlib
SimpleITK
itk (with Elastix support)
UKAT
```

The processing and initial registration branches additionally require authorised access to the relevant RESPECT helper modules.

### BRAINSFit registration

- 3D Slicer with the BRAINSFit/General Registration (BRAINS) command-line module;
- Windows for the supplied `.bat` launcher, or an equivalent platform-specific command invoking Slicer with `--python-script`.

The Slicer Python environment supplies the `slicer` module; installing the unrelated PyPI package named `slicer` is not a substitute.

### R

The final R scripts require:

```r
install.packages(c("tidyverse", "irr"))
```

## Local folder layout

The scripts assume the following study-specific layout:

```text
Everything related to data and analysis/
+-- Original data/                         # restricted DICOM examinations
+-- Masks/                                 # restricted kidney masks
+-- Processed - Results/                   # generated examination-level outputs
+-- code or notebooks/
    +-- Automate/
    |   +-- functions2.py
    +-- Registration/
    |   +-- Rigid_registration.py
    |   +-- brainsfit_batch_registration.py
    |   +-- run_brainsfit_registration.bat
    |   +-- Registration.ipynb
    +-- Analysis/
        +-- BRAINSFit_Registration/         # final analysis branch
        |   +-- Analysis.py
        |   +-- Descriptive_roi_values.py
        |   +-- stat_T1T2_by_group.R
        |   +-- distribution_plots_separate.R
        +-- Separate/                       # earlier Elastix-based branch
        +-- T1_2_parameter/                 # superseded exploratory branch
        +-- Together/                       # earlier pooled exploratory branch
```

The code was written for this layout. Absolute paths in the BRAINSFit batch files and relative paths in earlier scripts must be reviewed before use on another computer.

## Example processing call

Run from the `Automate` folder in a Python or Jupyter environment after configuring the required dependencies and data access:

```python
from functions2 import process_patient

patient_folder = "../../Original data/<examination_folder>"
process_patient(patient_folder)
```

The workflow skips already completed sequence outputs and records missing or failed sequences in `Processed - Results/Pipeline_Log.txt`.

## External projects and acknowledgements

### RESPECT Processing Module

The DICOM loading, organisation, and processing workflow was adapted from the RESPECT Processing Module. Access to the complete upstream module may be restricted.

- [RESPECT Processing Module](https://github.com/Computer-Assisted-Clinical-Medicine/RESPECT_Processing_Module)

### RESPECT Co-Registration Module

The initial ITK/Elastix branch uses helper functionality and conventions associated with the RESPECT Co-Registration Module. Restricted upstream components are not redistributed.

- [RESPECT Co-Registration Module](https://github.com/Computer-Assisted-Clinical-Medicine/RESPECT_Co-Registration_Module)

### UKAT

Quantitative mapping uses UKAT functionality, including the T1, T2 StimFit, B0, and utility components.

- [UKAT: UKRIN-MAPS](https://github.com/UKRIN-MAPS/ukat)

### 3D Slicer and BRAINSFit

The final quantitative-map registration was performed through 3D Slicer using BRAINSFit. Users should consult the 3D Slicer and BRAINSFit documentation for software citation and licence information.

- [3D Slicer](https://www.slicer.org/)
- [BRAINSFit documentation](https://www.slicer.org/wiki/Documentation/Nightly/Modules/BRAINSFit)

## Data availability and privacy

Raw DICOM data, kidney masks, participant identifiers, subject-level quantitative maps, registered maps, masked ROI maps, input CSV files, descriptive tables, and participant-level outputs are not included because they are subject to study-data access requirements and participant privacy restrictions.

Reproducing the complete workflow therefore requires authorised access to the study data and any restricted upstream modules. Generated CSV files and figures should be reviewed before publication because filenames and table contents may contain study identifiers.

## Citation

When reusing or adapting this workflow, cite the relevant RESPECT modules, UKAT, 3D Slicer/BRAINSFit, ITK/Elastix, the statistical methods, and the R packages used. Consult the thesis reference list and the upstream projects for the appropriate citations and licence terms.
