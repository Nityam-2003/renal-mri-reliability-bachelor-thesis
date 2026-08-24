# Renal MRI Reliability: Bachelor Thesis Code

This repository contains the thesis-specific Python and R scripts used to evaluate the reliability of quantitative renal MRI biomarkers in healthy volunteers and patients with chronic kidney disease (CKD).

The analysis covers respiratory-triggered T1 mapping, MOLLI T1 mapping, and respiratory-triggered T2 mapping. Healthy volunteers and CKD patients were analysed separately. Repeatability and reproducibility were assessed using the coefficient of variation (CoV), intraclass correlation coefficient (ICC), and Bland--Altman analysis.

> **Code provenance:** This is a study-specific research workflow rather than a standalone clinical application. It adapts and connects functionality from the RESPECT processing and co-registration workflows, UKAT, 3D Slicer BRAINSFit, ITK/Elastix, and established statistical methods. Restricted upstream code and study data are not redistributed.

## Included files

| File | Purpose |
|---|---|
| `functions2.py` | Quantitative-map processing for respiratory-triggered T1, MOLLI T1, T2 StimFit, B0, B1, B1 correction, and quality-control outputs. |
| `Rigid_registration.py` | Initial ITK/Elastix registration workflow used during the thesis. |
| `brainsfit_batch_registration.py` | Final batch registration of quantitative maps to the prepared anatomical T1-weighted reference images using 3D Slicer BRAINSFit. |
| `Analysis.py` | Applies whole-kidney masks to the BRAINSFit-registered maps and generates the four repeatability and reproducibility input tables. |
| `Descriptive_roi_values.py` | Generates examination-level and centre-summary descriptive whole-kidney biomarker tables from all available processed measurements. |
| `stat_T1T2_by_group.R` | Calculates CoV, ICC, and Bland--Altman statistics for healthy volunteers and CKD patients. |
| `distribution_plots_separate.R` | Generates centre-wise T1, MOLLI T1, and T2 distribution plots separately for healthy volunteers and CKD patients. |

Only the thesis-relevant scripts are included. Study data, generated results, temporary development files, and upstream helper modules are excluded.

## Workflow

```text
Restricted DICOM examinations
        |
        v
Quantitative-map generation and quality control
(functions2.py)
        |
        +-- respiratory-triggered T1: three-parameter fit and B1 correction
        +-- MOLLI T1: three-parameter fit and B1 correction
        +-- respiratory-triggered T2: UKAT StimFit
        +-- B0 and B1 maps
        |
        v
Registration to the anatomical T1-weighted reference
        |
        +-- initial ITK/Elastix workflow (Rigid_registration.py)
        +-- final BRAINSFit workflow (brainsfit_batch_registration.py)
        |
        v
Whole-kidney ROI extraction and reliability input tables
(Analysis.py)
        |
        +-------------------------------+
        |                               |
        v                               v
Reliability analysis                Descriptive analysis
(stat_T1T2_by_group.R)              (Descriptive_roi_values.py)
                                        |
                                        v
                                   Distribution plots
                                   (distribution_plots_separate.R)
```

The initial ITK/Elastix registration is retained because it formed part of the thesis workflow. The final reported ROI measurements were obtained from the BRAINSFit-registered maps. Registration-method comparison was a technical quality-control step rather than the main purpose of the study.

## Quantitative processing

`functions2.py` identifies the available sequences using study-specific series names and DICOM metadata. Missing or unprocessable sequences are skipped and recorded in the processing log.

The final processing choices were:

- **Respiratory-triggered T1:** three-parameter inversion-recovery fitting followed by B1 correction.
- **MOLLI T1:** three-parameter fitting followed by B1 correction.
- **Respiratory-triggered T2:** UKAT StimFit using vendor-specific models. StimFit accounts for stimulated-echo and B1-related refocusing effects, so no separate post-fit B1 correction was applied to T2.
- **B0:** maps were generated and stored but were not used in the final reliability analysis.
- **Quality control:** parameter maps, fit-quality outputs, numerical summaries, and processing logs were generated where applicable.

The script depends on authorised access to the required RESPECT processing helpers and UKAT functionality. These upstream modules are not included in this repository.

## Registration

### ITK/Elastix

`Rigid_registration.py` contains the initial thesis-specific rigid-registration workflow. It estimates transformations from the sequence-specific source images and applies them to the corresponding quantitative maps. The script relies on helper functionality associated with the RESPECT Co-Registration Module, which is not redistributed here.

### 3D Slicer BRAINSFit

`brainsfit_batch_registration.py` contains the final registration workflow. It must be executed through 3D Slicer rather than a standard Python installation. Each unmasked quantitative map is registered directly to the prepared T1-weighted fixed image for the same examination.

The BRAINSFit settings were:

- geometry-based initialisation;
- sampling percentage of `0.002`;
- rigid stage enabled;
- global-scale stage enabled;
- affine and B-spline stages disabled;
- linear interpolation;
- background fill value of `0`.

This is a linear rigid-plus-scale registration, not deformable registration. Existing registered outputs are protected by default, missing inputs are skipped, and the result of every parameter/examination attempt is written to `BRAINSFit_batch_registration_log.csv`.

Expected inputs within each processed examination are:

```text
Registration/T1/fixed.nii.gz
Registration/T2/fixed.nii.gz
Registration/MOLLI/fixed.nii.gz

T1_RespTrig/nifti/t1_map_b1corr.nii.gz
T2_RespTrig/nifti/stimfit_t2_map.nii.gz
MOLLI/nifti/molli_t1_map_b1corr.nii.gz
```

The final outputs are written separately from the Elastix results:

```text
Registration_BRAINS/T1/t1_map_registered_brains.nii.gz
Registration_BRAINS/T2/stimfit_t2_map_registered_brains.nii.gz
Registration_BRAINS/MOLLI/molli_t1_map_registered_brains.nii.gz
```

Software completion alone does not establish anatomical validity. The registered maps and transferred masks should therefore be visually checked across all retained slices, especially the outer slices.

## ROI extraction and analysis tables

`Analysis.py` reads the BRAINSFit-registered maps, resamples each whole-kidney mask to the corresponding map geometry using nearest-neighbour interpolation, and saves the masked maps as:

```text
t1_roi_brainsfit.nii.gz
molli_roi_brainsfit.nii.gz
t2_roi_brainsfit.nii.gz
```

Whole-kidney means are calculated from finite voxels within the final validity ranges:

```text
T1:        500-2500 ms
MOLLI T1:  500-2500 ms
T2:         15-150 ms
```

The script generates four private long-format tables:

```text
Repeatability_HealthyVolunteers_input.csv
Reproducibility_HealthyVolunteers_input.csv
Repeatability_Patients_input.csv
Reproducibility_Patients_input.csv
```

Centre-specific visit labels are mapped to the common visit columns required by the R analysis. Measurements from a single available visit may remain in an input table, but `stat_T1T2_by_group.R` constructs complete visit pairs and excludes unmatched visits before calculating reliability statistics.

## Descriptive analysis and distribution plots

`Descriptive_roi_values.py` is independent of the paired reliability analysis. It reads all available BRAINSFit ROI maps and retains repeated examinations intentionally. It generates:

```text
Descriptive_HealthyVolunteers_master.csv
Descriptive_CKDPatients_master.csv
Descriptive_HealthyVolunteers_centre_summary.csv
Descriptive_CKDPatients_centre_summary.csv
```

The master tables contain the participant, population, centre, vendor, raw visit, parameter, examination-level mean, valid-voxel count, and inclusion status. The summary tables report the number of measurements, number of unique participants, mean, standard deviation, median, minimum, and maximum for each centre and parameter, together with an `All included centres` row.

These values describe all available processed measurements and are not single-visit population reference values. `distribution_plots_separate.R` reads the two descriptive master tables rather than the four paired reliability tables, preventing the repeatability and reproducibility inputs from duplicating measurements in the distributions.

## Statistical analysis

`stat_T1T2_by_group.R` constructs complete measurement pairs and calculates:

- overall, centre-specific, and vendor-specific CoV;
- overall, centre-specific, and vendor-specific ICC with confidence intervals;
- Bland--Altman mean relative differences and 95% limits of agreement;
- corresponding tables and figures for healthy volunteers and CKD patients separately.

The CoV implementation follows the repeated-measurement formulation described by de Boer et al. and was adapted to the renal MRI data structure used in this thesis.

## Software requirements

### Python

The Python scripts use the following principal packages and software:

```text
numpy
pandas
nibabel
pydicom
matplotlib
SimpleITK
ITK with Elastix support
UKAT
3D Slicer with BRAINSFit
```

`functions2.py` and `Rigid_registration.py` also rely on authorised upstream RESPECT helper functionality that is not included here. The `slicer` module used by `brainsfit_batch_registration.py` is supplied by 3D Slicer; installing the unrelated PyPI package named `slicer` is not a substitute.

### R

The R scripts require:

```r
install.packages(c("tidyverse", "irr"))
```

## External projects

### RESPECT Processing Module

The quantitative-processing workflow was adapted from the RESPECT Processing Module.

- [RESPECT Processing Module](https://github.com/Computer-Assisted-Clinical-Medicine/RESPECT_Processing_Module)

### RESPECT Co-Registration Module

The initial ITK/Elastix workflow uses helper functionality and conventions associated with the RESPECT Co-Registration Module.

- [RESPECT Co-Registration Module](https://github.com/Computer-Assisted-Clinical-Medicine/RESPECT_Co-Registration_Module)

### UKAT

Quantitative mapping uses UKAT functionality, including T1 mapping, T2 StimFit, B0 mapping, and supporting utilities.

- [UKAT: UKRIN-MAPS](https://github.com/UKRIN-MAPS/ukat)

### 3D Slicer and BRAINSFit

The final quantitative-map registration was performed through 3D Slicer using BRAINSFit.

- [3D Slicer](https://www.slicer.org/)
- [BRAINSFit documentation](https://www.slicer.org/wiki/Documentation/Nightly/Modules/BRAINSFit)

## Data availability and privacy

Raw DICOM data, kidney masks, participant identifiers, subject-level maps, registered images, ROI maps, input CSV files, descriptive tables, and generated figures are not included because they are subject to study-data access requirements and participant privacy restrictions.

Reproducing the complete workflow requires authorised access to the study data and the relevant upstream modules.

## Citation and acknowledgement

When reusing or adapting this workflow, cite the relevant RESPECT modules, UKAT, 3D Slicer/BRAINSFit, ITK/Elastix, the statistical methods, and the R packages used. Consult the thesis reference list and upstream projects for the appropriate citations and licence terms.
