# Renal MRI Reliability: Bachelor Thesis Code

This repository contains thesis-specific Python and R scripts used to
evaluate quantitative renal MRI biomarkers in healthy volunteers and
patients with chronic kidney disease (CKD). It is a study-specific
research workflow, not a standalone clinical application.

The analysis covers respiratory-triggered T1 mapping, MOLLI T1 mapping,
and respiratory-triggered T2 mapping. Repeatability and reproducibility
are assessed using the coefficient of variation (CoV), intraclass
correlation coefficient (ICC), and Bland--Altman analysis. Healthy
volunteers and CKD patients are analysed separately.

## Repository contents

``` text
.
├── README.md
├── requirements.txt
├── functions2.py
├── Rigid_registration.py
├── brainsfit_batch_registration.py
├── Analysis.py
├── Descriptive_roi_values.py
├── stat_T1T2_by_group.R
└── distribution_plots_separate.R
```

  -----------------------------------------------------------------------
  Script                              Role
  ----------------------------------- -----------------------------------
  `functions2.py`                     Quantitative-map processing for
                                      respiratory-triggered T1, MOLLI T1,
                                      T2 StimFit, B0/B1 maps, B1
                                      correction, and QC outputs.

  `Rigid_registration.py`             Initial ITK/Elastix registration
                                      workflow retained from thesis
                                      development.

  `brainsfit_batch_registration.py`   Final batch registration of
                                      quantitative maps to prepared
                                      anatomical T1-weighted references
                                      using 3D Slicer BRAINSFit.

  `Analysis.py`                       Applies whole-kidney masks to
                                      registered maps and creates four
                                      paired reliability-analysis input
                                      tables.

  `Descriptive_roi_values.py`         Generates examination-level and
                                      centre-summary descriptive
                                      whole-kidney tables from all
                                      available processed measurements.

  `stat_T1T2_by_group.R`              Calculates CoV, ICC, and
                                      Bland--Altman statistics separately
                                      for healthy volunteers and CKD
                                      patients.

  `distribution_plots_separate.R`     Creates centre-wise parameter
                                      distribution plots for the two
                                      populations.
  -----------------------------------------------------------------------

Only thesis-relevant scripts are included. Study data, generated
results, temporary development files, and upstream helper modules are
excluded.

## Execution guide

The complete workflow is:

``` text
Authorised DICOM data + upstream helper modules
                       |
                       v
                 functions2.py
          Map generation, B1 correction, QC
                       |
                       v
             Registration preparation
       (anatomical T1-weighted references)
                       |
            +----------+-----------+
            |                      |
            v                      v
   Rigid_registration.py   brainsfit_batch_registration.py
   Initial ITK/Elastix     Final 3D Slicer BRAINSFit
            |                      |
            +----------+-----------+
                       |
                       v
                  Analysis.py
       Whole-kidney ROI extraction and paired CSVs
                       |
            +----------+-----------+
            |                      |
            v                      v
   stat_T1T2_by_group.R   Descriptive_roi_values.py
    Reliability results    All-available measurements
                                    |
                                    v
                       distribution_plots_separate.R
```

The data directories are study-specific and are not included. Adapt
local paths to the authorised dataset and the conventions expected by
each script.

### 1. Prepare the environment and inputs

-   Install the Python packages listed in `requirements.txt`.
-   Install R and the packages listed under [Software
    requirements](#software-requirements).
-   Install 3D Slicer for final BRAINSFit registration.
-   Ensure authorised access to RESPECT data and the required upstream
    helper modules. These modules are not redistributed.
-   Prepare the local dataset using the path and filename conventions
    expected by the scripts.

### 2. Generate maps and QC

Run `functions2.py` in the configured processing environment. It
identifies available sequences using study-specific series names and
DICOM metadata. Missing or unprocessable sequences are skipped and
recorded in the processing log.

Final processing: - **Respiratory-triggered T1:** three-parameter
inversion-recovery fitting followed by B1 correction. - **MOLLI T1:**
three-parameter fitting followed by B1 correction. -
**Respiratory-triggered T2:** UKAT StimFit with vendor-specific models.
StimFit accounts for stimulated-echo and B1-related refocusing effects;
no separate post-fit B1 correction is applied. - **B0:** maps are
generated and stored but are not used in the final reliability
analysis. - **QC:** parameter maps, fit-quality outputs, numerical
summaries, and processing logs are generated where applicable.

### 3. Register quantitative maps

`Rigid_registration.py` contains the initial ITK/Elastix workflow and
relies on helper functionality associated with the RESPECT
Co-Registration Module. It is retained for provenance; final reported
ROI measurements use BRAINSFit-registered maps.

Run `brainsfit_batch_registration.py` inside 3D Slicer, not standard
Python. It registers each unmasked quantitative map to the prepared
T1-weighted fixed image for the same examination. Settings:
geometry-based initialisation, sampling percentage `0.002`, rigid and
global-scale stages enabled, affine and B-spline stages disabled, linear
interpolation, and background fill value `0`. This is rigid-plus-scale
linear registration, not deformable registration.

Expected inputs within each processed examination:

``` text
Registration/T1/fixed.nii.gz
Registration/T2/fixed.nii.gz
Registration/MOLLI/fixed.nii.gz

T1_RespTrig/nifti/t1_map_b1corr.nii.gz
T2_RespTrig/nifti/stimfit_t2_map.nii.gz
MOLLI/nifti/molli_t1_map_b1corr.nii.gz
```

Final outputs:

``` text
Registration_BRAINS/T1/t1_map_registered_brains.nii.gz
Registration_BRAINS/T2/stimfit_t2_map_registered_brains.nii.gz
Registration_BRAINS/MOLLI/molli_t1_map_registered_brains.nii.gz
```

Existing registered outputs are protected by default, missing inputs are
skipped, and each parameter/examination attempt is recorded in
`BRAINSFit_batch_registration_log.csv`.

**QC checkpoint:** visually inspect registered maps and transferred
masks across all retained slices, especially outer slices. Successful
execution alone does not establish anatomical validity.

### 4. Extract whole-kidney ROIs

Run `Analysis.py` after registration and mask preparation. It resamples
each whole-kidney mask to map geometry using nearest-neighbour
interpolation and saves:

``` text
t1_roi_brainsfit.nii.gz
molli_roi_brainsfit.nii.gz
t2_roi_brainsfit.nii.gz
```

Whole-kidney means use finite voxels within these validity ranges:

``` text
T1:        500–2500 ms
MOLLI T1:  500–2500 ms
T2:         15–150 ms
```

The script generates:

``` text
Repeatability_HealthyVolunteers_input.csv
Reproducibility_HealthyVolunteers_input.csv
Repeatability_Patients_input.csv
Reproducibility_Patients_input.csv
```

Centre-specific visit labels are mapped to the common visit columns
expected by R. A table may retain a single available visit, but the
reliability script forms complete pairs and excludes unmatched visits
from reliability calculations.

### 5. Run reliability statistics

Run `stat_T1T2_by_group.R` using the four input tables. It calculates
overall, centre-specific, and vendor-specific CoV and ICC with
confidence intervals, plus Bland--Altman mean relative differences and
95% limits of agreement. Tables and figures are generated separately for
healthy volunteers and CKD patients.

The CoV implementation follows the repeated-measurement formulation
described by de Boer et al., adapted to this thesis dataset.

### 6. Generate descriptive tables and distributions

Run `Descriptive_roi_values.py` independently of paired reliability
analysis. It reads all available BRAINSFit ROI maps and intentionally
retains repeated examinations. Outputs:

``` text
Descriptive_HealthyVolunteers_master.csv
Descriptive_CKDPatients_master.csv
Descriptive_HealthyVolunteers_centre_summary.csv
Descriptive_CKDPatients_centre_summary.csv
```

Master tables contain participant, population, centre, vendor, raw
visit, parameter, examination-level mean, valid-voxel count, and
inclusion status. Summary tables report measurement and
unique-participant counts, mean, standard deviation, median, minimum,
and maximum by centre and parameter, including an `All included centres`
row.

These summarise all available processed measurements, not single-visit
population reference values. Run `distribution_plots_separate.R` using
the two descriptive master tables; it does not use the paired
reliability inputs.

## Software requirements

### Python

``` bash
pip install -r requirements.txt
```

The workflow additionally requires authorised RESPECT processing and
co-registration helper functionality, UKAT, and ITK with Elastix support
for the initial registration. Final registration requires 3D Slicer with
BRAINSFit.

The `slicer` module used by `brainsfit_batch_registration.py` is
supplied by 3D Slicer; the unrelated PyPI package named `slicer` is not
a substitute.

### R

``` r
install.packages(c("tidyverse", "irr"))
```

## External projects

-   [RESPECT Processing
    Module](https://github.com/Computer-Assisted-Clinical-Medicine/RESPECT_Processing_Module)
-   [RESPECT Co-Registration
    Module](https://github.com/Computer-Assisted-Clinical-Medicine/RESPECT_Co-Registration_Module)
-   [UKAT: UKRIN-MAPS](https://github.com/UKRIN-MAPS/ukat)
-   [3D Slicer](https://www.slicer.org/)
-   [BRAINSFit
    documentation](https://www.slicer.org/wiki/Documentation/Nightly/Modules/BRAINSFit)

The workflow adapts and connects functionality from these projects and
established statistical methods. Restricted upstream code and study data
are not redistributed.

## Data availability and privacy

Raw DICOM data, kidney masks, participant identifiers, subject-level
maps, registered images, ROI maps, input CSV files, descriptive tables,
and generated figures are excluded because they are subject to
study-data access requirements and participant privacy restrictions.
Reproducing the complete workflow requires authorised study-data access
and relevant upstream modules.

## Citation and acknowledgement

When reusing or adapting this workflow, cite the relevant RESPECT
modules, UKAT, 3D Slicer/BRAINSFit, ITK/Elastix, statistical methods,
and R packages. Consult the thesis reference list and upstream projects
for appropriate citations and licence terms.
