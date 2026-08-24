"""Apply kidney masks to BRAINSFit-registered maps and build input tables.

This is an independent analysis branch for the final 3D Slicer BRAINSFit
registration. It never overwrites the original Elastix ROI maps or the CSV
tables stored in ``Analysis/Separate``.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk


# ============================================================
# COMMON SETTINGS
# ============================================================

def find_project_root():
    """Locate the data-analysis root in script and notebook execution modes."""
    if "__file__" in globals():
        start_folder = Path(__file__).resolve().parent
    else:
        start_folder = Path.cwd().resolve()

    for candidate in (start_folder, *start_folder.parents):
        if (
            (candidate / "Processed - Results").is_dir()
            and (candidate / "Masks").is_dir()
        ):
            return candidate

    raise FileNotFoundError(
        "Could not locate the data-analysis root containing both "
        "'Processed - Results' and 'Masks'."
    )


ROOT = find_project_root()

PROCESSED_RESULTS = ROOT / "Processed - Results"
MASKS_FOLDER = ROOT / "Masks"

OUTPUT_FOLDER = (
    ROOT / "code or notebooks" / "Analysis" / "BRAINSFit_Registration"
)

OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)


# ============================================================
# PART 1
# APPLY KIDNEY MASKS TO REGISTERED PARAMETER MAPS
# ============================================================

REGISTERED_MAPS = {
    "T1": (
        "Registration_BRAINS/T1/"
        "t1_map_registered_brains.nii.gz"
    ),
    "T2": (
        "Registration_BRAINS/T2/"
        "stimfit_t2_map_registered_brains.nii.gz"
    ),
    "MOLLI": (
        "Registration_BRAINS/MOLLI/"
        "molli_t1_map_registered_brains.nii.gz"
    ),
}


ROI_OUTPUT_NAMES = {
    "T1": "t1_roi_brainsfit.nii.gz",
    "T2": "t2_roi_brainsfit.nii.gz",
    "MOLLI": "molli_roi_brainsfit.nii.gz",
}


def apply_masks():

    print("-" * 70)
    print("PART 1: APPLYING KIDNEY MASKS")
    print("-" * 70)

    masked_maps = 0
    missing_masks = 0
    missing_maps = 0

    for patient_folder in sorted(
        PROCESSED_RESULTS.iterdir()
    ):

        if not patient_folder.is_dir():
            continue

        patient_id = patient_folder.name

        # ----------------------------------------------------
        # Find kidney mask
        # ----------------------------------------------------

        mask_path = (
            MASKS_FOLDER /
            f"{patient_id}.nii"
        )

        if not mask_path.exists():

            # Optional support for .nii.gz masks
            alternative_mask = (
                MASKS_FOLDER /
                f"{patient_id}.nii.gz"
            )

            if alternative_mask.exists():
                mask_path = alternative_mask

            else:
                print(
                    f"{patient_id}: mask not found."
                )

                missing_masks += 1
                continue

        # ----------------------------------------------------
        # Analysis output folder
        # ----------------------------------------------------

        analysis_folder = (
            patient_folder /
            "Analysis"
        )

        analysis_folder.mkdir(
            parents=True,
            exist_ok=True
        )

        print(f"\nMasking {patient_id}")

        # Read original kidney mask once
        mask = sitk.ReadImage(
            str(mask_path)
        )

        # ----------------------------------------------------
        # Apply to all available registered parameter maps
        # ----------------------------------------------------

        for parameter, relative_path in (
            REGISTERED_MAPS.items()
        ):

            map_path = (
                patient_folder /
                relative_path
            )

            if not map_path.exists():

                print(
                    f"  {parameter}: "
                    "registered map not found."
                )

                missing_maps += 1
                continue

            # Read registered quantitative map
            image = sitk.ReadImage(
                str(map_path)
            )

            # ------------------------------------------------
            # Resample kidney mask into map geometry
            # ------------------------------------------------

            resampler = (
                sitk.ResampleImageFilter()
            )

            resampler.SetReferenceImage(
                image
            )

            resampler.SetInterpolator(
                sitk.sitkNearestNeighbor
            )

            resampler.SetDefaultPixelValue(
                0
            )

            mask_resampled = (
                resampler.Execute(mask)
            )

            # ------------------------------------------------
            # Apply binary whole-kidney mask
            # ------------------------------------------------

            image_array = (
                sitk.GetArrayFromImage(
                    image
                )
            )

            mask_array = (
                sitk.GetArrayFromImage(
                    mask_resampled
                )
            )

            roi_array = np.where(
                mask_array > 0,
                image_array,
                0
            )

            roi_image = (
                sitk.GetImageFromArray(
                    roi_array
                )
            )

            roi_image.CopyInformation(
                image
            )

            output_path = (
                analysis_folder /
                ROI_OUTPUT_NAMES[
                    parameter
                ]
            )

            sitk.WriteImage(
                roi_image,
                str(output_path)
            )

            masked_maps += 1

            print(
                f"  {parameter}: "
                f"saved {output_path.name}"
            )

    print("\n" + "=" * 70)
    print("MASK APPLICATION FINISHED")
    print("=" * 70)

    print(
        f"ROI maps created : "
        f"{masked_maps}"
    )

    print(
        f"Missing masks    : "
        f"{missing_masks}"
    )

    print(
        f"Missing maps     : "
        f"{missing_maps}"
    )


# ============================================================
# PART 2
# EXTRACT ROI MEANS AND BUILD THE FOUR ANALYSIS TABLES
# ============================================================


ROI_FILENAMES = {
    "T1": "t1_roi_brainsfit.nii.gz",
    "T2": "t2_roi_brainsfit.nii.gz",
    "MOLLI": "molli_roi_brainsfit.nii.gz",
}


VALUE_RANGES = {
    "T1": (500, 2500),
    "MOLLI": (500, 2500),
    "T2": (15, 150),
}


VISIT_MAPPINGS = {

    "01": {
        "repeatability": {
            "v2a": "v2a",
            "v2b": "v3",
        },
        "reproducibility": {
            "v2a": "v2a",
            "v3": "v3",
        },
    },

    "02": {
        "repeatability": {
            "01": "v2a",
            "02": "v3",
        },
        "reproducibility": {
            "01": "v2a",
            "03": "v3",
        },
    },

    "03": {
        "repeatability": {
            "v1a": "v2a",
            "v2a": "v3",
        },
        "reproducibility": {
            "v1a": "v2a",
            "v3a": "v3",
        },
    },

    "04": {
        "repeatability": {
            "v2a": "v2a",
            "v2b": "v3",
        },
        "reproducibility": {
            "v2a": "v2a",
            "v3": "v3",
        },
    },
}


TABLE_FILENAMES = {

    (
        "repeatability",
        "HealthyVolunteers",
    ):
        "Repeatability_HealthyVolunteers_input.csv",

    (
        "reproducibility",
        "HealthyVolunteers",
    ):
        "Reproducibility_HealthyVolunteers_input.csv",

    (
        "repeatability",
        "Patients",
    ):
        "Repeatability_Patients_input.csv",

    (
        "reproducibility",
        "Patients",
    ):
        "Reproducibility_Patients_input.csv",
}


TABLE_COLUMNS = [
    "ID",
    "site",
    "visit",
    "anatomy",
    "vendor",
    "parameter",
    "value",
]


def parse_patient_folder(
    folder_name
):

    if "-" in folder_name:

        subject, raw_visit = (
            folder_name.split(
                "_",
                1
            )
        )

        subject_number, site = (
            subject.split(
                "-",
                1
            )
        )

    else:

        subject_number, site, raw_visit = (
            folder_name.split(
                "_",
                2
            )
        )

        subject = (
            f"{subject_number}-{site}"
        )

    if site not in VISIT_MAPPINGS:

        raise ValueError(
            f"Unsupported centre: {site}"
        )

    group = (
        "HealthyVolunteers"
        if int(subject_number) >= 500
        else "Patients"
    )

    vendor = (
        "GE"
        if site in {"01", "03"}
        else "Siemens"
    )

    return (
        subject,
        site,
        raw_visit,
        group,
        vendor,
    )


def extract_mean_roi_value(
    roi_image,
    parameter
):

    image_array = (
        sitk.GetArrayFromImage(
            roi_image
        )
    )

    minimum, maximum = (
        VALUE_RANGES[
            parameter
        ]
    )

    values = image_array[
        np.isfinite(
            image_array
        )
        &
        (
            image_array
            >= minimum
        )
        &
        (
            image_array
            <= maximum
        )
    ]

    if values.size == 0:
        return None

    return float(
        np.mean(values)
    )


def add_measurement(
    tables,
    analysis,
    group,
    subject,
    site,
    raw_visit,
    vendor,
    parameter,
    value,
):

    standardised_visit = (
        VISIT_MAPPINGS[
            site
        ][
            analysis
        ].get(
            raw_visit
        )
    )

    # This visit does not belong to
    # this particular analysis.
    if standardised_visit is None:
        return

    # IMPORTANT:
    # No paired-visit check here.
    # A single available visit is retained.
    tables[
        (
            analysis,
            group,
        )
    ].append(
        {
            "ID": subject,
            "site": int(site),
            "visit": standardised_visit,
            "anatomy": "whole",
            "vendor": vendor,
            "parameter": parameter,
            "value": value,
        }
    )


def build_tables():

    print("=" * 70)
    print("PART 2: BUILDING ANALYSIS TABLES")
    print("=" * 70)

    tables = {
        table_key: []
        for table_key
        in TABLE_FILENAMES
    }

    for patient_folder in sorted(
        PROCESSED_RESULTS.iterdir()
    ):

        if not patient_folder.is_dir():
            continue

        try:

            (
                subject,
                site,
                raw_visit,
                group,
                vendor,
            ) = parse_patient_folder(
                patient_folder.name
            )

        except (
            ValueError,
            IndexError,
        ):

            print(
                f"{patient_folder.name}: "
                "unrecognised folder name, "
                "skipping."
            )

            continue

        analysis_folder = (
            patient_folder /
            "Analysis"
        )

        if not analysis_folder.exists():
            continue

        print(
            f"Processing "
            f"{patient_folder.name}"
        )

        for (
            parameter,
            roi_filename,
        ) in ROI_FILENAMES.items():

            roi_path = (
                analysis_folder /
                roi_filename
            )

            if not roi_path.exists():
                continue

            roi_image = (
                sitk.ReadImage(
                    str(roi_path)
                )
            )

            mean_value = (
                extract_mean_roi_value(
                    roi_image,
                    parameter,
                )
            )

            if mean_value is None:

                print(
                    f"  {parameter}: "
                    "no valid voxels."
                )

                continue

            for analysis in (
                "repeatability",
                "reproducibility",
            ):

                add_measurement(
                    tables=tables,
                    analysis=analysis,
                    group=group,
                    subject=subject,
                    site=site,
                    raw_visit=raw_visit,
                    vendor=vendor,
                    parameter=parameter,
                    value=mean_value,
                )

            print(
                f"  {parameter}: "
                f"{mean_value:.2f} ms"
            )

    # --------------------------------------------------------
    # Save four tables
    # --------------------------------------------------------

    for (
        table_key,
        filename,
    ) in TABLE_FILENAMES.items():

        table = pd.DataFrame(
            tables[
                table_key
            ],
            columns=TABLE_COLUMNS,
        )

        table = (
            table
            .sort_values(
                [
                    "site",
                    "ID",
                    "parameter",
                    "visit",
                ],
                kind="stable",
            )
            .reset_index(
                drop=True
            )
        )

        # ----------------------------------------------------
        # Duplicate protection
        # ----------------------------------------------------

        duplicate_mask = (
            table.duplicated(
                subset=[
                    "ID",
                    "site",
                    "visit",
                    "parameter",
                ],
                keep=False,
            )
        )

        if duplicate_mask.any():

            duplicates = (
                table.loc[
                    duplicate_mask,
                    [
                        "ID",
                        "site",
                        "visit",
                        "parameter",
                    ],
                ]
            )

            raise RuntimeError(
                f"Duplicate measurements "
                f"generated for "
                f"{filename}:\n"
                f"{duplicates}"
            )

        output_path = (
            OUTPUT_FOLDER /
            filename
        )

        table.to_csv(
            output_path,
            index=False,
        )

        print(
            f"Saved {filename}: "
            f"{len(table)} rows"
        )

    print("\nFour tables created.")


def main():
    apply_masks()
    build_tables()


if __name__ == "__main__":
    main()
