"""
IMPORTANT: Thesis-specific preparation of final analysis input tables.

Run this script from ``Analysis/Separate`` after processing and registration
have been completed and the matching private kidney masks are available. It
uses the folder conventions described in the repository README. The script
reads restricted maps and masks, and must not be run against or committed with
public data.

Create masked kidney ROI maps and final analysis input tables.

The script processes all available registered quantitative maps, applies the
matching kidney mask, and writes four separate long-format input tables:

* Repeatability_HealthyVolunteers_input.csv
* Reproducibility_HealthyVolunteers_input.csv
* Repeatability_Patients_input.csv
* Reproducibility_Patients_input.csv

The tables are regenerated from the available maps each time this script runs.
It requires Python with NumPy, pandas, and SimpleITK.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk


# This script is stored in: code or notebooks/Analysis/Separate/
PROJECT_ROOT = Path(__file__).resolve().parents[3]
PROCESSED_RESULTS = PROJECT_ROOT / "Processed - Results"
MASKS_FOLDER = PROJECT_ROOT / "Masks"
OUTPUT_FOLDER = Path(__file__).resolve().parent


MAPS = {
    "T1": "Registration/T1/t1_map_registered.nii.gz",
    "T2": "Registration/T2/stimfit_t2_map_registered.nii.gz",
    "MOLLI": "Registration/MOLLI/molli_map_registered.nii.gz",
}

ROI_FILENAMES = {
    "T1": "t1_roi.nii.gz",
    "T2": "t2_roi.nii.gz",
    "MOLLI": "molli_roi.nii.gz",
}

VALUE_RANGES = {
    "T1": (500, 3000),
    "MOLLI": (500, 3000),
    "T2": (15, 150),
}

# Raw folder visits are converted to the common visit labels used by the
# repeatability and reproducibility analyses.
VISIT_MAPPINGS = {
    "01": {
        "repeatability": {"v2a": "v2a", "v2b": "v3"},
        "reproducibility": {"v2a": "v2a", "v3": "v3"},
    },
    "02": {
        "repeatability": {"01": "v2a", "02": "v3"},
        "reproducibility": {"01": "v2a", "03": "v3"},
    },
    "03": {
        "repeatability": {"v1a": "v2a", "v2a": "v3"},
        "reproducibility": {"v1a": "v2a", "v3a": "v3"},
    },
    "04": {
        "repeatability": {"v2a": "v2a", "v2b": "v3"},
        "reproducibility": {"v2a": "v2a", "v3": "v3"},
    },
}

TABLE_FILENAMES = {
    ("repeatability", "HealthyVolunteers"):
        "Repeatability_HealthyVolunteers_input.csv",
    ("reproducibility", "HealthyVolunteers"):
        "Reproducibility_HealthyVolunteers_input.csv",
    ("repeatability", "Patients"):
        "Repeatability_Patients_input.csv",
    ("reproducibility", "Patients"):
        "Reproducibility_Patients_input.csv",
}

TABLE_COLUMNS = ["ID", "site", "visit", "anatomy", "vendor", "parameter", "value"]


def parse_patient_folder(folder_name):
    """Parse centre-specific processed-result folder names.

    Centres 01, 03, and 04 use names such as ``001-01_v2a``. Centre 02
    uses names such as ``001_02_01``. IDs starting with 5 denote healthy
    volunteers; the remaining IDs denote CKD patients.
    """
    if "-" in folder_name:
        subject, raw_visit = folder_name.split("_", 1)
        subject_number, site = subject.split("-", 1)
    else:
        subject_number, site, raw_visit = folder_name.split("_", 2)
        subject = f"{subject_number}-{site}"

    if site not in VISIT_MAPPINGS:
        raise ValueError(f"Unsupported centre: {site}")

    group = "HealthyVolunteers" if int(subject_number) >= 500 else "Patients"
    vendor = "GE" if site in {"01", "03"} else "Siemens"

    return subject, site, raw_visit, group, vendor


def apply_mask(image, mask):
    """Resample a label mask to map geometry and apply it to the map."""
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(image)
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    resampler.SetDefaultPixelValue(0)

    resampled_mask = resampler.Execute(mask)
    image_array = sitk.GetArrayFromImage(image)
    mask_array = sitk.GetArrayFromImage(resampled_mask)

    roi_array = np.where(mask_array > 0, image_array, 0)
    roi_image = sitk.GetImageFromArray(roi_array)
    roi_image.CopyInformation(image)

    return roi_image, image_array, mask_array


def extract_mean_roi_value(image_array, mask_array, parameter):
    """Return the mean value within the labelled kidney ROI and valid range."""
    minimum, maximum = VALUE_RANGES[parameter]
    values = image_array[
        (mask_array > 0)
        & np.isfinite(image_array)
        & (image_array >= minimum)
        & (image_array <= maximum)
    ]

    if values.size == 0:
        return None

    return float(np.mean(values))


def add_measurement(tables, analysis, group, subject, site, raw_visit, vendor, parameter, value):
    """Add a measurement only if its raw visit belongs to an analysis pair."""
    standardised_visit = VISIT_MAPPINGS[site][analysis].get(raw_visit)
    if standardised_visit is None:
        return

    tables[(analysis, group)].append(
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


def main():
    tables = {table_key: [] for table_key in TABLE_FILENAMES}

    for patient_folder in sorted(PROCESSED_RESULTS.iterdir()):
        if not patient_folder.is_dir():
            continue

        try:
            subject, site, raw_visit, group, vendor = parse_patient_folder(patient_folder.name)
        except (ValueError, IndexError):
            print(f"{patient_folder.name}: unrecognised folder name, skipping.")
            continue

        mask_path = MASKS_FOLDER / f"{patient_folder.name}.nii"
        if not mask_path.exists():
            print(f"{patient_folder.name}: mask not found, skipping.")
            continue

        mask = sitk.ReadImage(str(mask_path))
        analysis_folder = patient_folder / "Analysis"
        analysis_folder.mkdir(exist_ok=True)

        print(f"Processing {patient_folder.name}")

        for parameter, relative_map_path in MAPS.items():
            map_path = patient_folder / relative_map_path
            if not map_path.exists():
                print(f"  {parameter}: registered map not found.")
                continue

            image = sitk.ReadImage(str(map_path))
            roi_image, image_array, mask_array = apply_mask(image, mask)

            roi_path = analysis_folder / ROI_FILENAMES[parameter]
            sitk.WriteImage(roi_image, str(roi_path))

            mean_value = extract_mean_roi_value(image_array, mask_array, parameter)
            if mean_value is None:
                print(f"  {parameter}: no valid voxels in the masked ROI.")
                continue

            for analysis in ("repeatability", "reproducibility"):
                add_measurement(
                    tables,
                    analysis,
                    group,
                    subject,
                    site,
                    raw_visit,
                    vendor,
                    parameter,
                    mean_value,
                )

            print(f"  {parameter}: ROI map saved; mean value = {mean_value:.2f}")

    for table_key, filename in TABLE_FILENAMES.items():
        table = pd.DataFrame(tables[table_key], columns=TABLE_COLUMNS)
        table = table.sort_values(["ID", "parameter", "visit"], kind="stable")
        output_path = OUTPUT_FOLDER / filename
        table.to_csv(output_path, index=False)
        print(f"Saved {filename}: {len(table)} rows")


if __name__ == "__main__":
    main()
