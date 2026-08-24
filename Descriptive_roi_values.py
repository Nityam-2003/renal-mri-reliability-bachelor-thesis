"""
Create descriptive whole-kidney biomarker-value tables.

This thesis-specific script is independent of the repeatability and
reproducibility analyses. It scans every available masked ROI map in
``Processed - Results`` and therefore does not use visit-pair definitions or
the four reliability-analysis input tables.

For each processed examination and available parameter, the script calculates
the mean of all valid voxels in the already masked whole-kidney ROI map. It
writes the following private data files to ``Analysis/BRAINSFit_Registration``:

* Descriptive_HealthyVolunteers_master.csv
* Descriptive_CKDPatients_master.csv
* Descriptive_HealthyVolunteers_centre_summary.csv
* Descriptive_CKDPatients_centre_summary.csv

The master files retain every available examination and include the raw visit
label and valid-voxel count. The summary files provide one row per imaging
centre and parameter, followed by an ``All included centres`` row for each
parameter. Repeated visits are intentionally retained: these tables describe
all available processed measurements and are not population reference values
based on one examination per participant.

Run this script from ``Analysis/BRAINSFit_Registration`` after ROI map
generation. It reads
restricted thesis data and must not be committed to a public repository.

Requirements: Python, NumPy, pandas, and SimpleITK.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk


def find_project_root():
    """Find the data-analysis root in script and notebook execution modes."""
    if "__file__" in globals():
        start_folder = Path(__file__).resolve().parent
    else:
        start_folder = Path.cwd().resolve()

    for candidate in (start_folder, *start_folder.parents):
        if (candidate / "Processed - Results").is_dir():
            return candidate

    raise FileNotFoundError(
        "Could not locate the project root. Start the notebook from within "
        "the 'Everything related to data and analysis' folder or one of its "
        "subfolders."
    )


PROJECT_ROOT = find_project_root()
PROCESSED_RESULTS = PROJECT_ROOT / "Processed - Results"
OUTPUT_FOLDER = (
    PROJECT_ROOT
    / "code or notebooks"
    / "Analysis"
    / "BRAINSFit_Registration"
)

ROI_FILENAMES = {
    "T1": "t1_roi_brainsfit.nii.gz",
    "MOLLI": "molli_roi_brainsfit.nii.gz",
    "T2": "t2_roi_brainsfit.nii.gz",
}

# These are the same validity ranges used to rebuild the reliability input
# tables. Zero-valued background voxels in masked ROI maps are excluded.
VALUE_RANGES = {
    "T1": (500, 2500),
    "MOLLI": (500, 2500),
    "T2": (15, 150),
}

MASTER_COLUMNS = [
    "ID",
    "population",
    "site",
    "vendor",
    "raw_visit",
    "parameter",
    "mean_roi_value_ms",
    "valid_voxel_count",
    "status",
]


def parse_patient_folder(folder_name):
    """Parse names such as ``001-01_v2a`` and ``001_02_01``.

    IDs beginning with 5 are healthy volunteers. All other IDs are CKD
    patients. Centres 01 and 03 used GE scanners; Centres 02 and 04 used
    Siemens scanners.
    """
    if "-" in folder_name:
        subject, raw_visit = folder_name.split("_", 1)
        subject_number, site = subject.split("-", 1)
    else:
        subject_number, site, raw_visit = folder_name.split("_", 2)
        subject = f"{subject_number}-{site}"

    if site not in {"01", "02", "03", "04"}:
        raise ValueError(f"Unsupported centre: {site}")

    population = "HealthyVolunteers" if int(subject_number) >= 500 else "CKDPatients"
    vendor = "GE" if site in {"01", "03"} else "Siemens"

    return subject, int(site), raw_visit, population, vendor


def extract_mean_roi_value(roi_path, parameter):
    """Calculate a mean from valid voxels in an existing masked ROI map."""
    roi_image = sitk.ReadImage(str(roi_path))
    image_array = sitk.GetArrayFromImage(roi_image)
    minimum, maximum = VALUE_RANGES[parameter]
    values = image_array[
        np.isfinite(image_array)
        & (image_array >= minimum)
        & (image_array <= maximum)
    ]

    if values.size == 0:
        return np.nan, 0, "no_valid_voxels"

    return float(np.mean(values)), int(values.size), "included"


def build_master_table():
    """Read every available ROI map, without pairing or visit selection."""
    rows = []

    for examination_folder in sorted(PROCESSED_RESULTS.iterdir()):
        if not examination_folder.is_dir():
            continue

        try:
            subject, site, raw_visit, population, vendor = parse_patient_folder(
                examination_folder.name
            )
        except (ValueError, IndexError):
            print(f"{examination_folder.name}: unrecognised folder name, skipping.")
            continue

        analysis_folder = examination_folder / "Analysis"
        if not analysis_folder.is_dir():
            continue

        for parameter, roi_filename in ROI_FILENAMES.items():
            roi_path = analysis_folder / roi_filename
            if not roi_path.exists():
                continue

            mean_value, valid_voxel_count, status = extract_mean_roi_value(
                roi_path, parameter
            )
            rows.append(
                {
                    "ID": subject,
                    "population": population,
                    "site": site,
                    "vendor": vendor,
                    "raw_visit": raw_visit,
                    "parameter": parameter,
                    "mean_roi_value_ms": mean_value,
                    "valid_voxel_count": valid_voxel_count,
                    "status": status,
                }
            )

    master_table = pd.DataFrame(rows, columns=MASTER_COLUMNS)
    if master_table.empty:
        raise RuntimeError("No masked ROI maps were found in Processed - Results.")

    return master_table.sort_values(
        ["population", "site", "ID", "raw_visit", "parameter"], kind="stable"
    ).reset_index(drop=True)


def summarise_population(master_table, population):
    """Summarise all valid examination-level values by centre and parameter."""
    valid_measurements = master_table.loc[
        (master_table["population"] == population)
        & (master_table["status"] == "included")
    ].copy()

    summary_rows = []
    for (site, vendor, parameter), group in valid_measurements.groupby(
        ["site", "vendor", "parameter"], sort=True
    ):
        summary_rows.append(
            {
                "centre": f"{site:02d}",
                "vendor": vendor,
                "parameter": parameter,
                "measurements_n": len(group),
                "participants_n": group["ID"].nunique(),
                "mean_ms": group["mean_roi_value_ms"].mean(),
                "sd_ms": group["mean_roi_value_ms"].std(ddof=1),
                "median_ms": group["mean_roi_value_ms"].median(),
                "minimum_ms": group["mean_roi_value_ms"].min(),
                "maximum_ms": group["mean_roi_value_ms"].max(),
            }
        )

    for parameter, group in valid_measurements.groupby("parameter", sort=False):
        summary_rows.append(
            {
                "centre": "All included centres",
                "vendor": "Mixed",
                "parameter": parameter,
                "measurements_n": len(group),
                "participants_n": group["ID"].nunique(),
                "mean_ms": group["mean_roi_value_ms"].mean(),
                "sd_ms": group["mean_roi_value_ms"].std(ddof=1),
                "median_ms": group["mean_roi_value_ms"].median(),
                "minimum_ms": group["mean_roi_value_ms"].min(),
                "maximum_ms": group["mean_roi_value_ms"].max(),
            }
        )

    summary = pd.DataFrame(summary_rows)
    parameter_order = pd.CategoricalDtype(["T1", "MOLLI", "T2"], ordered=True)
    summary["parameter"] = summary["parameter"].astype(parameter_order)
    summary["_all_centres"] = summary["centre"].eq("All included centres")
    summary = summary.sort_values(
        ["_all_centres", "centre", "parameter"], kind="stable"
    ).drop(columns="_all_centres")

    return summary.reset_index(drop=True)


def main():
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
    master_table = build_master_table()

    for population in ("HealthyVolunteers", "CKDPatients"):
        population_master = master_table.loc[
            master_table["population"] == population
        ].copy()
        summary_table = summarise_population(master_table, population)

        master_path = OUTPUT_FOLDER / f"Descriptive_{population}_master.csv"
        summary_path = OUTPUT_FOLDER / f"Descriptive_{population}_centre_summary.csv"
        population_master.to_csv(master_path, index=False)
        summary_table.to_csv(summary_path, index=False, float_format="%.3f")

        included = (population_master["status"] == "included").sum()
        skipped = (population_master["status"] != "included").sum()
        print(
            f"{population}: saved {len(population_master)} available ROI maps "
            f"({included} included, {skipped} without valid voxels)."
        )
        print(f"  Master table: {master_path}")
        print(f"  Centre summary: {summary_path}")


if __name__ == "__main__":
    main()
