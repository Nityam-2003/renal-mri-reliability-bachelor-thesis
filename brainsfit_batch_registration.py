"""Batch-register renal quantitative maps with 3D Slicer BRAINSFit.

This script must be executed by 3D Slicer, not by a normal Python or Jupyter
installation. It registers each unmasked quantitative map directly to the
already prepared fixed T1-weighted image for the same dataset.

Existing registration results are never overwritten. New results are written
under each dataset's Registration_BRAINS directory.
"""

from __future__ import annotations

import csv
import traceback
from datetime import datetime
from pathlib import Path

import slicer


PROCESSED_RESULTS = Path(
    r"C:\Users\cr714\Downloads\Bachelorarbeit\Everything related to data and analysis\Processed - Results"
)

# These settings reproduce the visible General Registration (BRAINS) settings
# supplied by the user. The stages run in the listed order: rigid, then global
# scale. Change these only after validating a representative group of cases.
SAMPLING_PERCENTAGE = 0.002
USE_RIGID = True
USE_GLOBAL_SCALE = True
INITIALISATION_MODE = "useGeometryAlign"

# False protects an already completed BRAINSFit output if the script is rerun.
OVERWRITE_EXISTING = False


PARAMETERS = {
    "T1": {
        "fixed": Path("Registration/T1/fixed.nii.gz"),
        "map": Path("T1_RespTrig/nifti/t1_map_b1corr.nii.gz"),
        "output_dir": Path("Registration_BRAINS/T1"),
        "output_map": "t1_map_registered_brains.nii.gz",
    },
    "T2": {
        "fixed": Path("Registration/T2/fixed.nii.gz"),
        "map": Path("T2_RespTrig/nifti/stimfit_t2_map.nii.gz"),
        "output_dir": Path("Registration_BRAINS/T2"),
        "output_map": "stimfit_t2_map_registered_brains.nii.gz",
    },
    "MOLLI": {
        "fixed": Path("Registration/MOLLI/fixed.nii.gz"),
        "map": Path("MOLLI/nifti/molli_t1_map_b1corr.nii.gz"),
        "output_dir": Path("Registration_BRAINS/MOLLI"),
        "output_map": "molli_t1_map_registered_brains.nii.gz",
    },
}


LOG_FIELDS = [
    "timestamp",
    "dataset",
    "parameter",
    "status",
    "fixed_image",
    "moving_map",
    "registered_map",
    "transform",
    "brainsfit_status",
    "message",
]


def append_log(log_path: Path, record: dict[str, str]) -> None:
    """Append one result immediately so progress survives an interrupted run."""
    new_file = not log_path.exists()
    with log_path.open("a", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=LOG_FIELDS)
        if new_file:
            writer.writeheader()
        writer.writerow(record)


def remove_nodes(*nodes) -> None:
    """Remove temporary MRML nodes and keep Slicer memory usage bounded."""
    for node in nodes:
        if node is not None and node.GetScene() is not None:
            slicer.mrmlScene.RemoveNode(node)


def register_one(dataset_dir: Path, parameter: str, specification: dict) -> dict[str, str]:
    fixed_path = dataset_dir / specification["fixed"]
    moving_path = dataset_dir / specification["map"]
    output_dir = dataset_dir / specification["output_dir"]
    output_path = output_dir / specification["output_map"]
    transform_path = output_dir / "brainsfit_transform.h5"

    base_record = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "dataset": dataset_dir.name,
        "parameter": parameter,
        "fixed_image": str(fixed_path),
        "moving_map": str(moving_path),
        "registered_map": str(output_path),
        "transform": str(transform_path),
    }

    if not fixed_path.exists() or not moving_path.exists():
        missing = []
        if not fixed_path.exists():
            missing.append("fixed image")
        if not moving_path.exists():
            missing.append("quantitative map")
        return {
            **base_record,
            "status": "SKIPPED_MISSING_INPUT",
            "brainsfit_status": "",
            "message": "Missing " + " and ".join(missing),
        }

    if output_path.exists() and not OVERWRITE_EXISTING:
        return {
            **base_record,
            "status": "SKIPPED_EXISTING",
            "brainsfit_status": "",
            "message": "Existing BRAINSFit output protected",
        }

    output_dir.mkdir(parents=True, exist_ok=True)

    fixed_node = None
    moving_node = None
    output_node = None
    transform_node = None
    cli_node = None

    try:
        fixed_node = slicer.util.loadVolume(str(fixed_path))
        moving_node = slicer.util.loadVolume(str(moving_path))
        if fixed_node is None or moving_node is None:
            raise RuntimeError("Slicer could not load one or both input volumes")

        output_node = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLScalarVolumeNode",
            f"BRAINS_{dataset_dir.name}_{parameter}_registered",
        )
        transform_node = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLLinearTransformNode",
            f"BRAINS_{dataset_dir.name}_{parameter}_transform",
        )

        cli_parameters = {
            "fixedVolume": fixed_node.GetID(),
            "movingVolume": moving_node.GetID(),
            "outputVolume": output_node.GetID(),
            "linearTransform": transform_node.GetID(),
            "initializeTransformMode": INITIALISATION_MODE,
            "samplingPercentage": SAMPLING_PERCENTAGE,
            "useRigid": USE_RIGID,
            "useScaleVersor3D": USE_GLOBAL_SCALE,
            "useScaleSkewVersor3D": False,
            "useAffine": False,
            "useBSpline": False,
            "interpolationMode": "Linear",
            "backgroundFillValue": 0,
        }

        cli_node = slicer.cli.runSync(
            slicer.modules.brainsfit,
            None,
            cli_parameters,
        )
        cli_status = cli_node.GetStatusString()

        if cli_status != "Completed":
            error_text = cli_node.GetErrorText() if hasattr(cli_node, "GetErrorText") else ""
            raise RuntimeError(f"BRAINSFit status: {cli_status}. {error_text}".strip())

        if output_node.GetImageData() is None:
            raise RuntimeError("BRAINSFit completed without producing an output volume")

        if not slicer.util.saveNode(output_node, str(output_path)):
            raise RuntimeError("Slicer could not save the registered quantitative map")
        if not slicer.util.saveNode(transform_node, str(transform_path)):
            raise RuntimeError("Slicer could not save the BRAINSFit transformation")

        return {
            **base_record,
            "status": "COMPLETED",
            "brainsfit_status": cli_status,
            "message": "",
        }

    except Exception as error:
        return {
            **base_record,
            "status": "FAILED",
            "brainsfit_status": cli_node.GetStatusString() if cli_node is not None else "",
            "message": f"{type(error).__name__}: {error}",
        }
    finally:
        remove_nodes(cli_node, transform_node, output_node, moving_node, fixed_node)
        slicer.app.processEvents()


def main() -> int:
    if not PROCESSED_RESULTS.is_dir():
        print(f"ERROR: Processed-results directory does not exist: {PROCESSED_RESULTS}")
        return 1

    log_path = PROCESSED_RESULTS / "BRAINSFit_batch_registration_log.csv"
    dataset_dirs = sorted(path for path in PROCESSED_RESULTS.iterdir() if path.is_dir())

    attempted = 0
    completed = 0
    failed = 0
    skipped_missing = 0
    skipped_existing = 0

    print("\nBRAINSFit renal-map batch registration")
    print(f"Processed-results folder: {PROCESSED_RESULTS}")
    print(f"Datasets found: {len(dataset_dirs)}")
    print(f"Log: {log_path}")
    print("Existing outputs will not be overwritten.\n")

    for dataset_index, dataset_dir in enumerate(dataset_dirs, start=1):
        for parameter, specification in PARAMETERS.items():
            attempted += 1
            print(
                f"[{dataset_index}/{len(dataset_dirs)}] "
                f"{dataset_dir.name} | {parameter}"
            )

            record = register_one(dataset_dir, parameter, specification)
            append_log(log_path, record)
            print(f"  {record['status']} {record['message']}")

            if record["status"] == "COMPLETED":
                completed += 1
            elif record["status"] == "FAILED":
                failed += 1
            elif record["status"] == "SKIPPED_MISSING_INPUT":
                skipped_missing += 1
            elif record["status"] == "SKIPPED_EXISTING":
                skipped_existing += 1

    print("\nBatch registration finished")
    print(f"Parameter/dataset combinations checked: {attempted}")
    print(f"Completed: {completed}")
    print(f"Failed: {failed}")
    print(f"Skipped because inputs were missing: {skipped_missing}")
    print(f"Skipped because outputs already existed: {skipped_existing}")
    print(f"Full log: {log_path}")

    return 0 if failed == 0 else 2


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception:
        traceback.print_exc()
        exit_code = 3

    slicer.app.processEvents()
    slicer.app.exit(exit_code)

