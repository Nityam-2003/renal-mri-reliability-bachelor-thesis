"""
IMPORTANT: Thesis-specific rigid registration of quantitative renal MRI maps.

This script registers available T1, T2, and MOLLI maps to the relevant
anatomical reference image for each processed dataset. It is an adaptation
around ITK/Elastix and helper functionality associated with the restricted
RESPECT Co-Registration Module; upstream helper code is not redistributed.

Run this file from the ``Registration`` folder using
``python Rigid_registration.py``. The relative paths beginning with ``../../``
assume the project layout described in the repository README. Existing
registered maps are skipped individually.

This script reads restricted data and creates outputs in ``Processed -
Results``. Do not commit raw data, masks, participant identifiers, or generated
outputs to a public repository.
"""

from pathlib import Path
import pydicom
import shutil
import itk
import util


T1_NAMES = ["Cor_T1map_RespTrig","T1mapping_RespTrig","T1MAPPING_RESPTRIG","T1MAPPING_FB"]
T2_NAMES = ["Cor_T2map_RespTrig","T2_mapping_RespTrig","T2_MAPPING_RESPTRIG"]
T1W_NAMES = ["Cor_T1weighted_BH", "T1W_BH", "T1W BH"]

def log_issue(patient_id, section, message):
    logfile = "../../Processed - Results/Pipeline_Log.txt"
    print(f"\n[{patient_id}] {section}: {message}\n\n")

    with open(logfile, "a") as f:
        f.write(
            f"{patient_id} | "
            f"{section} | "
            f"{message}\n"
            f"{'-'*125}\n")


def register_t1(patient_folder):
    print("\nT1 REGISTRATION")
    
    t1_folders = []
    for folder in patient_folder.iterdir():
        if folder.is_dir() and any(name in folder.name for name in T1_NAMES):
            t1_folders.append(folder)

    moving_folder = max(t1_folders, key=lambda folder: len(list(folder.iterdir())))
    print(f"Moving folder: {moving_folder}")
    
    registration_dir = processed_results / "Registration" / "T1"
    registration_dir.mkdir(parents=True, exist_ok=True)

    moving_5_files = registration_dir / "moving_5_files"
    if moving_5_files.exists():
        shutil.rmtree(moving_5_files)
    moving_5_files.mkdir()

    moving_files = sorted(f for f in moving_folder.iterdir() if f.is_file())

    for file in moving_files[:5]:
        shutil.copy2(file, moving_5_files / file.name)

    util.crop_image_to_match_fov(str(fixed_input), str(moving_5_files), str(registration_dir))

    fixed = itk.imread(str(registration_dir / "fixed.nii.gz"), itk.F)
    moving = itk.imread(str(registration_dir / "moving.nii.gz"), itk.F)

    parameter_object = itk.ParameterObject.New()
    parameter_object.AddParameterMap(parameter_object.GetDefaultParameterMap("rigid"))

    result_image, result_transform_parameters = itk.elastix_registration_method(fixed, moving, parameter_object=parameter_object, log_to_console=True)

    itk.imwrite(result_image, str(registration_dir / "t1_rigid.nii.gz"))

    t1_map = itk.imread(str(processed_results / "T1_RespTrig" / "nifti" / "t1_map_b1corr.nii.gz"), itk.F)

    registered_map = itk.transformix_filter(t1_map, result_transform_parameters)

    itk.imwrite(registered_map, str(registration_dir / "t1_map_registered.nii.gz"))

    if moving_5_files.exists():
        shutil.rmtree(moving_5_files)
        
    print("Registered T1 map created.")
    print("T1 registration completed.")

#--------------------------------------------------------------------------------------------------------------------------------------------------------

def register_t2(patient_folder):
    print("\nT2 REGISTRATION")
    
    patient_id = patient_folder.name

    t2_folders = []
    for folder in patient_folder.iterdir():
        if folder.is_dir() and any(name in folder.name for name in T2_NAMES):
            t2_folders.append(folder)

    moving_folder = max(t2_folders, key=lambda folder: len(list(folder.iterdir())))
    print(f"Moving folder: {moving_folder}")
    
    registration_dir = processed_results / "Registration" / "T2"
    registration_dir.mkdir(parents=True, exist_ok=True)

    moving_5_files = registration_dir / "moving_5_files"
    if moving_5_files.exists():
        shutil.rmtree(moving_5_files)
    moving_5_files.mkdir()

    moving_files = sorted(f for f in moving_folder.iterdir() if f.is_file())

    if "-01" in patient_id or "-03" in patient_id:
        for file in moving_files[::10]:
            shutil.copy2(file, moving_5_files / file.name)

    elif "-02" in patient_id or "_02" in patient_id or "-04" in patient_id:
        for file in moving_files[:5]:
            shutil.copy2(file, moving_5_files / file.name)

    util.crop_image_to_match_fov(str(fixed_input), str(moving_5_files), str(registration_dir))

    fixed = itk.imread(str(registration_dir / "fixed.nii.gz"), itk.F)
    moving = itk.imread(str(registration_dir / "moving.nii.gz"), itk.F)

    parameter_object = itk.ParameterObject.New()
    parameter_object.AddParameterMap(parameter_object.GetDefaultParameterMap("rigid"))

    result_image, result_transform_parameters = itk.elastix_registration_method(fixed, moving, parameter_object=parameter_object, log_to_console=True)

    itk.imwrite(result_image, str(registration_dir / "t2_rigid.nii.gz"))

    t2_map = itk.imread(str(processed_results / "T2_RespTrig" / "nifti" / "stimfit_t2_map.nii.gz"), itk.F)

    registered_map = itk.transformix_filter(t2_map, result_transform_parameters)

    itk.imwrite(registered_map, str(registration_dir / "stimfit_t2_map_registered.nii.gz"))

    if moving_5_files.exists():
        shutil.rmtree(moving_5_files)
        
    print("Registered T2 map created.")
    print("T2 registration completed.")

#--------------------------------------------------------------------------------------------------------------------------------------------------------  

def register_molli(patient_folder):

    print("\nMOLLI REGISTRATION")

    patient_id = patient_folder.name

    molli_folders = []

    for folder in patient_folder.iterdir():
        if folder.is_dir() and "MOLLI" in folder.name:
            molli_folders.append(folder)

    registration_dir = processed_results / "Registration" / "MOLLI"
    registration_dir.mkdir(parents=True, exist_ok=True)

    if "-01" in patient_id or "-03" in patient_id:
        molli_folders = sorted(molli_folders, key=lambda folder: pydicom.dcmread(str(sorted(folder.iterdir())[0])).SliceLocation)
    
        moving_5_files = registration_dir / "moving_5_files"
        if moving_5_files.exists():
            shutil.rmtree(moving_5_files)
        moving_5_files.mkdir()
    
        for i, folder in enumerate(molli_folders):
            first_file = sorted(f for f in folder.iterdir() if f.is_file())[0]
            shutil.copy2(first_file, moving_5_files / f"{i+1}.dcm")
        
        moving_input = [str(f) for f in sorted(moving_5_files.iterdir())]
            
        print(f"Moving folders:{[folder.name for folder in molli_folders]}")

    elif "-04" in patient_id:
        moving_folder = max(molli_folders, key=lambda folder: len(list(folder.iterdir())))

        moving_5_files = registration_dir / "moving_5_files"
        if moving_5_files.exists():
            shutil.rmtree(moving_5_files)
        moving_5_files.mkdir()

        moving_files = sorted(f for f in moving_folder.iterdir() if f.is_file())

        for file in moving_files[:5]:
            shutil.copy2(file, moving_5_files / file.name)

        moving_input = str(moving_5_files)

        print(f"Moving folder: {moving_folder.name}")

    util.crop_image_to_match_fov(str(fixed_input),moving_input,str(registration_dir))
    
    # Reuse the fixed reference generated during T2 registration.
    fixed = itk.imread(str(registration_dir/ "fixed.nii.gz"), itk.F)
    moving = itk.imread(str(registration_dir / "moving.nii.gz"), itk.F)

    parameter_object = itk.ParameterObject.New()
    parameter_object.AddParameterMap(parameter_object.GetDefaultParameterMap("rigid"))
    
    result_image, result_transform_parameters = itk.elastix_registration_method(
    fixed,
    moving,
    parameter_object=parameter_object,
    log_to_console=True)

    itk.imwrite(result_image, str(registration_dir / "molli_rigid.nii.gz"))

    molli_map = itk.imread(str(processed_results / "MOLLI" / "nifti" / "molli_t1_map_b1corr.nii.gz"),itk.F,)

    registered_map = itk.transformix_filter(molli_map, result_transform_parameters)

    itk.imwrite(registered_map,str(registration_dir / "molli_map_registered.nii.gz"),)

    if moving_5_files.exists():
        shutil.rmtree(moving_5_files)

    temporary_17_files = processed_results / "Registration" / "17_files"
    if temporary_17_files.exists():
        shutil.rmtree(temporary_17_files)

    print("Registered MOLLI map created.")
    print("MOLLI registration completed.\n")

#--------------------------------------------------------------------------------------------------------------------------------------------------------

def find_fixed_folder(patient_folder):
    patient_id = patient_folder.name
    registration_dir = processed_results / "Registration"

    t1w_folders = []

    for folder in patient_folder.iterdir():
        if folder.is_dir() and any(name in folder.name for name in T1W_NAMES):
            t1w_folders.append(folder)

    if "-01" in patient_id or "-03" in patient_id:
        t1w_folders = [folder for folder in t1w_folders if "BigFOV" not in folder.name]

    fixed_folder = t1w_folders[0]

    fixed_files = sorted(f for f in fixed_folder.iterdir() if f.is_file())

    fixed_input = registration_dir / "17_files"

    if fixed_input.exists():
        shutil.rmtree(fixed_input)

    if len(fixed_files) > 17:
        print(f"{patient_id} has more than 17 reference files.")

        fixed_input.mkdir(parents=True, exist_ok=True)

        for file in fixed_files[2:]:
            shutil.copy2(file, fixed_input / file.name)

        print("Using the relevant 17 files!")

    else:
        fixed_input = fixed_folder
        
    print(f"Fixed folder: {fixed_input}")
    return fixed_input

original_data = Path("../../Original data")
processed_data = Path("../../Processed - Results")

for patient_folder in sorted(original_data.iterdir()):
    if not patient_folder.is_dir():
        continue

    patient_id = patient_folder.name
    processed_results = processed_data / patient_id

    if not processed_results.exists():
        print(f"{patient_id}: no processed-results folder found, skipping.")
        continue

    t1_map = processed_results / "T1_RespTrig" / "nifti" / "t1_map_b1corr.nii.gz"
    t2_map = processed_results / "T2_RespTrig" / "nifti" / "stimfit_t2_map.nii.gz"
    molli_map = processed_results / "MOLLI" / "nifti" / "molli_t1_map_b1corr.nii.gz"

    t1_registered = processed_results / "Registration" / "T1" / "t1_map_registered.nii.gz"
    t2_registered = processed_results / "Registration" / "T2" / "stimfit_t2_map_registered.nii.gz"
    molli_registered = processed_results / "Registration" / "MOLLI" / "molli_map_registered.nii.gz"

    if not any((t1_map.exists(), t2_map.exists(), molli_map.exists())):
        print(f"{patient_id}: no quantitative maps available for registration, skipping.")
        continue

    print(f"Registering {patient_id}...")
    fixed_input = find_fixed_folder(patient_folder)
    
    if t1_map.exists() and not t1_registered.exists():
        try:
            register_t1(patient_folder)
        except Exception as e:
            log_issue(patient_id, "T1 Registering", f"{type(e).__name__}: {e}")
            
    if t2_map.exists() and not t2_registered.exists():
        try:
            register_t2(patient_folder)
        except Exception as e:
            log_issue(patient_id, "T2 Registering", f"{type(e).__name__}: {e}")
            
    if molli_map.exists() and not molli_registered.exists():
        try:
            register_molli(patient_folder)
        except Exception as e:
            log_issue(patient_id, "MOLLI Registering", f"{type(e).__name__}: {e}")
            
    print(f"{patient_id} registration check completed.")
    print("----" * 38)

print("Registration completed.")
