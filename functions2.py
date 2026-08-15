from fetch import get_t1_data
from fetch import get_t2_data
from fetch import get_b1_data
from fetch import get_b0_data

from ukat.mapping.t1 import T1
from ukat.mapping.t2_stimfit import StimFitModel
from ukat.mapping.t2_stimfit import T2StimFit
from ukat.mapping.b0 import B0

from util import resample
from ukat.utils.tools import convert_to_pi_range
import nibabel as nib
import numpy as np
import pydicom
import matplotlib.pyplot as plt
import os
import re
import SimpleITK as sitk
import sys
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
B1_NAMES = ["Cor_B1_map_BH","B1_mapping_BH","B1_MAPPING_BH"]

B0_NAMES = ["Cor_B0map_BH","B0_map_BH","B0_MAP_BH","B0 map BH"]

T1_NAMES = ["Cor_T1map_RespTrig","T1mapping_RespTrig","T1MAPPING_RESPTRIG","T1MAPPING_FB"]

T2_NAMES = ["Cor_T2map_RespTrig","T2_mapping_RespTrig","T2_MAPPING_RESPTRIG"]
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def log_issue(patient_id, section, message):
    logfile = "../../Processed - Results/Pipeline_Log.txt"
    print(f"\n[{patient_id}] {section}: {message}\n\n")

    with open(logfile, "a") as f:
        f.write(
            f"{patient_id} | "
            f"{section} | "
            f"{message}\n"
            f"{'-'*125}\n")
#--------------------------------------------------------------------------------------------------------------------------------------------------------------

def process_t1(patient_folder):
    patient_id = os.path.basename(patient_folder)
    
    t1_folders = []

    for f in os.listdir(patient_folder):
        if any(name in f for name in T1_NAMES):
            t1_folders.append(f)

    if len(t1_folders) == 0:
        if "-04" in patient_id:
            return None
        else:
            log_issue(patient_id,"T1","No T1 folder found, T1 skipped.")
            return None

    folder_name = max(t1_folders,key=lambda x: len(os.listdir(os.path.join(patient_folder, x))))
    n_files = len(os.listdir(os.path.join(patient_folder, folder_name)))
    print("T1 processing started...")
    if len(t1_folders) > 1:
        log_issue(patient_id,"T1",f"Multiple T1 folders found. Using "f"{folder_name} ({n_files} files)")
    
    folder = os.path.join(patient_folder,folder_name)
    
    first_file = os.listdir(folder)[0]
    ds = pydicom.dcmread(os.path.join(folder, first_file))

    vendor = ds.Manufacturer
    
    if vendor == "GE MEDICAL SYSTEMS" and n_files < 100:
        log_issue(patient_id,"T1",f"Skipping GE T1: only {n_files} files found")
        return None

    if vendor == "SIEMENS" and n_files < 75:
        log_issue(patient_id,"T1",f"Skipping SIEMENS T1: only {n_files} files found")
        return None
        
    series = re.sub(r'[<>:"/\\|?*]', '_', ds.SeriesDescription)
    magnitude, phase, affine, ti, tss = get_t1_data(folder, vendor, series)
    
    print("\n===== DATASET INFO =====")
    print(f"Patient ID: {patient_id}")
    print(f"Vendor: {ds.Manufacturer}")
    print(f"Scanner: {ds.ManufacturerModelName}")
    print(f"Series: {ds.SeriesDescription}")
    
    print(f"Field Strength: {ds.MagneticFieldStrength}T")
        
    print(f"Flip Angle: {ds.FlipAngle}°")
        
    print(f"Data Structure: {magnitude.shape}")    
    print("========================")

    if magnitude.shape[2] > 5:
        log_issue(patient_id,"T1",f"{magnitude.shape[2]} Slices found!")
        
    mid = magnitude.shape[2]//2
    
    phase = convert_to_pi_range(phase)
    
    #complex_data = magnitude * (np.cos(phase) + 1j * np.sin(phase)) # convert magnitude and phase into complex data
    ti = np.array(ti) * 1000  # convert TIs to ms
    tss *= 1000 # convert tss into ms
    
    mapper = T1(magnitude, ti, affine=affine, parameters=3, tss=tss, multithread=False)
    
    
    os.makedirs(f"../../Processed - Results/{patient_id}/T1_RespTrig/Plots", exist_ok=True)
    os.makedirs(f"../../Processed - Results/{patient_id}/T1_RespTrig/nifti", exist_ok=True)
    plot_folder = f"../../Processed - Results/{patient_id}/T1_RespTrig/Plots"
    nifti_folder = f"../../Processed - Results/{patient_id}/T1_RespTrig/nifti"

     # Alle nifti .nii.gz Files speichern
    mapper.to_nifti(
        output_directory=nifti_folder,
        base_file_name="",
        maps="all")
    
    #Remove _ in filename
    for filename in os.listdir(nifti_folder):
        if filename.endswith(".nii.gz") and filename!="rawdata.nii.gz":
            os.rename(os.path.join(nifti_folder, filename), os.path.join(nifti_folder, filename[1:]))
    
    # T1 Map
    plt.figure()
    plt.imshow(mapper.t1_map[:, :, mid], cmap="hot", vmin=500, vmax=2500)
    plt.colorbar(label="T1 (ms)")
    plt.title(f"T1 Map - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/T1_Map_Slice{mid+1}.png")
    plt.close()
    
    # T1 Error Map
    plt.figure()
    plt.imshow(mapper.t1_err[:, :, mid], cmap="viridis", vmin=0, vmax=500)
    plt.colorbar(label="T1 Error")
    plt.title(f"T1 Fit Error - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/T1_err_Slice{mid+1}.png")
    plt.close()
    
    # R1 Map
    r1 = mapper.r1_map() * 1000
    plt.figure()
    plt.imshow(r1[:, :, mid], cmap="hot", vmin=0, vmax=4)
    plt.colorbar(label="R1")
    plt.title(f"R1 Map - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/R^1_Slice{mid+1}.png")
    plt.close()
    
    # R^2 Map
    plt.figure()
    plt.imshow(mapper.r2[:, :, mid], cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(label="R²")
    plt.title(f"R² Fit - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/R^2_Slice{mid+1}.png")
    plt.close()
    
    # M0 Map
    plt.figure()
    plt.imshow(mapper.m0_map[:, :, mid], cmap="gray",vmin=np.percentile(mapper.m0_map[:,:,mid], 1), vmax=np.percentile(mapper.m0_map[:,:,mid], 99.88))
    plt.colorbar(label="M0")
    plt.title(f"M0 Map - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/M0_Map_Slice{mid+1}.png")
    plt.close()
    
    # M0 Error Map
    plt.figure()
    plt.imshow(mapper.m0_err[:, :, mid], cmap="viridis",vmin=np.percentile(mapper.m0_err[:,:,mid], 1), vmax=np.percentile(mapper.m0_err[:,:,mid], 99.88))
    plt.colorbar(label="M0 error")
    plt.title(f"M0 Error - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/M0_err_Slice{mid+1}.png")
    plt.close()
    
    # ==============================
    # QC MASKS
    # ==============================
    
    qc_mask = ((mapper.t1_map > 200) &(mapper.t1_map < 4000) &(mapper.r2 >= 0) &(mapper.r2 <= 1))
    
    # ==============================
    # T1 QC
    # ==============================
    
    mean_t1 = np.mean(mapper.t1_map[qc_mask])
    median_t1 = np.median(mapper.t1_map[qc_mask])
    std_t1 = np.std(mapper.t1_map[qc_mask])
    
    t1_min = np.min(mapper.t1_map[qc_mask])
    t1_max = np.max(mapper.t1_map[qc_mask])
    
    # print("===== T1 QC =====")
    
    # print(f"Mean T1: {mean_t1:.2f} ms")
    # print(f"Median T1: {median_t1:.2f} ms")
    # print(f"Std T1: {std_t1:.2f} ms")
    # print(f"T1 range: {t1_min:.2f} - {t1_max:.2f} ms")
    
    # ==============================
    # INVALID T1 VOXELS
    # ==============================
    
    invalid_t1 = np.sum((~np.isfinite(mapper.t1_map)) |(mapper.t1_map <= 0))
    total_voxels = mapper.t1_map.size
    # print(f"Invalid T1 voxels: {100*invalid_t1/total_voxels:.2f}%")
    
    # ==============================
    # R² QC
    # ==============================
    
    valid_r2 = mapper.r2[qc_mask]
    mean_r2 = np.mean(valid_r2)
    median_r2 = np.median(valid_r2)
    
    # print("\n===== R² QC =====")
    
    # print(f"Mean R²: {mean_r2:.4f}")
    # print(f"Median R²: {median_r2:.4f}")
    
    # ==============================
    # POOR FIT VOXELS
    # ==============================
    
    valid_t1 = ((mapper.t1_map > 200) &(mapper.t1_map < 4000))
    
    poor_fit = np.sum((mapper.r2 < 0.8) &valid_t1)
    
    total_valid = np.sum(valid_t1)
    
    # print(
    #     f"Poor-fit valid voxels (<0.8): "
    #     f"{100*poor_fit/total_valid:.2f}%")
    
    # ==============================
    # NUMERICALLY INVALID R²
    # ==============================
    invalid_r2 = np.sum(
        (mapper.r2 < 0) |
        (mapper.r2 > 1))
    
    # print(
    #     f"Numerically invalid R² voxels: "
    #     f"{100*invalid_r2/total_voxels:.2f}%")

    # ==============================
    # SAVE QC REPORT
    # ==============================

    qc_report = f"""
    ===== DATASET INFO =====
    
    Patient ID: {patient_id}
    Vendor: {ds.Manufacturer}
    Scanner: {ds.ManufacturerModelName}
    Series: {ds.SeriesDescription}
    
    Field Strength: {ds.MagneticFieldStrength} T
    
    Image Size: {ds.Rows} x {ds.Columns}
    Pixel Spacing: {ds.PixelSpacing}
    Slice Thickness: {ds.SliceThickness} mm
    
    TR: {ds.RepetitionTime} ms
    TE: {ds.EchoTime} ms
    Flip Angle: {ds.FlipAngle}°
    
    Images in Acquisition: {len(os.listdir(folder))}
    
    Data Structure: {magnitude.shape}
    TI Shape: {np.array(ti).shape}
    
    ===== T1 QC =====
    
    Mean T1: {mean_t1:.2f} ms
    Median T1: {median_t1:.2f} ms
    Std T1: {std_t1:.2f} ms
    T1 range: {t1_min:.2f} - {t1_max:.2f} ms
    
    Invalid T1 voxels: {100*invalid_t1/total_voxels:.2f}%
    
    ===== R² QC =====
    
    Mean R²: {mean_r2:.4f}
    Median R²: {median_r2:.4f}
    
    Poor-fit valid voxels (<0.8): {100*poor_fit/total_valid:.2f}%
    
    Numerically invalid R² voxels:
    {100*invalid_r2/total_voxels:.2f}%
    """
    
    # Save report
    with open(f"../../Processed - Results/{patient_id}/T1_RespTrig/QC_Report_{patient_id}_T1_RespTrig.txt","w") as f:
        f.write(qc_report)
    print("\nT1 processing completed and QC report saved")
    return mapper
#--------------------------------------------------------------------------------------------------------------------------------------------------------------

def process_t2(patient_folder):
    print("\nT2 processing started...\n")
    patient_id = os.path.basename(patient_folder)

    t2_folders = []
    
    for f in os.listdir(patient_folder):
        if any(name in f for name in T2_NAMES):
            t2_folders.append(f)
    
    if len(t2_folders) == 0:
        log_issue(patient_id,"T2","No T2map folder found, T2 skipped.")
        return None
    
    folder_name = max(t2_folders,key=lambda x: len(os.listdir(os.path.join(patient_folder, x))))

    folder = os.path.join(patient_folder,folder_name)
    
    if len(os.listdir(folder)) == 0:
        log_issue(patient_id,"T2","T2 folder is empty, T2 skipped.")
        return None
        
    first_file = os.listdir(folder)[0]
    ds = pydicom.dcmread(os.path.join(folder, first_file))
        
    vendor = ds.Manufacturer
    series = re.sub(r'[<>:"/\\|?*]', '_', ds.SeriesDescription)
    image, affine, te = get_t2_data(folder, vendor)  # load files and extract image, affine and echo times
    
    te = te * 1000  # convert TE to ms
    mask = image[..., 0] > 50 # 20000 # Generate a mask based on the signal intensity of the first echo
    
    print("\n===== DATASET INFO =====")
    
    print(f"Patient ID: {patient_id}")
    print(f"Vendor: {ds.Manufacturer}")
    print(f"Scanner: {ds.ManufacturerModelName}")
    print(f"Series: {ds.SeriesDescription}")
    
    print(f"Field Strength: {ds.MagneticFieldStrength} T")
    
    print(f"Flip Angle: {ds.FlipAngle}°")
        
    print(f"Data Structure: {image.shape}")
    print(f"TE Values: {te}")
    
    print("========================")

    if image.shape[2] > 5:
        log_issue(patient_id,"T2",f"{image.shape[2]} Slices found!")
    
    mid = image.shape[2]//2

    if vendor == "GE MEDICAL SYSTEMS":
        stimfit_vendor ="ge"
    elif vendor == "SIEMENS":
        stimfit_vendor = "siemens"
    else:
        raise ValueError(f"Unsupported Vendor: {vendor}")
        
    image = image.astype(np.float32)
    model = StimFitModel(ukrin_vendor = stimfit_vendor)
    mapper_2p = T2StimFit(image, affine, model, mask=mask, multithread=False)
    
    os.makedirs(f"../../Processed - Results/{patient_id}/T2_RespTrig/Plots", exist_ok=True)
    os.makedirs(f"../../Processed - Results/{patient_id}/T2_RespTrig/nifti", exist_ok=True)
    plot_folder = f"../../Processed - Results/{patient_id}/T2_RespTrig/Plots"
    nifti_folder = f"../../Processed - Results/{patient_id}/T2_RespTrig/nifti"

    # Save output maps to Nifti
    mapper_2p.to_nifti(output_directory=nifti_folder, base_file_name='stimfit', maps='all')
    
    # T2 Map
    plt.figure()
    plt.imshow(mapper_2p.t2_map[:, :, mid], cmap="hot",vmin=10, vmax=300)
    plt.colorbar(label="T2 (ms)")
    plt.title(f"T2 Map - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/T2_Map_Slice{mid+1}.png")
    plt.close()
    
    # R^2 Map
    plt.figure()
    plt.imshow(mapper_2p.r2_map[:, :, mid], cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(label="R²")
    plt.title(f"R² Fit Quality - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/R^2_Slice{mid+1}.png")
    plt.close()
    
    # M0 Map
    plt.figure()
    plt.imshow(mapper_2p.m0_map[:, :, mid], cmap="gray",vmin=np.percentile(mapper_2p.m0_map[:,:,mid], 1), vmax=np.percentile(mapper_2p.m0_map[:,:,mid], 99.88))
    plt.colorbar(label="M0")
    plt.title(f"M0 Map - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/M0_Map_Slice{mid+1}.png")
    plt.close()
    
    # ==============================
    # QC MASKS
    # ==============================
    
    qc_mask = (
        (mapper_2p.t2_map > 10) &
        (mapper_2p.t2_map < 300) &
        (mapper_2p.r2_map >= 0) &
        (mapper_2p.r2_map <= 1))
    
    # ==============================
    # T2 QC
    # ==============================
    
    mean_t2 = np.mean(mapper_2p.t2_map[qc_mask])
    median_t2 = np.median(mapper_2p.t2_map[qc_mask])
    std_t2 = np.std(mapper_2p.t2_map[qc_mask])
    
    t2_min = np.min(mapper_2p.t2_map[qc_mask])
    t2_max = np.max(mapper_2p.t2_map[qc_mask])
    
    # print("===== T2 QC =====")
    
    # print(f"Mean T2: {mean_t2:.2f} ms")
    # print(f"Median T2: {median_t2:.2f} ms")
    # print(f"Std T2: {std_t2:.2f} ms")
    # print(f"T2 range: {t2_min:.2f} - {t2_max:.2f} ms")
    
    # ==============================
    # INVALID T2 VOXELS
    # ==============================
    
    invalid_t2 = np.sum(
        (~np.isfinite(mapper_2p.t2_map)) |
        (mapper_2p.t2_map <= 0))
    
    total_voxels = mapper_2p.t2_map.size
    
    # print(f"Invalid T2 voxels: {100 * invalid_t2 / total_voxels:.2f}%")
    
    # ==============================
    # R² QC
    # ==============================
    
    valid_r2 = mapper_2p.r2_map[qc_mask]
    
    mean_r2 = np.mean(valid_r2)
    median_r2 = np.median(valid_r2)
    
    # print("\n===== R² QC =====")
    
    # print(f"Mean R²: {mean_r2:.4f}")
    # print(f"Median R²: {median_r2:.4f}")
    
    # ==============================
    # POOR-FIT VOXELS
    # ==============================
    
    valid_t2 = (
        (mapper_2p.t2_map > 10) &
        (mapper_2p.t2_map < 300)
    )
    
    poor_fit = np.sum(
        (mapper_2p.r2_map < 0.8) &
        valid_t2
    )
    
    total_valid = np.sum(valid_t2)
    
    # print(
    #     f"Poor-fit valid voxels (<0.8): "
    #     f"{100 * poor_fit / total_valid:.2f}%")

    # ==============================
    # NUMERICALLY INVALID R²
    # ==============================
    
    invalid_r2 = np.sum(
        (mapper_2p.r2_map < 0) |
        (mapper_2p.r2_map > 1)
    )
    
    # print(
    #     f"Numerically invalid R² voxels: "
    #     f"{100 * invalid_r2 / total_voxels:.2f}%")

    # ==============================
    # SAVE QC REPORT
    # ==============================
    
    qc_report = f"""
    ===== DATASET INFO =====
    
    Patient ID: {patient_id}
    Vendor: {ds.Manufacturer}
    Scanner: {ds.ManufacturerModelName}
    Series: {ds.SeriesDescription}
    
    Field Strength: {ds.MagneticFieldStrength} T
    
    Image Size: {ds.Rows} x {ds.Columns}
    Pixel Spacing: {ds.PixelSpacing}
    Slice Thickness: {ds.SliceThickness} mm
    
    TR: {ds.RepetitionTime} ms
    TE: {ds.EchoTime} ms
    Flip Angle: {ds.FlipAngle}°
    
    Images in Acquisition: {len(os.listdir(folder))}
    
    Data Structure: {image.shape}
    TE Shape: {np.array(te).shape}
    TE Values: {te}
    
    ===== T2 QC =====
    
    Mean T2: {mean_t2:.2f} ms
    Median T2: {median_t2:.2f} ms
    Std T2: {std_t2:.2f} ms
    
    T2 range: {t2_min:.2f} - {t2_max:.2f} ms
    
    Invalid T2 voxels: {100 * invalid_t2 / total_voxels:.2f}%
    
    ===== R² QC =====
    
    Mean R²: {mean_r2:.4f}
    Median R²: {median_r2:.4f}
    
    Poor-fit valid voxels (<0.8): {100 * poor_fit / total_valid:.2f}%
    
    Numerically invalid R² voxels:
    {100 * invalid_r2 / total_voxels:.2f}%
    """
    
    # Save report
    with open(
        f"../../Processed - Results/"
        f"{patient_id}/T2_RespTrig/"
        f"QC_Report_{patient_id}_T2_RespTrig.txt",
        "w"
    ) as f:
        f.write(qc_report)
    
    print("\nT2 processing completed and QC report saved\n")
    return mapper_2p
#--------------------------------------------------------------------------------------------------------------------------------------------------------------

def process_molli(patient_folder):
    patient_id = os.path.basename(patient_folder)

    new_molli_map = os.path.join("../../Processed - Results", patient_id, "MOLLI","nifti","t1_map.nii.gz")
    if os.path.exists(new_molli_map):
        print(f"New MOLLI already processed for {patient_id}, skipping...")
        return None
    
    molli = []
    for f in os.listdir(patient_folder):
        if "MOLLI" in f:
            folder = os.path.join(patient_folder, f)
            if len(os.listdir(folder)) > 5:
                molli.append(folder)
    
    if len(molli) == 0:
        if "-02" in patient_id or "_02" in patient_id:
            return None
        else:
            log_issue(patient_id,"MOLLI","No non-empty MOLLI folder found, MOLLI skipped.")
            return None
        
    molli.sort()
    first_file = os.listdir(molli[0])[0]
    ds = pydicom.dcmread(os.path.join(molli[0], first_file))
    
    print("MOLLI processing started...")

    vendor = ds.Manufacturer
    series = re.sub(r'[<>:"/\\|?*]', '_', ds.SeriesDescription)

    if vendor == "SIEMENS":
        magnitude, phase, affine, ti, tss = get_t1_data(molli[0], vendor, series)
    else:
        magnitude, phase, affine, ti, tss = get_t1_data(molli, vendor, series)
        
    
    print("\n===== DATASET INFO =====")
    print(f"Patient ID: {patient_id}")
    print(f"Vendor: {ds.Manufacturer}")
    print(f"Scanner: {ds.ManufacturerModelName}")
    print(f"Series: {ds.SeriesDescription}")
    
    print(f"Field Strength: {ds.MagneticFieldStrength}T")
        
    print(f"Flip Angle: {ds.FlipAngle}°")
        
    if vendor == "SIEMENS":
        print(f"Data Structure: {len(magnitude)} slices, {magnitude[0].shape}")
        mid = len(magnitude) // 2
    else:
        print(f"Data Structure: {magnitude.shape}")
        mid = magnitude.shape[2] // 2    
    print("========================")
        
    if vendor == "SIEMENS":  # Siemens MOLLI data have slice-specific inversion times
    
        mapper_list = []
        for i in range(len(magnitude)):
            mapper = []
            phase[i] = convert_to_pi_range(phase[i])
            #complex_data = magnitude * (np.cos(phase) + 1j * np.sin(phase)) # convert magnitude and phase into complex data
            ti[i] = np.array(ti[i]) * 1000  # convert TIs to ms
            tss *= 1000  # convert tss into ms
            if np.all(magnitude[i] == 0):
                mask = np.zeros(magnitude[i].shape[:3], dtype=np.uint8)
                mapper = T1(magnitude[i], ti[i], affine=affine, mask=mask, parameters=3, molli=True, tss=tss, multithread=False)
            else:
                mapper = T1(magnitude[i], ti[i], affine=affine, parameters=3, molli=True, tss=tss, multithread=False)
            mapper_list.append(mapper)
    
        # Append the parameters from the objects to a list
        mapper.t1_map = np.squeeze(np.moveaxis(np.array([maps.t1_map for maps in mapper_list]), 0, -1))
        mapper.t1_err = np.squeeze(np.moveaxis(np.array([maps.t1_err for maps in mapper_list]), 0, -1))
        mapper.m0_map = np.squeeze(np.moveaxis(np.array([maps.m0_map for maps in mapper_list]), 0, -1))
        mapper.m0_err = np.squeeze(np.moveaxis(np.array([maps.m0_err for maps in mapper_list]), 0, -1))
        try:
            mapper.eff_map = np.squeeze(np.moveaxis(np.array([maps.eff_map for maps in mapper_list]), 0, -1))
            mapper.eff_err = np.squeeze(np.moveaxis(np.array([maps.eff_err for maps in mapper_list]), 0, -1))
        except:
            pass
        mapper.r2 = np.squeeze(np.moveaxis(np.array([maps.r2 for maps in mapper_list]), 0, -1))
        mapper.mask = np.squeeze(np.moveaxis(np.array([maps.mask for maps in mapper_list]), 0, -1))
    
    else:
        phase = convert_to_pi_range(phase)
    
        ##complex_data = magnitude * (np.cos(phase) + 1j * np.sin(phase)) # convert magnitude and phase into complex data
        ti = np.array(ti) * 1000  # convert TIs to ms
        tss *= 1000 # convert tss into ms
    
        mapper = T1(magnitude, ti, affine=affine, parameters=3,molli=True, tss=tss, multithread=False)
    
    
    os.makedirs(f"../../Processed - Results/{patient_id}/MOLLI/Plots", exist_ok=True)
    os.makedirs(f"../../Processed - Results/{patient_id}/MOLLI/nifti", exist_ok=True)
    plot_folder = f"../../Processed - Results/{patient_id}/MOLLI/Plots"
    nifti_folder = f"../../Processed - Results/{patient_id}/MOLLI/nifti"

    # Alle nifti .nii.gz Files speichern
    mapper.to_nifti(
        output_directory=nifti_folder,
        base_file_name="",
        maps="all")
    
    #Remove _ in filename
    for filename in os.listdir(nifti_folder):
        if filename.endswith(".nii.gz") and filename!="rawdata.nii.gz":
            os.rename(os.path.join(nifti_folder, filename), os.path.join(nifti_folder, filename[1:]))
    
    # MOLLI T1 Map
    plt.figure()
    plt.imshow(mapper.t1_map[:, :, mid], cmap="hot", vmin=500, vmax=3200)
    plt.colorbar(label="MOLLI T1 (ms)")
    plt.title(f"MOLLI T1 Map - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/MOLLI_T1_Map_Slice{mid+1}.png")
    plt.close()
    
    # MOLLI T1 Error Map
    plt.figure()
    plt.imshow(mapper.t1_err[:, :, mid], cmap="viridis", vmin=0, vmax=500)
    plt.colorbar(label="MOLLI T1 Error")
    plt.title(f"MOLLI T1 Fit Error - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/MOLLI_T1_err_Slice{mid+1}.png")
    plt.close()
    
    # R1 Map
    r1 = mapper.r1_map() * 1000
    plt.figure()
    plt.imshow(r1[:, :, mid], cmap="hot", vmin=0, vmax=4)
    plt.colorbar(label="R1")
    plt.title(f"R1 Map - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/R^1_Slice{mid+1}.png")
    plt.close()
    
    # R^2 Map
    plt.figure()
    plt.imshow(mapper.r2[:, :, mid], cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(label="R²")
    plt.title(f"R² Fit - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/R^2_Slice{mid+1}.png")
    plt.close()
    
    # M0 Map
    plt.figure()
    plt.imshow(mapper.m0_map[:, :, mid], cmap="gray",vmin=np.percentile(mapper.m0_map[:,:,mid], 1), vmax=np.percentile(mapper.m0_map[:,:,mid], 99.88))
    plt.colorbar(label="M0")
    plt.title(f"M0 Map - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/M0_Map_Slice{mid+1}.png")
    plt.close()
    
    # M0 Error Map
    plt.figure()
    plt.imshow(mapper.m0_err[:, :, mid], cmap="viridis",vmin=np.percentile(mapper.m0_err[:,:,mid], 1), vmax=np.percentile(mapper.m0_err[:,:,mid], 99.88))
    plt.colorbar(label="M0 error")
    plt.title(f"M0 Error - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/M0_err_Slice{mid+1}.png")
    plt.close()
        
    # ==============================
    # QC MASKS
    # ==============================
    
    qc_mask = ((mapper.t1_map > 200) &(mapper.t1_map < 4000) &(mapper.r2 >= 0) &(mapper.r2 <= 1))
    
    # ==============================
    # MOLLI T1 QC
    # ==============================
    
    mean_t1 = np.mean(mapper.t1_map[qc_mask])
    median_t1 = np.median(mapper.t1_map[qc_mask])
    std_t1 = np.std(mapper.t1_map[qc_mask])
    
    t1_min = np.min(mapper.t1_map[qc_mask])
    t1_max = np.max(mapper.t1_map[qc_mask])
    
    # print("===== MOLLI T1 QC =====")
    
    # print(f"Mean MOLLI T1: {mean_t1:.2f} ms")
    # print(f"Median MOLLI T1: {median_t1:.2f} ms")
    # print(f"Std MOLLI T1: {std_t1:.2f} ms")
    # print(f"MOLLI T1 range: {t1_min:.2f} - {t1_max:.2f} ms")
    
    # ==============================
    # INVALID MOLLI T1 VOXELS
    # ==============================
    
    invalid_t1 = np.sum((~np.isfinite(mapper.t1_map)) |(mapper.t1_map <= 0))
    total_voxels = mapper.t1_map.size
    # print(f"Invalid MOLLI T1 voxels: {100*invalid_t1/total_voxels:.2f}%")
    
    # ==============================
    # R² QC
    # ==============================
    
    valid_r2 = mapper.r2[qc_mask]
    mean_r2 = np.mean(valid_r2)
    median_r2 = np.median(valid_r2)
    
    # print("\n===== R² QC =====")
    
    # print(f"Mean R²: {mean_r2:.4f}")
    # print(f"Median R²: {median_r2:.4f}")
    
    # ==============================
    # POOR FIT VOXELS
    # ==============================
    
    valid_t1 = ((mapper.t1_map > 200) &(mapper.t1_map < 4000))
    
    poor_fit = np.sum((mapper.r2 < 0.8) &valid_t1)
    
    total_valid = np.sum(valid_t1)
    
    # print(
    #     f"Poor-fit valid voxels (<0.8): "
    #     f"{100*poor_fit/total_valid:.2f}%")
    
    # ==============================
    # NUMERICALLY INVALID R²
    # ==============================
    invalid_r2 = np.sum(
        (mapper.r2 < 0) |
        (mapper.r2 > 1))
    
    # print(
    #     f"Numerically invalid R² voxels: "
    #     f"{100*invalid_r2/total_voxels:.2f}%")

    # ==============================
    # SAVE QC REPORT
    # ==============================

    if vendor == "SIEMENS":
        data_structure = f"{len(magnitude)} slices, {magnitude[0].shape}"
    else:
        data_structure = str(magnitude.shape)
    
    qc_report = f"""
    ===== DATASET INFO =====
    
    Patient ID: {patient_id}
    Vendor: {ds.Manufacturer}
    Scanner: {ds.ManufacturerModelName}
    Series: {ds.SeriesDescription}
    
    Field Strength: {ds.MagneticFieldStrength} T
    
    Image Size: {ds.Rows} x {ds.Columns}
    Pixel Spacing: {ds.PixelSpacing}
    Slice Thickness: {ds.SliceThickness} mm
    
    TR: {ds.RepetitionTime} ms
    TE: {ds.EchoTime} ms
    Flip Angle: {ds.FlipAngle}°
    
    Images in Acquisition: {len(os.listdir(folder))}
    
    Data Structure: {data_structure}
    TI Shape: {np.array(ti).shape}
    
    ===== MOLLI T1 QC =====
    
    Mean MOLLI T1: {mean_t1:.2f} ms
    Median MOLLI T1: {median_t1:.2f} ms
    Std MOLLI T1: {std_t1:.2f} ms
    MOLLI T1 range: {t1_min:.2f} - {t1_max:.2f} ms
    
    Invalid MOLLI T1 voxels: {100*invalid_t1/total_voxels:.2f}%
    
    ===== R² QC =====
    
    Mean R²: {mean_r2:.4f}
    Median R²: {median_r2:.4f}
    
    Poor-fit valid voxels (<0.8): {100*poor_fit/total_valid:.2f}%
    
    Numerically invalid R² voxels:
    {100*invalid_r2/total_voxels:.2f}%
    """
    
    # Save report
    with open(f"../../Processed - Results/{patient_id}/MOLLI/QC_Report_{patient_id}_MOLLI.txt","w") as f:
        f.write(qc_report)
    print()
    print("\nMOLLI processing completed and QC report saved")
    return mapper

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def process_b0(patient_folder):
    patient_id = os.path.basename(patient_folder)
    
    b0_folders = []
    for f in os.listdir(patient_folder):
        if any(name in f for name in B0_NAMES):
            b0_folders.append(os.path.join(patient_folder, f))

    if len(b0_folders) == 0:
        log_issue(patient_id,"B0","No B0 folder found, B0 skipped.")
        return None

    # ONE FOLDER (GE)
    if len(b0_folders) == 1:
        file_folder = b0_folders[0]
        first_file = os.listdir(file_folder)[0]
        ds = pydicom.dcmread(os.path.join(file_folder, first_file),stop_before_pixels=True)

    # TWO FOLDERS
    elif len(b0_folders) == 2:
        folder_info = []
        for folder in b0_folders:
            first_file = os.listdir(folder)[0]
            ds = pydicom.dcmread(os.path.join(folder, first_file),stop_before_pixels=True)

            folder_info.append((folder, ds))

        vendor = folder_info[0][1].Manufacturer

        # Siemens
        if vendor == "SIEMENS":
            magnitude_folder = None
            phase_folder = None
            for folder, ds in folder_info:
                image_type = ds.ImageType
                if "P" in image_type:
                    phase_folder = folder
                elif "M" in image_type:
                    magnitude_folder = folder

            if magnitude_folder is None or phase_folder is None:
                log_issue(patient_id,"B0","Could not identify Siemens magnitude/phase folders.")
                return None

            file_folder = [magnitude_folder,phase_folder]
            print(
            f"Siemens B0 folders identified:\n"
            f"  Magnitude: {os.path.basename(magnitude_folder)}\n"
            f"  Phase: {os.path.basename(phase_folder)}")
            ds = pydicom.dcmread(os.path.join(magnitude_folder,os.listdir(magnitude_folder)[0]),stop_before_pixels=True)

        # GE duplicate acquisition
        else:
            candidates = []
            for folder, ds in folder_info:
                n_files = len(os.listdir(folder))
                candidates.append((n_files,ds.SeriesNumber,folder))

            candidates.sort(key=lambda x: (x[0], x[1]),reverse=True)
            best_files, best_series, best_folder = candidates[0]
            file_folder = best_folder

            log_issue(patient_id,"B0",
                f"Multiple B0 folders found. Using "
                f"{os.path.basename(best_folder)} "
                f"({best_files} files, Series {best_series})")

            ds = pydicom.dcmread(os.path.join(best_folder,os.listdir(best_folder)[0]),stop_before_pixels=True)

    else:
        log_issue(patient_id,"B0",f"Unexpected number of B0 folders: {len(b0_folders)}")
        return None

    # LOAD B0 DATA
    vendor = ds.Manufacturer
    magnitude, phase, affine, te = get_b0_data(file_folder,vendor)

    te = te * 1000

    mapper = B0(phase,te,affine=affine,unwrap=True)

    # SAVE RESULTS
    result_folder = (f"../../Processed - Results/{patient_id}/B0")

    os.makedirs(result_folder, exist_ok=True)

    mapper.to_nifti(output_directory=result_folder,base_file_name="")

    # remove leading "_"
    for filename in os.listdir(result_folder):
        if (
            filename.endswith(".nii.gz")
            and filename != "rawdata_magnitude.nii.gz"
            and filename != "rawdata_phase.nii.gz"):

            os.rename(os.path.join(result_folder, filename),os.path.join(result_folder, filename[1:]))

    mid_slice = mapper.b0_map.shape[2] // 2
    plt.figure()
    plt.imshow(mapper.b0_map[:, :, mid_slice],cmap="jet")
    plt.colorbar(label="B0 (Hz)")
    plt.title(f"B0 Map - Slice {mid_slice + 1}")
    plt.tight_layout()
    plt.savefig(os.path.join(result_folder,f"B0_Map_Slice{mid_slice + 1}.png"),dpi=300)
    plt.close()

    # print("\n===== B0 QC =====\n")

    # print(f"Mean B0: {np.mean(mapper.b0_map):.3f}")
    # print(f"Std B0: {np.std(mapper.b0_map):.3f}")
    # print(f"Min B0: {np.min(mapper.b0_map):.3f}")
    # print(f"Max B0: {np.max(mapper.b0_map):.3f}")

    # QC
    qc_file = os.path.join(result_folder,f"B0_QC_{patient_id}.txt")

    with open(qc_file, "w") as f:
        f.write("===== B0 QC =====\n\n")

        f.write(
            f"Mean B0: "
            f"{np.mean(mapper.b0_map):.3f}\n"
        )

        f.write(
            f"Std B0: "
            f"{np.std(mapper.b0_map):.3f}\n"
        )

        f.write(
            f"Min B0: "
            f"{np.min(mapper.b0_map):.3f}\n"
        )

        f.write(
            f"Max B0: "
            f"{np.max(mapper.b0_map):.3f}\n"
        )

    print("\nB0 processing completed and QC saved")

    return mapper
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def process_b1(patient_folder):
    patient_id = os.path.basename(patient_folder)

    b1_folders = []
    for f in os.listdir(patient_folder):
        if any(name in f for name in B1_NAMES):
            b1_folders.append(os.path.join(patient_folder, f))

    if len(b1_folders) == 0:
        log_issue(os.path.basename(patient_folder),"B1","No B1 folder found, B1 skipped.")
        return None

    # GE (1 folder)
    if len(b1_folders) == 1:
        file_folder = b1_folders[0]
        first_file = os.listdir(file_folder)[0]
        ds = pydicom.dcmread(os.path.join(file_folder, first_file),stop_before_pixels=True)

    elif len(b1_folders) == 2:
        folder_info = []
        
        for folder in b1_folders:
            first_file = os.listdir(folder)[0]
            ds = pydicom.dcmread(os.path.join(folder, first_file),stop_before_pixels=True)
            folder_info.append((folder, ds))
    
        vendor = folder_info[0][1].Manufacturer
    
        # SIEMENS: magnitude + flip angle map
        if vendor == "SIEMENS":
            magnitude_folder = None
            flipangle_folder = None
    
            for folder, ds in folder_info:
                image_type = ds.ImageType
    
                if "FLIP ANGLE MAP" in image_type:
                    flipangle_folder = folder
                elif "M" in image_type:
                    magnitude_folder = folder
    
            if magnitude_folder is None or flipangle_folder is None:
                log_issue(patient_id,"B1","Could not identify Siemens magnitude/flip-angle folders.")
                return None
    
            file_folder = [magnitude_folder, flipangle_folder]
    
            ds = pydicom.dcmread(os.path.join(magnitude_folder,os.listdir(magnitude_folder)[0]),stop_before_pixels=True)
    
        # GE duplicate acquisition
        else:
            candidates = []
            for folder, ds in folder_info:
                n_files = len(os.listdir(folder))
                candidates.append((n_files, ds.SeriesNumber, folder))
    
            candidates.sort(key=lambda x: (x[0], x[1]),reverse=True)
            best_files, best_series, best_folder = candidates[0]
            file_folder = best_folder
            log_issue(
                patient_id,
                "B1",
                f"Multiple B1 folders found. Using "
                f"{os.path.basename(best_folder)} "
                f"({best_files} files, Series {best_series})"
            )
    
            ds = pydicom.dcmread(os.path.join(best_folder,os.listdir(best_folder)[0]),stop_before_pixels=True)
    else:
        log_issue(patient_id,"B1",f"Unexpected number of B1 folders: {len(b1_folders)}")
        return None

    vendor = ds.Manufacturer
    magnitude, flipanglemap, affine, flip_nom = get_b1_data(file_folder,vendor)
    
    # CREATE OUTPUT FOLDER
    result_folder = (f"../../Processed - Results/{patient_id}/B1")
    os.makedirs(result_folder, exist_ok=True)

    # SAVE MAGNITUDE
    magnitude_img = nib.Nifti1Image(magnitude,affine=affine)
    nib.save(magnitude_img,os.path.join(result_folder,"b1_magnitude.nii.gz"))

    # SAVE FLIP ANGLE MAP
    flipanglemap_img = nib.Nifti1Image(flipanglemap,affine=affine)
    nib.save(flipanglemap_img,os.path.join(result_folder,"b1_flipanglemap.nii.gz"))
    
    # B1 FACTOR MAP
    if vendor == "SIEMENS":
        factormap = flipanglemap / (flip_nom * 10)
    
    else:  # GE
        flip_nom = ds[0x0018, 0x1314].value * 10
        factormap = flipanglemap / flip_nom
        
    factormap_img = nib.Nifti1Image(factormap,affine=affine)
    nib.save(factormap_img,os.path.join(result_folder,"b1_factormap.nii.gz"))

    print("\nB1 Calculation done and maps saved\n")

    return factormap

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
def apply_b1_to_t1(patient_folder):
    patient_id = os.path.basename(patient_folder)
  
    plot_folder = (f"../../Processed - Results/{patient_id}/T1_RespTrig/Plots")
    nifti_folder = (f"../../Processed - Results/{patient_id}/T1_RespTrig/nifti")

    # Load saved T1 map as SITK image
    t1_img = sitk.ReadImage(
        os.path.join(nifti_folder, "t1_map.nii.gz"),
        sitk.sitkFloat32)
    
    # Load saved B1 factor map
    b1_img = sitk.ReadImage(
        f"../../Processed - Results/{patient_id}/B1/b1_factormap.nii.gz",
        sitk.sitkFloat32)
    
    # Resample B1 to T1 geometry
    b1_img = resample(t1_img, b1_img)
    
    # Convert to arrays
    b1_factor_resampled = sitk.GetArrayFromImage(b1_img)
    t1_array = sitk.GetArrayFromImage(t1_img)
    
    # Apply correction
    t1_corr = t1_array * b1_factor_resampled
    t1_diff = t1_corr - t1_array

    t1_nifti = nib.load(os.path.join(nifti_folder, "t1_map.nii.gz"))  
    t1_corr_img = nib.Nifti1Image(t1_corr.astype(np.float32),affine=t1_nifti.affine)
    
    nib.save(t1_corr_img,os.path.join(nifti_folder,"t1_map_b1corr.nii.gz"))

    mid = t1_corr.shape[0] // 2

    plt.figure()
    plt.imshow(t1_diff[mid, :, :],cmap="bwr",vmin=-500,vmax=500)
    plt.colorbar(label="ΔT1 (ms)")
    plt.savefig(f"{plot_folder}/T1_B1corr_Difference_Slice{mid+1}.png")
    plt.close()

    # print("\n===== B1 CORRECTION FOR T1 QC =====\n")

    # print(f"Mean B1 factor: {np.mean(b1_factor_resampled):.3f}")
    # print(f"Std B1 factor: {np.std(b1_factor_resampled):.3f}")
    
    # print(f"Mean T1 before: {np.mean(t1_array):.2f} ms")
    # print(f"Mean T1 after : {np.mean(t1_corr):.2f} ms")
    print("B1 Correction for T1 done and corrected map saved")
    return t1_corr

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

def apply_b1_to_molli(patient_folder):
    patient_id = os.path.basename(patient_folder)

    plot_folder = f"../../Processed - Results/{patient_id}/MOLLI/Plots"
    nifti_folder = f"../../Processed - Results/{patient_id}/MOLLI/nifti"

    # Load saved MOLLI T1 map
    molli_img = sitk.ReadImage(
        os.path.join(nifti_folder, "t1_map.nii.gz"),
        sitk.sitkFloat32)

    # Load saved B1 factor map
    b1_img = sitk.ReadImage(
        f"../../Processed - Results/{patient_id}/B1/b1_factormap.nii.gz",
        sitk.sitkFloat32)

    # Resample B1 to MOLLI geometry
    b1_img = resample(molli_img, b1_img)

    # Convert to arrays
    b1_factor_resampled = sitk.GetArrayFromImage(b1_img)
    molli_array = sitk.GetArrayFromImage(molli_img)

    # Apply correction
    molli_corr = molli_array * b1_factor_resampled

    molli_diff = molli_corr - molli_array

    # Save corrected map
    molli_corr_img = sitk.GetImageFromArray(molli_corr.astype(np.float32))

    molli_corr_img.CopyInformation(molli_img)

    sitk.WriteImage(molli_corr_img,os.path.join(nifti_folder, "molli_t1_map_b1corr.nii.gz"))

    # Middle slice for QC plot
    mid = molli_corr.shape[0] // 2

    plt.figure()
    plt.imshow(molli_diff[mid, :, :],cmap="bwr",vmin=-500,vmax=500)
    plt.colorbar(label="ΔT1 (ms)")
    plt.title(f"MOLLI T1 Correction Difference - Slice {mid+1}")
    plt.savefig(f"{plot_folder}/MOLLI_B1corr_Difference_Slice{mid+1}.png")
    plt.close()

    # print("\n===== B1 CORRECTION FOR MOLLI QC =====\n")

    # print(f"Mean MOLLI T1 before: {np.mean(molli_array):.2f} ms")
    # print(f"Mean MOLLI T1 after : {np.mean(molli_corr):.2f} ms")
    print("B1 Correction for MOLLI done and corrected map saved")

    return molli_corr

#--------------------------------------------------------------------------------------------------------------------------------------------------------------

def process_patient(patient_folder):
    patient_id = os.path.basename(patient_folder)
    result_folder = f"../../Processed - Results/{patient_id}"
    
    if all([os.path.exists(f"{result_folder}/T1_RespTrig"),
           os.path.exists(f"{result_folder}/T2_RespTrig"),
           os.path.exists(f"{result_folder}/MOLLI"),
           os.path.exists(f"{result_folder}/B0"),
           os.path.exists(f"{result_folder}/B1")]):
        return None

    print(f"Started processing {patient_id}\n")

    if not os.path.exists(f"{result_folder}/B1"):
        try:
            b1_factor = process_b1(patient_folder)
        except Exception as e:
            log_issue(patient_id,"B1",f"{type(e).__name__}: {e}")

    if not os.path.exists(f"{result_folder}/T1_RespTrig"):
        try:
            t1_mapper = process_t1(patient_folder)
            if t1_mapper is not None:
                apply_b1_to_t1(patient_folder)
        except Exception as e:
            log_issue(patient_id,"T1",f"{type(e).__name__}: {e}")

    if not os.path.exists(f"{result_folder}/T2_RespTrig"):
        try:
            process_t2(patient_folder)
        except Exception as e:
            log_issue(patient_id,"T2",f"{type(e).__name__}: {e}")

    if not os.path.exists(f"{result_folder}/MOLLI"):
        try:
            molli_mapper = process_molli(patient_folder)
            if molli_mapper is not None:
                apply_b1_to_molli(patient_folder)
        except Exception as e:
            log_issue(patient_id,"MOLLI",f"{type(e).__name__}: {e}")

    if not os.path.exists(f"{result_folder}/B0"):
        try:
            process_b0(patient_folder)
        except Exception as e:
            log_issue(patient_id,"B0",f"{type(e).__name__}: {e}")

    # if b1_factor is not None:
    #     save_b1_corr_qc(patient_id,b1_factor,t1_mapper,t1_corr,t2_mapper,t2_corr,molli_mapper,molli_corr)

    print(f"\n{patient_id} fully processed")
    print("----" * 37)
    print()

    return
