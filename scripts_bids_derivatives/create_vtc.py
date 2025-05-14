"""
VTC Creation Script for BrainVoyager BIDS data — Supports MNI, ACPC, Talairach, and Native Space

This script creates VTC files from preprocessed and coregistered functional data
in BrainVoyager, based on the corresponding intra-session VMRs. The resulting VTCs can be generated in
Native, ACPC, Talairach, or MNI space, depending on user configuration.

The script is structured for multi-subject, multi-session, and multi-run pipelines
using a BIDS-like directory layout. It locates FMR and transformation
files depending on the configuration section in the script, applies optional 
bounding box constraints in Native or ACPC space, and logs each action for 
troubleshooting purposes.

Key Features:
- Loads the relevant FMR and transformation files
- Supports intra-session coregistration or run-specific alignment
- Optional bounding box definition for native or ACPC VTCs
- Configurable resolution and interpolation settings
- Skips VTC creation if the output already exists (optional)
- Logs every step and tracks errors for review

Example Inputs:
- Anatomical VMR: derivatives/workflow_id-1_type-2_name-anat-preprocessing/sub-01/ses-01/anat/sub-01_ses-01_T1w_IIHC.vmr
- Functional FMR: derivatives/workflow_id-3_type-1_name-func-preprocessing/sub-01/ses-01/func/sub-01_ses-01_task-Localizer_run-02_bold_3DMCTS_SCCTBL.fmr
- Transformation files: *_IA.trf, *_FA.trf, *_ACPC.trf, *_aACPC.tal, *_TO_MNI_a12_adjBBX.trf

Example Outputs:
- Subject and Session-Specific Logfile: derivatives/workflow_id-5_type-5_name-func-normalization/sub-01/ses-01/func/sub-01_ses-01_task-Localizer_vtc_creation_20250513_181945.log
- Native VTC: derivatives/workflow_id-5_type-5_name-func-normalization/sub-01/ses-01/func/sub-01_ses-01_task-Localizer_run-02_bold_3DMCTS_SCCTBL_256_sinc3_2x1_NATIVE.vtc
- ACPC VTC:  ..._ACPC.vtc
- TAL VTC:   ..._TAL.vtc
- MNI VTC:   ..._MNI.vtc

Configuration:
- You must adapt the configuration section at the top of this script to match:
    - Your project folder structure
    - Workflow folder names (anat/func/coregistration)
    - Filename suffixes 
    - List of subject IDs, session IDs, and runs
    - settings for the VTC creation

"""

import os
from datetime import datetime

# === CONFIGURATION ===

SKIP_IF_VTC_EXISTS = True  # Skip VTC creation if output already exists

project_path = '/Users/judithe/Documents/BrainVoyager/Projects/NeWBI4fMRI2020/derivatives/'
anat_prep = 'workflow_id-1_type-2_name-anat-preprocessing'
anat_norm = 'workflow_id-2_type-3_name-anat-normalization'
func_prep = 'workflow_id-3_type-1_name-func-preprocessing'
coreg = 'workflow_id-4_type-4_name-func-to-anat-coreg'
vtc_folder = 'workflow_id-5_type-5_name-func-normalization'

subjects = [1,2,3,4,5,6,7,8,9,10]
sessions = [1]
runs = [1,2]
task_name = 'task-Localizer'

use_intrasesscoreg = True
coreg_run = 1
coreg_task_name = 'task-Localizer'

res = 2
interpol = 2  # 0=nearest, 1=trilin, 2=sinc

fmr_suffix = '_3DMCTS_SCCTBL.fmr'
anat_suffix = '_IIHC.vmr'

coreg_BBR = True
fa_suffix = '_BBR_FA.trf' if coreg_BBR else '_FA.trf'
ia_suffix = '_IA.trf'

use_bbx = False
bbx = dict(from_x=1, to_x=174, from_y=1, to_y=245, from_z=1, to_z=196)
intensity_thresh = 100

CREATE_MNI_VTC = True
CREATE_ACPC_VTC = False
CREATE_TAL_VTC = False
extended_TAL_space = True
CREATE_NATIVE_VTC = False


# === INTERNAL SETTINGS ===
interp_map = {0: 'nearest', 1: 'trilin', 2: 'sinc3'}
interp_label = interp_map.get(interpol, 'invalid')
if interp_label == 'invalid':
    raise ValueError(f"Invalid interpolation value: {interpol}")
    
def apply_bounding_box(doc):
    if use_bbx:
        doc.vtc_creation_use_bounding_box = True
        doc.vtc_creation_bounding_box_from_x = bbx['from_x']
        doc.vtc_creation_bounding_box_to_x = bbx['to_x']
        doc.vtc_creation_bounding_box_from_y = bbx['from_y']
        doc.vtc_creation_bounding_box_to_y = bbx['to_y']
        doc.vtc_creation_bounding_box_from_z = bbx['from_z']
        doc.vtc_creation_bounding_box_to_z = bbx['to_z']
    else:
        doc.vtc_creation_use_bounding_box = False

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def pjoin(*args):
    return os.path.join(*args)

def log_message(log_file, msg):
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(f'{msg}\n')

# === MAIN ===
acpc_trf_suffixes = ["_aACPC.trf", "_ACPC.trf"]
tal_file_suffixes = ["_aACPC.tal", "_ACPC.tal"]
errors = []
bv.close_all()

for subj in subjects:
    sub_id = f'sub-{subj:02d}'
    for sess in sessions:
        ses_id = f'ses-{sess:02d}'
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        out_path = pjoin(project_path, vtc_folder, sub_id, ses_id, 'func')
        ensure_dir(out_path)
        log_file = pjoin(out_path, f'{sub_id}_{ses_id}_{task_name}_vtc_creation_{timestamp}.log')

        log_message(log_file, f'VTC Creation Log — {sub_id} {ses_id} {task_name}')
        log_message(log_file, f'Intra-Session Coregistration Used: {use_intrasesscoreg}')
        if use_intrasesscoreg:
            log_message(log_file, f'Intra-Session Coregistration Task Name: {coreg_task_name} |  Intra-Session Coregistration Run-ID: {coreg_run}')
        log_message(log_file, f'Resolution: {res} | Interpolation: {interp_label} | Intensity Threshold: {intensity_thresh} | Bounding Box: {use_bbx}')
        if use_bbx:
            log_message(log_file, f'Bounding Box for Native and/or ACPC VTCs: {bbx}')
        log_message(log_file, f'Spaces — NATIVE: {CREATE_NATIVE_VTC}, ACPC: {CREATE_ACPC_VTC}, TALAIRACH: {CREATE_TAL_VTC}, MNI: {CREATE_MNI_VTC}\n')

        try:
            vmr_file = f'{sub_id}_{ses_id}_T1w{anat_suffix}'
            vmr_path = pjoin(project_path, anat_prep, sub_id, ses_id, 'anat', vmr_file)

        except Exception as e:
            msg = f'ERROR: {sub_id} {ses_id}: Failed to locate VMR: {e}'
            print(msg)
            errors.append(msg)
            log_message(log_file, msg)
            continue

        for run in runs:
            run_id = f'run-{run:02d}'
            coreg_id = f'run-{coreg_run:02d}' if use_intrasesscoreg else run_id
            coreg_task = coreg_task_name if use_intrasesscoreg else task_name

            # Locate required files
            base = f'{sub_id}_{ses_id}_{task_name}_{run_id}_bold'
            fmr = pjoin(project_path, func_prep, sub_id, ses_id, 'func', f'{base}{fmr_suffix}')
            coreg_file = f'{sub_id}_{ses_id}_{coreg_task}_{coreg_id}_bold{fmr_suffix}'
            ia_trf = pjoin(project_path, coreg, sub_id, ses_id, 'func', f'{coreg_file[:-4]}-TO-{vmr_file[:-4]}{ia_suffix}')
            fa_trf = pjoin(project_path, coreg, sub_id, ses_id, 'func', f'{coreg_file[:-4]}-TO-{vmr_file[:-4]}{fa_suffix}')
            mni_trf = pjoin(project_path, anat_norm, sub_id, ses_id, 'anat', f'{vmr_file[:-4]}_TO_MNI_a12_adjBBX.trf')
            acpc_trf = next((pjoin(project_path, anat_norm, sub_id, ses_id, 'anat', f'{vmr_file[:-4]}{suffix}')
                             for suffix in acpc_trf_suffixes if os.path.exists(pjoin(project_path, anat_norm, sub_id, ses_id, 'anat', f'{vmr_file[:-4]}{suffix}'))), None)
            tal_file = next((pjoin(project_path, anat_norm, sub_id, ses_id, 'anat', f'{vmr_file[:-4]}{suffix}')
                             for suffix in tal_file_suffixes if os.path.exists(pjoin(project_path, anat_norm, sub_id, ses_id, 'anat', f'{vmr_file[:-4]}{suffix}'))), None)


            if CREATE_NATIVE_VTC:
                missing = [path for path in [fmr, ia_trf, fa_trf] if not os.path.exists(path)]
                if missing:
                    msg = f'Skipped Native VTC for {base}{fmr_suffix}→ missing: {", ".join(os.path.basename(p) for p in missing)}'
                    print(msg)
                    errors.append(msg)
                    log_message(log_file, msg)
                else:
                    doc_vmr = brainvoyager.open(vmr_path)
                    apply_bounding_box(doc_vmr)
                    vtc = pjoin(out_path, f'{base}{fmr_suffix[:-4]}_{doc_vmr.dim_z}_{interp_label}_{res}x{doc_vmr.voxelsize_x}_NATIVE.vtc')
                    if SKIP_IF_VTC_EXISTS and os.path.exists(vtc): 
                        msg = f'Skipped Native VTC → already exists: {vtc}'
                        print(msg)
                        errors.append(msg)
                        log_message(log_file, msg)
                    else:
                        doc_vmr.create_vtc_in_native_space(fmr, ia_trf, fa_trf, vtc, res, interpol, intensity_thresh)
                        log_message(log_file, f'Native VTC → {vtc}')
                    bv.close_all()

            if CREATE_ACPC_VTC:
                missing = [path for path in [fmr, ia_trf, fa_trf] if not os.path.exists(path)]
                if not acpc_trf: 
                    missing.append("ACPC_TRF")
                if missing:
                    msg = f'Skipped ACPC VTC for {base}{fmr_suffix} → missing: {", ".join(os.path.basename(p) if isinstance(p, str) else p for p in missing)}'
                    print(msg)
                    errors.append(msg)
                    log_message(log_file, msg)
                else:
                    doc_vmr = brainvoyager.open(vmr_path)
                    apply_bounding_box(doc_vmr)
                    vtc = pjoin(out_path, f'{base}{fmr_suffix[:-4]}_{doc_vmr.dim_z}_{interp_label}_{res}x{doc_vmr.voxelsize_x}_ACPC.vtc')
                    if SKIP_IF_VTC_EXISTS and os.path.exists(vtc): 
                        msg = f'Skipped ACPC VTC → already exists: {vtc}'
                        print(msg)
                        errors.append(msg)
                        log_message(log_file, msg)
                    else: 
                        doc_vmr.create_vtc_in_acpc_space(fmr, ia_trf, fa_trf, acpc_trf, vtc, res, interpol, intensity_thresh)
                        log_message(log_file, f'ACPC VTC → {vtc}')
                    bv.close_all()

            if CREATE_TAL_VTC:
                missing = [path for path in [fmr, ia_trf, fa_trf] if not os.path.exists(path)]
                if not acpc_trf: 
                    missing.append("ACPC_TRF")
                if not tal_file: 
                    missing.append("TAL_FILE")
                if missing:
                    msg = f'Skipped Talairach VTC for {base}{fmr_suffix} → missing: {", ".join(os.path.basename(p) if isinstance(p, str) else p for p in missing)}'
                    print(msg)
                    errors.append(msg)
                    log_message(log_file, msg)                    
                else:
                    doc_vmr = brainvoyager.open(vmr_path)
                    vtc = pjoin(out_path, f'{base}{fmr_suffix[:-4]}_{doc_vmr.dim_z}_{interp_label}_{res}x{doc_vmr.voxelsize_x}_TAL.vtc')
                    doc_vmr.vtc_creation_extended_tal_space = extended_TAL_space
                    if SKIP_IF_VTC_EXISTS and os.path.exists(vtc): 
                        msg = f'Skipped Talairach VTC → already exists: {vtc}'
                        print(msg)
                        errors.append(msg)
                        log_message(log_file, msg)
                    else: 
                        doc_vmr.create_vtc_in_tal_space(fmr, ia_trf, fa_trf, acpc_trf, tal_file, vtc, res, interpol, intensity_thresh)
                        log_message(log_file, f'Talairach VTC → {vtc}')
                    bv.close_all()

            if CREATE_MNI_VTC:
                missing = [path for path in [fmr, ia_trf, fa_trf, mni_trf] if not os.path.exists(path)]
                if missing:
                    msg = f'Skipped MNI VTC for {base}{fmr_suffix} → missing: {", ".join(os.path.basename(p) for p in missing)}'
                    print(msg)
                    errors.append(msg)
                    log_message(log_file, msg)
                else:
                    doc_vmr = brainvoyager.open(vmr_path)
                    vtc = pjoin(out_path, f'{base}{fmr_suffix[:-4]}_{doc_vmr.dim_z}_{interp_label}_{res}x{doc_vmr.voxelsize_x}_MNI.vtc')
                    if SKIP_IF_VTC_EXISTS and os.path.exists(vtc): 
                        msg = f'Skipped MNI VTC → already exists: {vtc}'
                        print(msg)
                        errors.append(msg)
                        log_message(log_file, msg)
                    else: 
                        doc_vmr.create_vtc_in_mni_space(fmr, ia_trf, fa_trf, mni_trf, vtc, res, interpol, intensity_thresh)
                        log_message(log_file, f'MNI VTC → {vtc}')
                    bv.close_all()

# === SUMMARY ===
if errors:
    print('\n=== SUMMARY OF ERRORS ===')
    for err in errors:
        print(f'- {err}')
else:
    print('\nAll VTC creations completed without errors.')
