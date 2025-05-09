"""
Anatomical Preprocessing for BrainVoyager Data Analysis File and Folder Structure (BIDS)
Updated 2025-05-09

This script performs modular anatomical preprocessing of T1-weighted VMR files
using BrainVoyager in a BIDS-like folder structure.

What it does:
- Converts raw NIfTI or VMR into a working VMR file
- Applies selected preprocessing steps (sagittal alignment, iso-voxel transformation, inhomogeneity correction)
- Saves each intermediate step in a new file
- Skips steps if output already exists (optional)
- Logs all actions, paths, and errors

Inputs:
- Raw data in NIfTI or BrainVoyager VMR format
- Stored in: rawdata/sub-XX/ses-YY/anat or derivatives/rawdata_bv/sub-XX/ses-YY/anat

Outputs:
- Preprocessed VMRs in: derivatives/workflow_id-1_type-2_name-anat-preprocessing/sub-XX/ses-YY/anat
- One log file per subject/session

Configuration:
- Use the configuration block below to define subjects, sessions, enabled steps, the preprocessing options and the preprocessing order
- Control output skipping with SKIP_IF_OUTPUT_EXISTS
"""

import os
from datetime import datetime
from shutil import copyfile

# === CONFIGURATION ===
project_path = '/Users/judithe/Documents/BrainVoyager/Projects/NeWBI4fMRI2020'
subjects = [1,2,3,4,5,6,7,8,9,10,11,12]
sessions = [1]
anat_prep = 'derivatives/workflow_id-1_type-2_name-anat-preprocessing/'

# Optional: define suffix of preprocessed file to resume from a preprocessing step (e.g. "_SAG_ISO-1.0") or "" for raw VMR
starting_suffix = ""

SKIP_IF_OUTPUT_EXISTS = True  # Set to True to skip steps if output already exists

# Preprocessing options
options = {
    "transform_sag": False,
    
    "transform_iso_voxel": False,
    "target_res": 1.0,
    "framing_cube_dim": 256,
    "interpolation_method": 3,
    
    "correct_intensity_inhomogeneities": True,
    "include_brain_extraction": True,
    "n_cycles": 3,
    "tissue_range_thresh": 0.25,
    "intensity_thresh": 0.3,
    "fit_polynom_order": 3,
}

# Pipeline execution order
pipeline_order = [
    "transform_sag",
    "transform_iso_voxel",
    "correct_intensity_inhomogeneities",
]

# === HELPER FUNCTIONS ===
def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def construct_vmr_filename(sub_id, ses_id, suffix=''):
    return f'{sub_id}_{ses_id}_T1w{suffix}.vmr'

def get_raw_nii_path(sub_id, ses_id):
    return os.path.join(project_path, 'rawdata', sub_id, ses_id, 'anat',
                        f'{sub_id}_{ses_id}_T1w.nii.gz')

def get_anat_prep_path(sub_id, ses_id):
    return os.path.join(project_path, anat_prep, sub_id, ses_id, 'anat')

def get_log_path(sub_id, ses_id):
    prep_path = get_anat_prep_path(sub_id, ses_id)
    ensure_dir(prep_path)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join(prep_path, f'{sub_id}_{ses_id}_preprocessing_{timestamp}.log')

def log_message(log_file, message):
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(f'{message}\n')

# === MAIN LOOP ===
errors = []

for subj in subjects:
    sub_id = f'sub-{subj:02d}'
    bv.close_all()
    for sess in sessions:
        ses_id = f'ses-{sess:02d}'
        log_file = get_log_path(sub_id, ses_id)
        log_message(log_file, f"=== Starting preprocessing for {sub_id} {ses_id} with starting_suffix '{starting_suffix}' ===")
        log_message(log_file, f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log_message(log_file, "Preprocessing options:")
        for key, value in options.items():
            log_message(log_file, f"  {key}: {value}")
        log_message(log_file, "")

        print(f'=== Starting preprocessing for {sub_id} {ses_id} ===')

        raw_nii_path = get_raw_nii_path(sub_id, ses_id)
        prep_path = get_anat_prep_path(sub_id, ses_id)
        ensure_dir(prep_path)

        raw_vmr_path = os.path.join(project_path, 'derivatives/rawdata_bv', sub_id, ses_id, 'anat',
                                    construct_vmr_filename(sub_id, ses_id))
        out_path = os.path.join(prep_path, construct_vmr_filename(sub_id, ses_id, starting_suffix))
        current_file = out_path

        # Step 0: Load or verify starting VMR
        if starting_suffix:
            if os.path.exists(out_path):
                log_message(log_file, f'Using existing preprocessed VMR → {out_path}')
            else:
                msg = f'ERROR: Starting VMR file not found: {out_path}. Skipping.'
                log_message(log_file, msg)
                print(msg)
                errors.append(msg)
                continue
        else:
            if os.path.exists(out_path):
                log_message(log_file, f'Using existing raw VMR in preprocessing folder → {out_path}')
            elif os.path.exists(raw_vmr_path):
                doc_vmr = brainvoyager.open(raw_vmr_path)
                doc_vmr.save_as(out_path)
                copyfile(raw_vmr_path[:-4] + '.v16', out_path[:-4] + '.v16')
                doc_vmr.close()
                log_message(log_file, f'Copied VMR from rawdata_bv → {out_path}')
            else:
                doc_vmr = brainvoyager.open(raw_nii_path)
                doc_vmr.save_as(out_path)
                doc_vmr.close()
                log_message(log_file, f'Converted NIfTI to VMR → {out_path}')

        # Step 1: Run preprocessing pipeline
        for step in pipeline_order:
            if not options.get(step, False):
                continue

            try:
                doc_vmr = brainvoyager.open(current_file)
                log_message(log_file, f'Starting step: {step} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

                if step == "transform_sag":
                    next_file = current_file[:-4] + "_SAG.vmr"
                    if SKIP_IF_OUTPUT_EXISTS and os.path.exists(next_file):
                        msg = f"Skipped {step} — output already exists: {next_file}"
                        log_message(log_file, msg)
                        print(f'{sub_id} {ses_id}: {msg}')
                        current_file = next_file
                        continue
                    doc_vmr.transform_to_std_sag(next_file)
                    current_file = next_file

                elif step == "transform_iso_voxel":
                    next_file = current_file[:-4] + f"_ISO-{options['target_res']}.vmr"
                    if SKIP_IF_OUTPUT_EXISTS and os.path.exists(next_file):
                        msg = f"Skipped {step} — output already exists: {next_file}"
                        log_message(log_file, msg)
                        print(f'{sub_id} {ses_id}: {msg}')
                        current_file = next_file
                        continue
                    doc_vmr.transform_to_isovoxel(options["target_res"], options["framing_cube_dim"], options["interpolation_method"], next_file)
                    current_file = next_file

                elif step == "correct_intensity_inhomogeneities":
                    next_file = current_file[:-4] + "_IIHC.vmr"
                    if SKIP_IF_OUTPUT_EXISTS and os.path.exists(next_file):
                        msg = f"Skipped {step} — output already exists: {next_file}"
                        log_message(log_file, msg)
                        print(f'{sub_id} {ses_id}: {msg}')
                        current_file = next_file
                        continue
                    
                    doc_vmr.correct_intensity_inhomogeneities_ext(
                        options["include_brain_extraction"],
                        options["n_cycles"],
                        options["tissue_range_thresh"],
                        options["intensity_thresh"],
                        options["fit_polynom_order"]
                    )
                    doc_vmr.save_as(doc_vmr.path_file_name)
                    current_file = current_file = next_file

                doc_vmr.close()
                log_message(log_file, f'Finished step: {step} → {current_file}')

            except Exception as e:
                msg = f'{sub_id} {ses_id}: ERROR in step "{step}": {e}'
                log_message(log_file, msg)
                print(msg)
                errors.append(msg)

        log_message(log_file, f'=== Finished preprocessing for {sub_id} {ses_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n')
        print(f'=== Finished preprocessing for {sub_id} {ses_id} ===\nLog saved to: {log_file}\n')

# === FINAL SUMMARY ===
if errors:
    print('\n=== SUMMARY OF ERRORS ===')
    for err in errors:
        print(f'- {err}')
else:
    print('\nAll preprocessing steps completed without errors.')
