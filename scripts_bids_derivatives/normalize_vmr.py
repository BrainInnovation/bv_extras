"""
Anatomical Normalization for BrainVoyager Data (BIDS Structure)

This script performs anatomical normalization on skull-stripped VMR files
for multiple subjects and sessions in a BIDS-like directory structure.

What it does:
- Copies a preprocessed skull-stripped VMR and its associated .v16 file
- Applies Talairach and/or MNI transformations using BrainVoyager, based on the
  configuration flags set below
- Skips transformation steps if the expected output file already exists and skipping is enabled
- Saves output files in a dedicated normalization workflow folder
- Logs every step and captures any errors for later review

Inputs:
- Preprocessed, skull-stripped anatomical VMR files in:
    derivatives/
        └── workflow_id-1_type-2_name-anat-preprocessing/
            └── sub-01/ses-01/anat/sub-01_ses-01_T1w_IIHC.vmr

Outputs:
- Normalized VMR files saved in:
    derivatives/
        └── workflow_id-2_type-3_name-anat-normalization/
            └── sub-01/ses-01/anat/
                ├── sub-01_ses-01_T1w_IIHC.vmr
                ├── sub-01_ses-01_T1w_IIHC_MNI.vmr (if MNI)
                ├── sub-01_ses-01_T1w_IIHC_aTal.vmr (if Talairach)
                └── log file

Configuration:
- Adapt the configuration section below to match your project structure
- Choose whether to apply MNI and/or Talairach normalization by setting:
      MNI = True / False
      Talairach = True / False
- Set SKIP_IF_OUTPUT_EXISTS = True to avoid reprocessing if output already exists
"""

import os
from datetime import datetime
from shutil import copyfile

# === CONFIGURATION ===
project_path = '/Users/judithe/Documents/BrainVoyager/Projects/NeWBI4fMRI2020/derivatives/'
subjects = [1,2,3,4,5,6,7,8,9,10,11,12]
sessions = [1]
anat_prep = 'workflow_id-1_type-2_name-anat-preprocessing'
anat_norm = 'workflow_id-2_type-3_name-anat-normalization'

MNI = True           # Apply MNI normalization
Talairach = True    # Apply Talairach transformation
SKIP_IF_OUTPUT_EXISTS = True  # Skip if output files already exist
suffix = '_IIHC'     # VMR should be at least skull-stripped

# === HELPER FUNCTIONS ===
def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def construct_vmr_filename(sub_id, ses_id, suffix):
    return f'{sub_id}_{ses_id}_T1w{suffix}.vmr'

def log_message(log_file, message):
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(f'{message}\n')

# === MAIN LOOP ===
errors = []
bv.close_all()

for subj in subjects:
    sub_id = f'sub-{subj:02d}'
    for sess in sessions:
        ses_id = f'ses-{sess:02d}'
        now_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        norm_path = os.path.join(project_path, anat_norm, sub_id, ses_id, 'anat')
        ensure_dir(norm_path)
        log_file = os.path.join(norm_path, f'{sub_id}_{ses_id}_anat-normalization_{timestamp}.log')

        print(f'=== Normalizing {sub_id} {ses_id} at {now_str} ===')
        log_message(log_file, f"=== Starting normalization for {sub_id} {ses_id} with suffix '{suffix}' ===")
        log_message(log_file, f"Start time: {now_str}\n")

        filename = construct_vmr_filename(sub_id, ses_id, suffix)
        input_vmr = os.path.join(project_path, anat_prep, sub_id, ses_id, 'anat', filename)
        input_v16 = input_vmr[:-4] + '.v16'
        norm_vmr = os.path.join(norm_path, filename)
        norm_v16 = norm_vmr[:-4] + '.v16'

        # Step 0: Copy files
        try:
            copyfile(input_vmr, norm_vmr)
            copyfile(input_v16, norm_v16)
            log_message(log_file, f'Copied VMR from {input_vmr} to {norm_vmr}\n')
        except Exception as e:
            msg = f'{sub_id} {ses_id}: Failed to copy input files: {e}'
            log_message(log_file, msg)
            errors.append(msg)
            continue

        # Step 1: Transformations
        try:
            brainvoyager.open(norm_vmr)

            if Talairach:
                tal_file = norm_vmr[:-4] + "_aTal.vmr"
                if SKIP_IF_OUTPUT_EXISTS and os.path.exists(tal_file):
                    msg = f"Skipped Talairach transformation — output already exists: {tal_file}"
                    log_message(log_file, msg)
                    print(f'{sub_id} {ses_id}: {msg}')
                else:
                    bv.active_document.auto_acpc_tal_transformation()
                    if os.path.exists(tal_file):
                        log_message(log_file, f'Finished Talairach transformation → {tal_file}')
                    else:
                        msg = f'{sub_id} {ses_id}: Talairach transformation failed — expected file not created: {tal_file}'
                        log_message(log_file, msg)
                        errors.append(msg)

            if MNI:
                mni_file = norm_vmr[:-4] + "_MNI.vmr"
                if SKIP_IF_OUTPUT_EXISTS and os.path.exists(mni_file):
                    msg = f"Skipped MNI normalization — output already exists: {mni_file}"
                    log_message(log_file, msg)
                    print(f'{sub_id} {ses_id}: {msg}')
                else:
                    bv.active_document.normalize_to_mni_space()
                    if os.path.exists(mni_file):
                        log_message(log_file, f'Finished MNI transformation → {mni_file}')
                    else:
                        msg = f'{sub_id} {ses_id}: MNI normalization failed — expected file not created: {mni_file}'
                        log_message(log_file, msg)
                        errors.append(msg)

            bv.close_all()

        except Exception as e:
            msg = f'{sub_id} {ses_id}: Normalization error: {e}'
            log_message(log_file, msg)
            errors.append(msg)
            continue

        end_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_message(log_file, f'\n=== Finished normalization for {sub_id} {ses_id} at {end_time} ===\n')
        print(f'=== Finished {sub_id} {ses_id} ===\nLog saved to: {log_file}\n')

# === SUMMARY ===

if errors:
    print('\n=== SUMMARY OF ERRORS ===')
    for err in errors:
        print(f'- {err}')
else:
    print('\nAll normalizations completed without errors.')
