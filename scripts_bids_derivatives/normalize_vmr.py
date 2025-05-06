"""
Anatomical Normalization for BrainVoyager Data (BIDS Structure)
Requires skull-stripped VMR
"""

import os
from datetime import datetime
from shutil import copyfile

# === CONFIGURATION ===
project_path = '/Users/judithe/Documents/BrainVoyager/Projects/NeWBI4fMRI2020/derivatives/'
subjects = [1,2,3,4,5,6,7,8,9,10,11,12]
sessions = [1]
anat_prep = 'workflow_id-1_type-2_name-anat-preprocessing/'
anat_norm = 'workflow_id-2_type-3_name-anat-normalization/'
MNI = True # set to True for MNI transformation
Talairach = False # set to True for Talairach transformation
suffix = '_IIHC'   # or _SAG_ISO-1.0_IIHC (depending on the preprocessing of the VMR that has been applied, minimum required is a skull-stripped VMR)

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
        norm_vmr = os.path.join(norm_path, filename)

        # Step 0: Copy files
        copyfile(input_vmr, norm_vmr)
        copyfile(input_vmr[:-4] + '.v16', norm_vmr[:-4] + '.v16')
        log_message(log_file, f'Copied VMR from {input_vmr} to {norm_vmr}\n')

        # Step 1: Transformations
        if Talairach:
            brainvoyager.open(norm_vmr)
            bv.active_document.auto_acpc_tal_transformation()
            log_message(log_file, f'Finished Talairach transformation → {bv.active_document.path_file_name}')
            bv.close_all()

        if MNI:
            brainvoyager.open(norm_vmr)
            bv.active_document.normalize_to_mni_space()
            log_message(log_file, f'Finished MNI transformation → {bv.active_document.path_file_name}')
            bv.close_all()

        end_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_message(log_file, f'\n=== Finished normalization for {sub_id} {ses_id} at {end_time} ===\n')
        print(f'=== Finished {sub_id} {ses_id} ===\nLog saved to: {log_file}\n')
 