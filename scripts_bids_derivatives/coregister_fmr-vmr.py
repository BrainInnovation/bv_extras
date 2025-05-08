"""
Functional-Anatomical Coregistration for BrainVoyager Data (BIDS Structure)

This script performs functional-to-anatomical coregistration using BrainVoyager.
It is designed for use with BIDS-like projects, where both anatomical (VMR) and
functional (FMR) data have been preprocessed and stored under a structured
`derivatives/` folder.

What it does:
- Loads a preprocessed VMR file for each subject/session.
- Loads preprocessed FMR files for specified runs.
- Performs coregistration (IA + FA transformation) from each FMR to the VMR.
- Uses either gradient-based or boundary-based registration (BBR).
- Saves all coregistration result files to a dedicated output workflow folder.
- Moves intermediate files (e.g. *_ETC-*, BBR files) to the same folder.
- Creates a log file per subject/session with timestamps and error messages.

Inputs:
- Preprocessed anatomical (VMR) and functional (FMR) files located in:
    derivatives/
        ├── workflow_id-1_type-2_name-anat-preprocessing/
        │   └── sub-03/ses-01/anat/sub-03_ses-01_T1w_IIHC.vmr
        ├── workflow_id-3_type-1_name-func-preprocessing/
        │   └── sub-03/ses-01/func/sub-03_ses-01_task-Localizer_run-01_bold_3DMCTS_SCCTBL.fmr

Output:
- Coregistration files (IA.trf and FA.trf) for each FMR are saved to:
    derivatives/
        └── workflow_id-4_type-4_name-func-to-anat-coreg/
            └── sub-03/ses-01/func/
                ├── sub-03_ses-01_task-Localizer_run-01_bold_..._IA.trf
                ├── sub-03_ses-01_task-Localizer_run-01_bold_..._FA.trf
                └── coregistration log file

Configuration:
- You must adapt the configuration section at the top of this script to match:
    - Your project folder structure
    - Workflow folder names (anat/func/coregistration)
    - Filename suffixes for VMR and FMR files
    - List of subject IDs, session IDs, and runs
    - Whether to use BBR or gradient-based registration

Dependencies:
- BrainVoyager Python API (accessible via `brainvoyager.open`)
- Standard Python 3 libraries only
"""

import os
from datetime import datetime

# === CONFIGURATION ===

PROJECT_PATH = '/Users/judithe/Documents/BrainVoyager/Projects/NeWBI4fMRI2020/derivatives/'
ANAT_PREP = 'workflow_id-1_type-2_name-anat-preprocessing'
FUNC_PREP = 'workflow_id-3_type-1_name-func-preprocessing'
COREG_OUT = 'workflow_id-4_type-4_name-func-to-anat-coreg'

ANAT_SUFFIX = '_IIHC'
FUNC_SUFFIX = '_3DMCTS_SCCTBL'

SUBJECTS = [1,2,3,4,5,6,7,8,9,10,11,12]
SESSIONS = [1]
RUNS = [1, 2]
TASK_NAME = 'task-Localizer'
USE_BBR = True


# === HELPER FUNCTIONS ===

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def construct_filename(sub_id, ses_id, task_id, suffix, run_id=None):
    if run_id:
        return f'{sub_id}_{ses_id}_{task_id}_{run_id}_bold{suffix}.fmr'
    else:
        return f'{sub_id}_{ses_id}_T1w{suffix}.vmr'

def log_message(log_file, message):
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(f'{message}\n')

def move_file(src, dst, run_id, log_file, errors):
    try:
        os.replace(src, dst)
    except Exception as e:
        msg = f'{run_id}: Failed to move file {os.path.basename(src)}: {e}'
        log_message(log_file, msg)
        errors.append(msg)

def run_coregistration(input_vmr, input_fmr, run_id, doc_vmr, fmr_name, coreg_path, fmr_path, log_file, errors):
    if USE_BBR:
        doc_vmr.coregister_fmr_to_vmr_using_bbr(input_fmr)
        ia_name = fmr_name[:-4] + "-TO-" + doc_vmr.file_name[:-4] + "_IA.trf"
        fa_name = fmr_name[:-4] + "-TO-" + doc_vmr.file_name[:-4] + "_BBR_FA.trf"
    else:
        doc_vmr.coregister_fmr_to_vmr(input_fmr, True)
        ia_name = fmr_name[:-4] + "-TO-" + doc_vmr.file_name[:-4] + "_IA.trf"
        fa_name = fmr_name[:-4] + "-TO-" + doc_vmr.file_name[:-4] + "_FA.trf"

    # Move transformation files
    move_file(os.path.join(fmr_path, ia_name), os.path.join(coreg_path, ia_name), run_id, log_file, errors)
    move_file(os.path.join(fmr_path, fa_name), os.path.join(coreg_path, fa_name), run_id, log_file, errors)

    if USE_BBR:
        for f in os.listdir(fmr_path):
            if 'bbr' in f.lower():
                move_file(os.path.join(fmr_path, f), os.path.join(coreg_path, f), run_id, log_file, errors)

    return ia_name, fa_name


# === MAIN SCRIPT ===

errors = []
bv.close_all()

for subj in SUBJECTS:
    sub_id = f'sub-{subj:02d}'
    for sess in SESSIONS:
        ses_id = f'ses-{sess:02d}'
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        # Prepare paths
        coreg_path = os.path.join(PROJECT_PATH, COREG_OUT, sub_id, ses_id, 'func')
        anat_path = os.path.join(PROJECT_PATH, ANAT_PREP, sub_id, ses_id, 'anat')
        ensure_dir(coreg_path)

        # Log file
        log_file = os.path.join(coreg_path, f'{sub_id}_{ses_id}_coregistration_{timestamp}.log')
        log_message(log_file, f"=== Starting coregistration for {sub_id} {ses_id} {TASK_NAME} ({'BBR' if USE_BBR else 'gradient'}) ===")
        log_message(log_file, f"Start time: {log_time}\n")

        # Input VMR
        vmr_file = construct_filename(sub_id, ses_id, '', ANAT_SUFFIX)
        input_vmr = os.path.join(anat_path, vmr_file)

        if not os.path.exists(input_vmr):
            msg = f'VMR file not found: {input_vmr}. Skipping {sub_id} {ses_id}.'
            log_message(log_file, msg)
            errors.append(msg)
            continue

        for run in RUNS:
            run_id = f'run-{run:02d}'
            func_path = os.path.join(PROJECT_PATH, FUNC_PREP, sub_id, ses_id, 'func')
            fmr_file = construct_filename(sub_id, ses_id, TASK_NAME, FUNC_SUFFIX, run_id)
            input_fmr = os.path.join(func_path, fmr_file)

            if not os.path.exists(input_fmr):
                msg = f'{run_id}: FMR file not found: {input_fmr}. Skipping.'
                log_message(log_file, msg)
                errors.append(msg)
                continue

            try:
                doc_vmr = brainvoyager.open(input_vmr)
                log_message(log_file, f'{run_id}: Starting at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

                ia, fa = run_coregistration(input_vmr, input_fmr, run_id, doc_vmr, fmr_file, coreg_path, func_path, log_file, errors)

                log_message(log_file,
                    f"{run_id}: Coregistration finished\n"
                    f"FMR → {input_fmr}\n"
                    f"IA  → {os.path.join(coreg_path, ia)}\n"
                    f"FA  → {os.path.join(coreg_path, fa)}\n"
                )
                print(f'Coregistration complete: {input_fmr}')
                bv.close_all()

            except Exception as e:
                msg = f'{run_id}: Coregistration error: {e}'
                log_message(log_file, msg)
                errors.append(msg)

        # Move ETC files if BBR
        if USE_BBR:
            for f in os.listdir(anat_path):
                if '_ETC-' in f:
                    move_file(os.path.join(anat_path, f), os.path.join(coreg_path, f), 'ETC', log_file, errors)

        log_message(log_file, f'=== Finished preprocessing for {sub_id} {ses_id} ===\n')
        print(f'=== Done with {sub_id} {ses_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n')

# === SUMMARY ===

if errors:
    print('\n=== SUMMARY OF ERRORS ===')
    for err in errors:
        print(f'- {err}')
else:
    print('\nAll coregistrations completed without errors.')
