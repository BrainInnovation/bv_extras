"""
fMRI Preprocessing Pipeline for BrainVoyager Data Analysis File and Folder Structure (BIDS)

This script performs preprocessing on fMRI data in a BIDS-like directory structure
using BrainVoyager. It supports multiple subjects, sessions, and runs, with flexible
control over preprocessing steps and logging.

What it does:
- Converts raw NIfTI files to FMR format or resumes preprocessing from an existing FMR with a specified suffix
- Applies preprocessing steps such as:
    - Motion correction (intrasession or single-run)
    - Slice scan time correction (SSC)
    - Temporal high-pass filtering
    - Temporal smoothing
    - Spatial smoothing
    - Mean intensity adjustment
- Supports configuration of step order, reference run/volume for motion correction
- Skips conversion steps if the expected preprocessed file already exists
- Saves outputs in a dedicated preprocessing workflow folder
- Logs every step of the process and captures any errors for later review

Inputs:
- Raw NIfTI files located in:
        rawdata/
            └── sub-XX/ses-YY/func/
                ├── sub-XX_ses-YY_task-XY_run-ZZ_bold.nii.gz
    OR raw FMR files located in
        derivatives/
            └── rawdata_bv/
                └── sub-XX/ses-YY/func/
                    ├── sub-XX_ses-YY_task-XY_run-ZZ_bold.fmr
    OR partly preprocessed FMR files located in
        derivatives/
            └── workflow_id-3_type-1_name-func-preprocessing/
                └── sub-XX/ses-YY/func/
                    ├── sub-XX_ses-YY_task-XY_run-ZZ_bold_*.fmr

Outputs:
- Preprocessed FMR files saved in:
    derivatives/
        └── workflow_id-3_type-1_name-func-preprocessing/
            └── sub-XX/ses-YY/func/
                ├── sub-XX_ses-YY_task-XY_run-ZZ_bold_*.fmr
                └── log file (timestamped)

Configuration:
- Adjust the `options` dictionary to enable/disable preprocessing steps
- Modify `pipeline_order` to control the order of applied steps
- Use `starting_suffix` to resume from an intermediate preprocessing file
- All options and paths can be customized to suit your project setup
"""

import os
from datetime import datetime

# === CONFIGURATION ===
project_path = '/Users/judithe/Documents/BrainVoyager/Projects/NeWBI4fMRI2020'
subjects = range(1,13)
sessions = [1]
runs = [1,2]
task_name = 'task-Localizer'
func_prep = 'derivatives/workflow_id-3_type-1_name-func-preprocessing/'

# Optional: define suffix of preprocessed file to resume from a preprocessing step (e.g. "_3DMCTS_SCCTBL") or "" for raw FMR
starting_suffix = ""

# Preprocessing options
options = {
    "adjust_mean_intensity": False,
    
    "motion_correction": True,
    "intrasession": True,
    "moco_ref_volume": "first", # "first" or "last"
    "moco_ref_run": "run-01",
    "moco_ref_task": "task-Localizer",

    "slice_timing": True,
    "ssc_interpolation_method": 1, # 0 = linear, 1 = cubic_spline, 2 = sinc

    "highpass_filter": True,
    "n_cycles": 3,

    "smooth_temporal": False,
    "temp_gauss_fwhm": 2,
    "temp_fwhm_unit": "dp", # "secs" or "dp"   

    "smooth_spatial": False,
    "gauss_fwhm": 4,
    "fwhm_unit": "mm"
}

# Full pipeline, change order to change the order of the preprocessing
pipeline_order = [
    "adjust_mean_intensity",
    "motion_correction",
    "slice_timing",
    "highpass_filter",
    "smooth_temporal",
    "smooth_spatial"
]

# === ERROR COLLECTION ===
all_errors = []

# === HELPER FUNCTIONS ===
def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def construct_fmr_filename(sub_id, ses_id, run_id, suffix=''):
    return f'{sub_id}_{ses_id}_{task_name}_{run_id}_bold{suffix}.fmr'

def get_raw_nii_path(sub_id, ses_id, run_id):
    return os.path.join(project_path, 'rawdata', sub_id, ses_id, 'func',
                        f'{sub_id}_{ses_id}_{task_name}_{run_id}_bold.nii.gz')

def get_func_prep_path(sub_id, ses_id):
    return os.path.join(project_path, func_prep, sub_id, ses_id, 'func')

def get_ref_run_path(sub_id, ses_id):
    return os.path.join(project_path, 'derivatives/rawdata_bv', sub_id, ses_id, 'func',
                        f'{sub_id}_{ses_id}_{options["moco_ref_task"]}_{options["moco_ref_run"]}_bold.fmr')

def get_log_path(sub_id, ses_id):
    prep_path = get_func_prep_path(sub_id, ses_id)
    ensure_dir(prep_path)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join(prep_path, f'{sub_id}_{ses_id}_preprocessing_{timestamp}.log')

def log_message(log_file, message):
    if 'ERROR' in message:
        all_errors.append(message)
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(f'{message}\n')

# === MAIN LOOP ===
for subj in subjects:
    sub_id = f'sub-{subj:02d}'
    for sess in sessions:
        ses_id = f'ses-{sess:02d}'
        log_file = get_log_path(sub_id, ses_id)
        log_message(log_file, f"=== Starting preprocessing for {sub_id} {ses_id} with starting_suffix '{starting_suffix}' ===")
        log_message(log_file, f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log_message(log_file, "Preprocessing options:")

        for key, value in options.items():
            log_message(log_file, f"  {key}: {value}")

        log_message(log_file, "")  # Blank line for readability
        
        print(f'=== Starting preprocessing for {sub_id} {ses_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n')

        # Reference volume setup
        if options["motion_correction"] and options["intrasession"]:
            ref_run = get_ref_run_path(sub_id, ses_id)
            
            # check that the reference FMR exists
            if not os.path.exists(ref_run):
                raw_nii_ref = os.path.join(
                    project_path, 'rawdata', sub_id, ses_id, 'func',
                    f'{sub_id}_{ses_id}_{options["moco_ref_task"]}_{options["moco_ref_run"]}_bold.nii.gz'
                    )

                if os.path.exists(raw_nii_ref):
                    doc_tmp = brainvoyager.open(raw_nii_ref)
                    doc_tmp.close()
                    log_message(log_file,
                        f'Reference run missing – converted NIfTI to FMR → {ref_run}')
                else:
                    msg = (f'ERROR: Reference FMR ({ref_run}) and NIfTI '
                           f'({raw_nii_ref}) not found.  ➜ Skipping {sub_id}, {ses_id}.')
                    log_message(log_file, msg)
                    print(msg)
                    continue                      # skip this session loop

            ref_vol = 0
            if options["moco_ref_volume"] == "last":
                doc_ref = brainvoyager.open(ref_run)
                ref_vol = doc_ref.n_volumes - 1
                doc_ref.close()

        for run in runs:
            run_id = f'run-{run:02d}'
            raw_nii_path = get_raw_nii_path(sub_id, ses_id, run_id)
            prep_path = get_func_prep_path(sub_id, ses_id)
            ensure_dir(prep_path)

            raw_fmr_path = os.path.join(project_path, 'derivatives/rawdata_bv', sub_id, ses_id, 'func',
                                        construct_fmr_filename(sub_id, ses_id, run_id))
            out_path = os.path.join(prep_path, construct_fmr_filename(sub_id, ses_id, run_id, starting_suffix))
            current_file = out_path

            # Step 0: Load or verify starting FMR file
            if starting_suffix:
                # If a preprocessed FMR is expected, check if it exists    
                if os.path.exists(out_path):
                    log_message(log_file, f'Run {task_name}_{run_id}: Using existing preprocessed FMR → {out_path}')
                else:
                    message = f'ERROR: Starting FMR file not found: {out_path}. Skipping this run.'
                    log_message(log_file, message)
                    print(message)
                    continue  # Skip to the next run
            else:
                # Otherwise start from raw NIfTI or from raw FMR in rawdata_bv or in preprocessing workflow folder
                if os.path.exists(out_path):
                    log_message(log_file, f'Run {task_name}_{run_id}: Using existing raw FMR in preprocessing folder → {out_path}')
                elif os.path.exists(raw_fmr_path):
                    doc_fmr = brainvoyager.open(raw_fmr_path)
                    doc_fmr.save_as(out_path)
                    doc_fmr.close()
                    log_message(log_file, f'Run {task_name}_{run_id}: Copied FMR from rawdata_bv → {out_path}')
                elif os.path.exists(raw_nii_path):
                    doc_fmr = brainvoyager.open(raw_nii_path)
                    doc_fmr.save_as(out_path)
                    doc_fmr.close()
                    log_message(log_file, f'Run {task_name}_{run_id}: Converted NIfTI to FMR → {out_path}')
                else:
                    message = (
                        f'ERROR: Could not find any valid input file for {sub_id}_{ses_id}_{task_name}_{run_id}: '
                        f'Checked:\n'
                        f'  - Preprocessed FMR: {out_path}\n'
                        f'  - Raw FMR: {raw_fmr_path}\n'
                        f'  - Raw NIfTI: {raw_nii_path}\n'
                        f'Skipping this run.'
                    )
                    log_message(log_file, message)
                    print(message)
                    continue  # Skip to the next run


            # Step 1: Run enabled steps in pipeline_order
            for step in pipeline_order:
                if not options.get(step, False):
                    continue  # Only run if enabled in options

                try:
                    doc_fmr = brainvoyager.open(current_file)
                    log_message(log_file, f'{task_name}_{run_id}: Starting step: {step} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

                    if step == "adjust_mean_intensity":
                        doc_fmr.adjust_mean_intensity()

                    elif step == "motion_correction":
                        if options["intrasession"]:
                            doc_fmr.correct_motion_to_run_ext(ref_run, ref_vol, 2, 1, 100, 1, 1)
                        else:
                            moco_ref_vol = 0 if options["moco_ref_volume"] == "first" else doc_fmr.n_volumes - 1
                            doc_fmr.correct_motion_ext(moco_ref_vol, 2, 1, 100, 1, 1)

                    elif step == "slice_timing":
                        doc_fmr.correct_slicetiming_using_timingtable(options["ssc_interpolation_method"])

                    elif step == "highpass_filter":
                        doc_fmr.filter_temporal_highpass_glm_fourier(options["n_cycles"])

                    elif step == "smooth_temporal":
                        doc_fmr.smooth_temporal(options["temp_gauss_fwhm"], options["temp_fwhm_unit"])

                    elif step == "smooth_spatial":
                        doc_fmr.smooth_spatial(options["gauss_fwhm"], options["fwhm_unit"])

                    current_file = doc_fmr.preprocessed_fmr_name
                    doc_fmr.close()

                    log_message(log_file, f'{task_name}_{run_id}: Finished step: {step} → {current_file}')
                    
                except Exception as e:
                    log_message(log_file, f'ERROR: {sub_id}_{ses_id}_{task_name}_{run_id}: ERROR in step "{step}": {str(e)}')
                    
            log_message(log_file, f'=== Finished preprocessing for {sub_id}_{ses_id}_{task_name}_{run_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n')
                    
        log_message(log_file, f'=== Finished preprocessing for {sub_id} {ses_id} ===\n')
        print(f'=== Finished preprocessing for {sub_id} {ses_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n\n')
        print(f'Log file saved: {log_file}')

# === FINAL ERROR SUMMARY ===
print("\n=== SUMMARY OF ERRORS ===")
if all_errors:
    for error in all_errors:
        print(error)
else:
    print("No errors encountered during preprocessing.")
