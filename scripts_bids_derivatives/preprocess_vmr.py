"""
Anatomical Preprocessing for BrainVoyager Data Analysis File and Folder Structure (BIDS)
Created on 2025-05-06
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

# Preprocessing options
options = {
    "transform_sag": False,
    
    "transform_iso_voxel": False,
    "target_res": 1.0, # specified in mm
    "framing_cube_dim": 256, # or 384 or 512 for sub-millimeter data
    "interpolation_method": 3, #1 = trilinear, 2 = cubic spline, 3 = sinc

    "correct_intensity_inhomogeneities": True,
    "include_brain_extraction": True, 
    "n_cycles": 3,
    "tissue_range_thresh": 0.25,
    "intensity_thresh": 0.3,
    "fit_polynom_order": 3,
}

# Full pipeline, change order to change the order of the preprocessing
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
        
        raw_nii_path = get_raw_nii_path(sub_id, ses_id)
        prep_path = get_anat_prep_path(sub_id, ses_id)
        ensure_dir(prep_path)

        raw_vmr_path = os.path.join(project_path, 'derivatives/rawdata_bv', sub_id, ses_id, 'anat',
                                    construct_vmr_filename(sub_id, ses_id))
        out_path = os.path.join(prep_path, construct_vmr_filename(sub_id, ses_id, starting_suffix))
        current_file = out_path

        
        # Step 0: Load or verify starting VMR file
        if starting_suffix:
            # If a preprocessed FMR is expected, check if it exists    
            if os.path.exists(out_path):
                log_message(log_file, f'Using existing preprocessed VMR → {out_path}')
            else:
                message = f'ERROR: Starting VMR file not found: {out_path}. Skipping this file.'
                log_message(log_file, message)
                print(message)
                continue  # Skip to the next run
        else:
            if os.path.exists(out_path):
                log_message(log_file, f'Using existing raw VMR in preprocessing folder → {out_path}')
            elif os.path.exists(raw_vmr_path):
                doc_vmr = brainvoyager.open(raw_vmr_path)
                doc_vmr.save_as(out_path)
                copyfile((raw_vmr_path[:-4] + '.v16'), (out_path[:-4] + '.v16'))
                doc_vmr.close()
                log_message(log_file, f'Copied VMR from rawdata_bv → {out_path}')
            else:
                doc_vmr = brainvoyager.open(raw_nii_path)
                doc_vmr.save_as(out_path)
                doc_vmr.close()
                log_message(log_file, f'Converted NIfTI to VMR → {out_path}')


        # Step 1: Run enabled steps in pipeline_order
        for step in pipeline_order:
            if not options.get(step, False):
                continue  # Only run if enabled in options

            try:
                doc_vmr = brainvoyager.open(current_file)
                log_message(log_file, f'Starting step: {step} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

                if step == "transform_sag":
                    current_file = current_file[:-4]+ '_SAG.vmr'
                    doc_vmr.transform_to_std_sag(current_file)

                elif step == "transform_iso_voxel":
                    current_file = current_file[:-4] + f"_ISO-{options['target_res']}.vmr"
                    doc_vmr.transform_to_isovoxel(options["target_res"], options["framing_cube_dim"], options["interpolation_method"], current_file) 
                                     
                elif step == "correct_intensity_inhomogeneities":
                    doc_vmr.correct_intensity_inhomogeneities_ext(options["include_brain_extraction"], options["n_cycles"], options["tissue_range_thresh"], options["intensity_thresh"], options["fit_polynom_order"])
                    doc_vmr.save_as(doc_vmr.path_file_name)
                    current_file = doc_vmr.path_file_name

                doc_vmr.close()
                log_message(log_file, f'Finished step: {step} → {current_file}')
                
            except Exception as e:
                log_message(log_file, f'ERROR in step "{step}": {str(e)}')
                    
        log_message(log_file, f'=== Finished preprocessing for {sub_id} {ses_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n')

        print(f'=== Finished preprocessing for {sub_id} {ses_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n\n')
        print(f'Log file saved: {log_file}')

