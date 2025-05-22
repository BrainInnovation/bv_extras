"""
VTC Preprocessing Pipeline for BrainVoyager Data Analysis File and Folder Structure (BIDS)

This script performs preprocessing on VTC data in a BIDS-like directory structure
using BrainVoyager. It supports multiple subjects, sessions, and runs, with flexible
control over preprocessing steps and logging.

- Applies preprocessing steps such as:
    - Temporal high-pass filtering (using Fourier GLM, DCT GLM or FFT)
    - Temporal smoothing
    - Spatial smoothing
- Saves outputs in a dedicated vtc preprocessing workflow folder
- Logs every step of the preprocessing and captures any errors for later review

Inputs:
- VTC files saved in:
    derivatives/
        └── workflow_id-5_type-5_name-func-normalization/
            └── sub-XX/ses-YY/func/
                ├── sub-XX_ses-YY_task-XY_run-ZZ_bold_*.vtc
                
Outputs:
- VTC files saved in:
    derivatives/
        └── workflow_id-6_type-9_name-vtc-preprocessing/
            └── sub-XX/ses-YY/func/
                ├── sub-XX_ses-YY_task-XY_run-ZZ_bold_*_THPGLMF3c_SD3DVSS4.00mm.vtc
                └── log file (timestamped)


Configuration
-------------
Edit the configuration block to match your project structure and enable/disable/adjust 
the performed preprocessing steps in the options dictionary

"""

import os
from datetime import datetime


# ============================ CONFIGURATION ==============================

project_path = '/Users/judithe/Documents/BrainVoyager/Projects/NeWBI4fMRI2020/derivatives/'
vtc_input_suffix = '_3DMCTS_SCCTBL_256_sinc3_2x1.0_MNI.vtc' 
anat_folder = 'workflow_id-2_type-3_name-anat-normalization/'
anat_suffix = '_IIHC_MNI.vmr' 
vtc_input_folder = 'workflow_id-5_type-5_name-func-normalization/'
vtc_output_folder = 'workflow_id-6_type-9_name-vtc-preprocessing/'

subjects = [1,2,3,4,5,6,7,8,9,10,11,12]
sessions = [1]
runs = [1,2]
task_name = 'task-Localizer'

# Preprocessing options
options = {
# ---- temporal high‑pass filtering ---------------------------------------
'highpass_filter': True,
'hp_method': 'fourier',  # 'fourier' | 'dct' | 'fft'

#  if hp_method == 'fourier', GLM with sine/cosine predictors as well as a linear (and a constant) predictor
'fourier_n_cycles': 3, # number of cycles (number of sine/cosine pairs of predictors)

#  if hp_method == 'dct', GLM containing discrete cosine predictors (plus a constant predictor))
'dct_n_basis': 3, # number of discrete cosine predictors (basis functions)

#  if hp_method == 'fft', removes low-frequency drifts using Fourier analysis
'hp_unit': 'Hz', # 'Hz', 'cycles'
'hp_cutoff': 0.01, # e.g. 3 cycles or 0.01 Hz

# ---- temporal Gaussian smoothing ----------------------------------------
'smooth_temporal': False,
'temp_fwhm': 3,
'temp_fwhm_unit':'dps',  # 'd or dps (data points) or s or secs (seconds)

# ---- spatial Gaussian smoothing -----------------------------------------
'smooth_spatial': True,
'spatial_fwhm': 4,
'spatial_fwhm_unit': 'mm'  # 'mm' or 'vx'
}

# Full pipeline, change order to change the order of the preprocessing
pipeline_order = [
    "highpass_filter",
    "smooth_temporal",
    "smooth_spatial"
]

# ========================================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def pjoin(*args):
    return os.path.join(*args)

def log_message(log_file, msg):
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(f'{msg}\n')
        
        
# ========================= MAIN PROCESS ==================================
 
errors = []
bv.close_all()

for subj in subjects:
    sub_id = f'sub-{subj:02d}'
    for sess in sessions:
        ses_id = f'ses-{sess:02d}'
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')    
        
        out_path = pjoin(project_path, vtc_output_folder, sub_id, ses_id, 'func')
        ensure_dir(out_path)
        log_file = pjoin(out_path, f'{sub_id}_{ses_id}_{task_name}_vtc_preprocessing_{timestamp}.log')

        log_message(log_file, f'VTC Preprocessing Log: {sub_id} {ses_id} {task_name}')
        log_message(log_file, f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log_message(log_file, "Preprocessing options:")

        for key, value in options.items():
            log_message(log_file, f"  {key}: {value}")

        log_message(log_file, "")  # Blank line for readability
        
        print(f'=== Starting preprocessing for {sub_id} {ses_id} {task_name} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n')
        
        try:
            vmr_file = f'{sub_id}_{ses_id}_T1w{anat_suffix}'
            vmr_path = pjoin(project_path, anat_folder, sub_id, ses_id, 'anat', vmr_file)
            doc_vmr = brainvoyager.open(vmr_path)

        except Exception as e:
            msg = f'ERROR: {sub_id} {ses_id}: Failed to locate VMR: {e}'
            print(msg)
            errors.append(msg)
            log_message(log_file, msg)
            continue

        for run in runs:
            run_id = f'run-{run:02d}'
            vtc_input = pjoin(project_path, vtc_input_folder, sub_id, ses_id, 'func', f'{sub_id}_{ses_id}_{task_name}_{run_id}_bold{vtc_input_suffix}')

            if not os.path.exists(vtc_input):
                msg = (
                    f'ERROR: Could not find VTC: {vtc_input} '
                    f'Skipping this run.'
                )
                print(msg)
                errors.append(msg)
                log_message(log_file, msg)
                continue  # Skip to the next run
             
            # copy vtc input file to vtc preprocessing folder
            success = doc_vmr.link_vtc(vtc_input)
            if success == False:
                msg = (
                    f'ERROR: Could not link VTC: {vtc_input} '
                    f'Skipping this run.'
                )
                print(msg)
                errors.append(msg)
                log_message(log_file, msg)
                continue  # Skip to the next run
                
            current_vtc = pjoin(project_path, vtc_output_folder, sub_id, ses_id, 'func', f'{sub_id}_{ses_id}_{task_name}_{run_id}_bold{vtc_input_suffix}')
            doc_vmr.save_vtc(current_vtc)
            
            for step in pipeline_order:
                if not options.get(step, False):
                    continue  # Only run if enabled in options

                try:
                    doc_vmr.link_vtc(current_vtc)
                    log_message(log_file, f'{task_name}_{run_id}: Starting step: {step} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
                    
                    if step == "highpass_filter":
                        if options['hp_method'] == 'fourier':
                            doc_vmr.filter_temporal_highpass_glm_fourier(options["fourier_n_cycles"])
                        elif options['hp_method'] == 'dct':
                            doc_vmr.filter_temporal_highpass_glm_dct(options["dct_n_basis"])
                        elif options['hp_method'] == 'fft':
                            doc_vmr.filter_temporal_highpass_fft(options["hp_cutoff"], options["hp_unit"])
                        else:
                            raise ValueError(f'Unknown Temporal Highpass Filtering Method: {options["hp_method"]}')
                            
                    elif step == "smooth_temporal":
                        doc_vmr.smooth_temporal(options["temp_fwhm"], options["temp_fwhm_unit"])
                    
                    elif step == "smooth_spatial":
                        doc_vmr.smooth_spatial(options["spatial_fwhm"], options["spatial_fwhm_unit"])
                        
                    current_vtc = doc_vmr.vtc_file
                    log_message(log_file, f'{task_name}_{run_id}: Finished step: {step} → {current_vtc}')
                    
                except Exception as e:
                    msg = f'ERROR: {sub_id}_{ses_id}_{task_name}_{run_id}: ERROR in step "{step}": {str(e)}'
                    errors.append(msg)
                    log_message(log_file, msg)
                                            
            log_message(log_file, f'=== Finished preprocessing for {sub_id}_{ses_id}_{task_name}_{run_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n')
        
            # delete vtc input file from vtc preprocessing folder
            try:
                os.remove(pjoin(project_path, vtc_output_folder, sub_id, ses_id, 'func', f'{sub_id}_{ses_id}_{task_name}_{run_id}_bold{vtc_input_suffix}')) 
            except FileNotFoundError:
                print("File already removed or never existed.")
            except PermissionError:
                print("You don’t have permission to delete this file.")   
                
        doc_vmr.close()
  
        log_message(log_file, f'=== Finished preprocessing for {sub_id} {ses_id} ===\n')
        print(f'=== Finished preprocessing for {sub_id} {ses_id} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ===\n\n')
        print(f'Log file saved: {log_file}')
        

# === FINAL ERROR SUMMARY ===
if errors:
    print('\n=== SUMMARY OF ERRORS ===')
    for err in errors:
        print(f'- {err}')
else:
    print('\nAll VTC preprocessing completed without errors.')


        
