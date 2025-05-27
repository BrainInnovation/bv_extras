
"""
Motion Qulaity Control and Visualization for BrainVoyager BIDS Data

The script is structured for multi-subject, multi-session, and multi-run pipelines
using a BIDS-like directory layout. It processes motion SDM (*3DMC.sdm) files to generate:
- Temporal derivatives, squared, and combined motion regressors
- Framewise displacement (FD) and spike motion predictors, based on a user-defined FD-threshold
- Quality control plots of motion, FD, optional PRT overlays, and carpet plots

Inputs per functional run:
    - *3DMC.sdm files (as the output from motion correction in BV)
    - optional: a corresponding preprocessed FMR file to visualize the effect of resdiual motion on the (z-scored) voxel time courses in a carpet plot
    - optional: a stimuluation protocol of the FMR (PRT) to visualize the subject motion in relation to the experimental stimulation
    
Outputs per functional run (saved in the same folder as the input *3DMC.sdm file):
    - *_zscored.sdm (original zscored)
    - *_derivative.sdm (zscored)
    - *_squared.sdm (zscored)
    - *_derivative_squared.sdm (zscored)
    - *_*model.sdm, combined zscored motion model (12 - original + temporal derivatives, 18 - original + temporal derivatives + squared, 24 - original + temporal derivatives + squared + squared tempral derivatives)
    - *_FDspikes.sdm
    - *_FDspikemodel.sdm
    - *_motion_plot.png
    - per-session-task logfile
    
Configuration:
- You must adapt the configuration section at the top of this script to match:
    - Your project folder structure
    - Workflow folder names
    - Filename suffixes 
    - List of subject IDs, session IDs, and runs
    - Spike framewise displacement threshold (fd_threshold) to identify functional volumes with spike motion (e.g. lenient threshold = 0.5, stringent threshold = 0.2)
    - mean voxel intensity threshold (intensity_thresh) to exclude background voxels from the carpet plot, if a carpet plot should be generated
    - the complexity of the generated motion model, i.e. the number of included motion regressors (model_type = 12,18,24)
    - !the prt file name has to be checked and possibly adapted by the user in the first occurence of prt_name in the script, as prt naming conventions vary!

Requires:
- bvbabel
- BIDS-style folder structure under `project_path`

some background literature:
    Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3), 2142?2154. https://doi.org/10.1016/J.NEUROIMAGE.2011.10.018
    Power, J. D., Lynch, C. J., Silver, B. M., Dubin, M. J., Martin, A., & Jones, R. M. (2019). Distinctions among real and apparent respiratory motions in human fMRI data. NeuroImage, 201, 116041. https://doi.org/10.1016/j.neuroimage.2019.116041
    Power, J. D., Mitra, A., Laumann, T. O., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2014). Methods to detect, characterize, and remove motion artifact in resting state fMRI. NeuroImage, 84, 320?341. https://doi.org/10.1016/J.NEUROIMAGE.2013.08.048
"""

import os
import math
import copy
import numpy as np
from datetime import datetime
from scipy.stats import zscore
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import scipy.signal as signal
import bvbabel


# === CONFIGURATION ===
project_path = '/Users/judithe/Documents/BrainVoyager/Projects/NeWBI4fMRI2020/derivatives'
motion_folder = 'workflow_id-3_type-1_name-func-preprocessing'
motion_suffix = '_3DMC.sdm' # (or _SCCTBL_3DMC.sdm depending on the FMR preprocessing order)
model_type = '12' # Options: 12', '18', '24'
fd_threshold = 0.2

subjects = [1,2,3,4,5,6,7,8,9,10,11,12]
sessions = [1]
runs = [1,2]
task_name = 'task-Localizer'

PLOTTING = True
USE_PRT = True
prt_folder = 'rawdata_bv/'
TR = 1000

USE_CARPET = True
fmr_folder = 'workflow_id-3_type-1_name-func-preprocessing'
fmr_suffix = '_3DMCTS.fmr' # (or the fully preprocessed fmr, including motion correction e.g _3DMCTS_SCCTBL_THPGLMF3c.fmr)
intensity_thresh = 1000

# === INTERNAL SETTINGS ===

# Transformation definitions for variants
variant_specs = [
    ("zscored",            lambda m, d, i: zscore(m[i]),               ' zscored'),
    ("derivative",         lambda m, d, i: zscore(d[i]),               ' derivative'),
    ("squared",            lambda m, d, i: zscore(np.square(m[i])),    ' squared'),
    ("derivative_squared", lambda m, d, i: zscore(np.square(d[i])),    ' derivative_squared'),
]

def deg2mm(deg, head_radius=50):
    '''
    convert rotational displacements from degrees to millimeters 
    by calculating displacement on the surface of a sphere of
    radius 50 mm, assuming a standard head size of 50 mm (see Power er al., 2012)
    '''
    return (deg * math.pi / 180) * head_radius

def compute_framewise_displacement(motion, head_radius=50):
    """
    Computes framewise displacement (FD) for 6 motion parameters.
    First 3: translations (mm), last 3: rotations (degrees → mm).
    """
    # subtract first time point of motion to reflect within-run motion
    motion_zero = motion - motion[:, [0]]

    # Convert rotations from degrees to mm
    motion_zero[3:] = deg2mm(motion_zero[3:], head_radius=head_radius)

    # Compute framewise displacement: sum of absolute differences across all 6 params
    fd = np.sum(np.abs(np.diff(motion_zero, axis=1)), axis=0)

    # Prepend zero to match time series length
    fd = np.insert(fd, 0, 0)

    return fd

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def log_message(log_file, msg):
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(f'{msg}\n')

def get_output_paths(motionfile):
    base = motionfile[:-4]
    paths = {
        'fdspikes': f'{base}_FDspikes.sdm',
        'fd': f'{base}_FD.sdm',
        'plot': f'{base}_motion_plot.png',
    }
    for name, _, _ in variant_specs:
        paths[name] = f'{base}_{name}.sdm'
    for k in ['12model', '18model', '24model']:
        paths[k] = f'{base}_{k}.sdm'
    return paths

def create_motion_plot(motion, fd, voxel_by_time_z, sdm_data, prt_data, prt_header, TR, rms, max_motion, n_spikes,
                       fd_threshold, rms_color, fmr_header, motion_file, base_name, motion_suffix, PRT=True, CARPET=True):
    
    fig = plt.figure(figsize=(10, 8))

    ax_motion = fig.add_subplot(311 if CARPET else 211)
    ax_motion.set_title('Motion: ' + base_name + motion_suffix)
    ax_motion.set_ylim(np.min(motion) - 0.1, np.max(motion) + 0.1)
    ax_motion.set_ylabel('BV Motion Parameters')
    ax_motion.grid(axis='y')

    for i in range(6):
        ax_motion.plot(motion[i], linewidth=1,
                       color=np.array(sdm_data[i]['ColorOfPredictor']) / 255,
                       label=sdm_data[i]['NameOfPredictor'])

    handles, _ = ax_motion.get_legend_handles_labels()

    if PRT:
        for cond in prt_data:
            for start, stop in zip(cond['Time start'], cond['Time stop']):
                start = start / TR if prt_header['ResolutionOfTime'] == 'msec' else start
                stop = stop / TR if prt_header['ResolutionOfTime'] == 'msec' else stop + 1
                duration = stop - start
                ax_motion.add_patch(Rectangle((start, np.min(motion) - 0.1), duration,
                                              np.max(motion) - np.min(motion) + 0.2,
                                              facecolor=np.array(cond['Color']) / 255,
                                              edgecolor='none', alpha=0.3))
        handles += [Patch(color=np.array(cond['Color']) / 255, label=cond['NameOfCondition']) for cond in prt_data]

    ax_motion.axes.xaxis.set_ticklabels([])
    x_limits = ax_motion.get_xlim()

    ax_fd = fig.add_subplot(312 if CARPET else 212)
    ax_fd.set_ylim(0, np.max(fd) + 0.1)
    ax_fd.set_ylabel('FD (mm)')
    ax_fd.grid(axis='y')
    ax_fd.plot(fd, linewidth=1, color='black', label='Framewise Displacement (FD)')
    ax_fd.plot(np.ones_like(fd) * fd_threshold, linestyle='dotted', linewidth=1, color='red',
               label=f'FD Threshold = {fd_threshold}')
    if not CARPET:
        ax_fd.set_xlabel('Volume #')
    else:
        ax_fd.axes.xaxis.set_ticklabels([])

    handles_fd, _ = ax_fd.get_legend_handles_labels()
    quality_info = [
        Patch(color='white', label='Motion Quality Measures:'),
        Patch(color='white', label=f'Mean FD = {np.mean(fd):.3f}'),
        Patch(color='white', label=f'# Spikes based on FD({fd_threshold}) = {n_spikes}'),
        Patch(color='white', label=f'Max Motion within Run = {max_motion:.3f} mm'),
        Patch(color='white', label=f'RMS Motion in Run = {rms:.3f}')
    ]
    handles += handles_fd + quality_info

    if CARPET and voxel_by_time_z is not None:
        ax_carpet = fig.add_subplot(313)
        ax_carpet.imshow(voxel_by_time_z, aspect='auto', cmap='gray')
        ax_carpet.set_xlabel('Volume #')
        ax_carpet.set_xlim(x_limits)
        ax_carpet.set_ylabel('Voxel #')
        ax_carpet.set_title(f'zscored {fmr_header["Prefix"]}.fmr')

    legend = fig.legend(handles=handles, loc="center right", frameon=False, framealpha=1, fontsize='small')
    for text in legend.get_texts():
        if 'RMS Motion' in text.get_text():
            text.set_color(rms_color)

    fig.tight_layout(rect=[0, 0, 0.75, 1])
    fig.savefig(f'{motion_file[:-4]}_motion_plot.png', dpi=300)
    plt.close(fig)


# === MAIN PROCESSING ===

errors = []
assert model_type in ('12', '18', '24'), f'Invalid MODEL_TYPE: {model_type}. Use "12", "18", or "24".'

for subj in subjects:
    sub_id = f'sub-{subj:02d}'
    for sess in sessions:
        ses_id = f'ses-{sess:02d}'
        
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = os.path.join(project_path, motion_folder, sub_id, ses_id, f'{sub_id}_{ses_id}_{task_name}_motionqc_{timestamp}.log')
        log_message(log_file, f'Motion Quality Control Log: {sub_id} {ses_id} {task_name}')
        log_message(log_file, f'Number of Motion Regressors in Motion Model: {model_type}')
        log_message(log_file, f'Framewise Displacement Spike Threshold: {fd_threshold} mm')
        log_message(log_file, f'Generate Motion Plots: {PLOTTING}')
        log_message(log_file, f'Use PRT for Motion Plots: {USE_PRT}')
        log_message(log_file, f'Create Carpet Plots for Motion Plots: {USE_CARPET}')
        
        
        for run in runs:
            run_id = f'run-{run:02d}'
            base_name = f'{sub_id}_{ses_id}_{task_name}_{run_id}_bold'
            
            motion_path = os.path.join(project_path, motion_folder, sub_id, ses_id, 'func')
            motion_file = os.path.join(motion_path, base_name + motion_suffix)
            
            log_message(log_file, f'\nProcessing {base_name} at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
            print(f'Processing {base_name}')
            
            
            if not os.path.exists(motion_file):
                msg = (
                    f'ERROR: Could not find Motion SDM: {motion_file} '
                    f'Skipping this run.'
                )
                print(msg)
                errors.append(msg)
                log_message(log_file, msg)
                continue  # Skip to the next run
            
            sdm_header, sdm_data = bvbabel.sdm.read_sdm(motion_file)
            out_files = get_output_paths(motion_file)
            
            carpet_available = USE_CARPET            
            if USE_CARPET:
                fmr_file = os.path.join(project_path, fmr_folder, sub_id, ses_id, 'func', base_name + fmr_suffix)
                if not os.path.exists(fmr_file):
                    msg = (
                        f'ERROR: Could not find FMR: {fmr_file} '
                        f'No carpet plot will be created for this run.'
                    )
                    print(msg)
                    errors.append(msg)
                    log_message(log_file, msg)
                    carpet_available= False
                else:
                    fmr_header, fmr_data = bvbabel.fmr.read_fmr(fmr_file)
                    # use a mean intensity threshold for excluding background noise
                    voxel_by_time = fmr_data.reshape(np.prod([fmr_header['ResolutionX'], fmr_header['ResolutionY'], fmr_header['NrOfSlices']]),fmr_header['NrOfVolumes'] )
                    voxel_by_time_thresh = voxel_by_time[np.mean(voxel_by_time, axis=1)>intensity_thresh]
                    voxel_by_time_z = np.nan_to_num(zscore(voxel_by_time_thresh, axis=-1))
                
            prt_available = USE_PRT
            if USE_PRT:
                prt_path = os.path.join(project_path, prt_folder, sub_id, ses_id, 'func')
                prt_name = f'{sub_id}_{ses_id}_{task_name[5:]}_run-{run:01d}.prt'
                prt_file = os.path.join(prt_path, prt_name)
                if not os.path.exists(prt_file):
                    msg = (
                        f'ERROR: Could not find PRT: {prt_file} '
                        f'The motion will be plotted without the experimental protocol.'
                    )
                    print(msg)
                    errors.append(msg)
                    log_message(log_file, msg)
                    prt_available = False
                else:
                    prt_header, prt_data = bvbabel.prt.read_prt(prt_file)
          
            
            # === WRITE SDM VARIANTS ===
          
            motion = np.vstack([sdm_data[i]['ValuesOfPredictor'] for i in range(6)])
            motion_d = np.insert(np.diff(motion, axis=1), 0, 0, axis=1)
          
            for name, func, postfix in variant_specs:
                    sdm_copy = copy.deepcopy(sdm_data)
                    for i in range(6):
                        sdm_copy[i]['ValuesOfPredictor'] = func(motion, motion_d, i)
                        sdm_copy[i]['NameOfPredictor'] += postfix
                    bvbabel.sdm.write_sdm(out_files[name], sdm_header, sdm_copy)
                    log_message(log_file, f'Wrote {out_files[name]}')
                    
            # === Build combined z-scored motion model (12, 18, 24 only) ===
            combined_model = []
            
            # Always include zscored original motion (6 predictors)
            for i in range(6):
                pred = copy.deepcopy(sdm_data[i])
                pred['ValuesOfPredictor'] = zscore(motion[i])
                pred['NameOfPredictor'] += ' zscored'
                combined_model.append(pred)
            
            # Add derivatives if model is 12 or higher
            if model_type in ('12', '18', '24'):
                for i in range(6):
                    pred = copy.deepcopy(sdm_data[i])
                    pred['ValuesOfPredictor'] = zscore(motion_d[i])
                    pred['NameOfPredictor'] += ' derivative'
                    combined_model.append(pred)
            
            # Add squared motion if model is 18 or higher
            if model_type in ('18', '24'):
                for i in range(6):
                    pred = copy.deepcopy(sdm_data[i])
                    pred['ValuesOfPredictor'] = zscore(np.square(motion[i]))
                    pred['NameOfPredictor'] += ' squared'
                    combined_model.append(pred)
            
            # Add squared derivatives if model is 24
            if model_type == '24':
                for i in range(6):
                    pred = copy.deepcopy(sdm_data[i])
                    pred['ValuesOfPredictor'] = zscore(np.square(motion_d[i]))
                    pred['NameOfPredictor'] += ' derivative_squared'
                    combined_model.append(pred)

            header = copy.deepcopy(sdm_header)
            header['NrOfPredictors'] = len(combined_model)
            output = out_files[f'{model_type}model']
            bvbabel.sdm.write_sdm(output, header, combined_model)
            
            log_message(log_file, f'Wrote {output}')
            

            # === Calculate the Framewise Displacement (FD) ===
 
            fd = compute_framewise_displacement(motion)
            # Save FD as a single regressor SDM
            fd_sdm = copy.deepcopy(sdm_data[:1])  # Use first predictor as template
            fd_sdm[0]['ValuesOfPredictor'] = fd
            fd_sdm[0]['NameOfPredictor'] = 'Framewise Displacement'
            fd_sdm[0]['ColorOfPredictor'] = [0, 0, 0]  # black

            fd_header = copy.deepcopy(sdm_header)
            fd_header['NrOfPredictors'] = 1
            fd_header['IncludesConstant'] = 0
            fd_header['FirstConfoundPredictor'] = 1

            bvbabel.sdm.write_sdm(out_files['fd'], fd_header, fd_sdm)
            log_message(log_file, f'Wrote {out_files["fd"]}')
            log_message(log_file,f'Mean Framewise Displacement (FD) = {np.mean(fd):.3f} mm')
            
            
            # === Calculate the Motion Spikes ===
            
            spikes = (fd > fd_threshold).astype(int)
            if spikes.any():
                spikes_sdm = []
                for i, vol in enumerate(np.where(spikes)[0]): 
                    p = np.zeros_like(spikes)
                    p[vol] = 1
                    spikes_sdm.append({
                        'ValuesOfPredictor': p,
                        'NameOfPredictor': f'Spike_{vol+1}',
                        'ColorOfPredictor': [255 - i*10 % 255, 0, 0]
                    })
                spikes_header = copy.deepcopy(sdm_header)
                spikes_header['NrOfPredictors'] = len(spikes_sdm)
                spikes_header['IncludesConstant'] = 0
                spikes_header['FirstConfoundPredictor'] = 1
                bvbabel.sdm.write_sdm(out_files['fdspikes'], spikes_header, spikes_sdm)
                n_spikes = int(np.sum(spikes))
            else:
                n_spikes = 0
            log_message(log_file, f'Number of Spikes based on FD spike threshold ({fd_threshold} mm) = {n_spikes} ')
                
                
            # === Calculate the Root Mean Square Displacement ===
            # Subtract the first timepoint from each predictor to zero-center and convert to mm
            motion_centered = motion - motion[:, [0]]
            motion_centered_mm= np.concatenate((motion_centered[0:3,:],deg2mm(motion_centered[3:])), axis =0)
            max_motion = np.max(np.abs(motion_centered_mm))
            motion_detrended = signal.detrend(motion_centered_mm, axis=1)
            rms = np.mean(np.sqrt(np.sum(np.square(motion_detrended), axis=0)))
            
            log_message(log_file, f'Maximum Motion within Run = {max_motion:.3f} mm')
            log_message(log_file, f'Root Mean Square Displacement (RMS) in Run = {rms:.3f}')
            
            # low motion subjects below 0.2 mm, moderate motion between 0.2 and 0.5 mm, high motion above 0.5 mm (used for color-coding the RMS legend as green, yellow, red)
            # see Ciric et al. (2017) DOI: 10.1016/j.neuroimage.2017.03.020 and Satterthwaite et al., (2013) DOI: 10.1016/j.neuroimage.2013.05.064
            if rms < 0.2:
                rms_color = [0.0, 0.6, 0.0]  # green
            elif rms < 0.5:
                rms_color = [0.85, 0.65, 0.0]  # yellow-orange
            else:
                rms_color = [0.8, 0.0, 0.0]  # red
 
            # === Generate the Motion Plots ===
            
            if PLOTTING:
                create_motion_plot(
                    motion=motion,
                    fd=fd,
                    voxel_by_time_z=voxel_by_time_z if carpet_available else None,
                    sdm_data=sdm_data,
                    prt_data=prt_data if prt_available else [],
                    prt_header=prt_header if prt_available else {},
                    TR=TR if prt_available else [],
                    rms=rms,
                    max_motion=max_motion,
                    n_spikes=n_spikes,
                    fd_threshold=fd_threshold,
                    rms_color=rms_color,
                    fmr_header=fmr_header if carpet_available else {},
                    motion_file=motion_file,
                    base_name=base_name,
                    motion_suffix=motion_suffix,
                    PRT=prt_available,
                    CARPET=carpet_available
                )
                log_message(log_file, f'Wrote {out_files["plot"]}')


# === FINAL ERROR SUMMARY ===
if errors:
    print('\n=== SUMMARY OF ERRORS ===')
    for err in errors:
        print(f'- {err}')
else:
    print('\nAll SDM motion files processed without errors. Please check the logfiles and the plots for detailed information.')