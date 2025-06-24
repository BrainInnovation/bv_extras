"""
Create a BrainVoyager SDM file from a PRT file using a canonical two-gamma HRF.

This script generates a single-study design matrix (SDM) based on conditions defined 
in a BrainVoyager protocol (PRT) file. A canonical two-gamma hemodynamic response 
function (HRF) is used to model the expected BOLD signal.

Inputs:
- A PRT file is required.
- An FMR or VTC file is optional and used to automatically extract the number of volumes (n_vols)
  and repetition time (TR). 
- If no functional file is provided, the user will be prompted to enter the TR (in milliseconds)
  and the number of volumes manually.

Features:
- Supports parametric modulation if the PRT contains parametric weights and 
  parametric SDM generation is enabled.
- Allows optional standardization of parametric weights (raw, demeaned, or z-scored).
- Optionally removes a specified "Rest" condition from the PRT.
- Scales predictors to unit amplitude if specified.

Output:
- A BrainVoyager-compatible SDM file saved alongside the input PRT.

"""

import numpy as np
from scipy.stats import gamma
from scipy.signal import convolve, resample_poly
import bvbabel

# =========================== CONFIGURATION ===================================

# to be defined only for parametric PRTs: choose whether you would like to add parametric predictors to the SDM
define_para_sdm = True
# to be defined only for parametric SDMs: use parametric weights as defined in PRT (0) or subtract mean of weights (1) or z-score weights (2)
standardize_weights = 0

# define if the PRT includes a Rest condition
RestCondInPRT = False
RestCondNumber = 0 # 0 - for first condition or -1 - for last condition

# Scale the amplitudes of the final SDM predictors to 1 (not valid for parametric predictors)
ScalePred1 = False

# read all necessary files (PRT, optionally FMR or VTC)
prt_file = '/Users/judithe/Documents/BrainVoyager/SampleData/BV_GSG/derivatives/rawdata_bv/sub-01/ses-04/func/sub-01_ses-04_task-blocked_run-01_bold.prt'
prt_header, prt_data = bvbabel.prt.read_prt(prt_file)
func_file = None # or path to FMR or VTC file


# =============================================================================

if func_file:
    if func_file.endswith(".fmr"):
        func_header, _  = bvbabel.fmr.read_fmr(func_file)
        n_vols = func_header['NrOfVolumes']
        TR = int(func_header['TR'])
    elif func_file.endswith(".vtc"):
        func_header, _  = bvbabel.vtc.read_vtc(func_file)
        n_vols = func_header['Nr time points']
        TR = int(func_header['TR (ms)'])
    else:
        raise ValueError("Unsupported file type. Only .fmr or .vtc are supported.")
else:
    # --- Prompt user for TR and number of volumes ---
    try:
        TR = int(input("Enter TR in milliseconds (e.g., 2000): "))
        n_vols = int(input("Enter number of volumes (e.g., 160): "))
    except ValueError:
        raise ValueError("Invalid input: TR and number of volumes must be integers.")
        

# define how parametric weights should be treated for SDM creation
def create_weights(a, b):
    '''
    create weight vector (unchanged, demeaned, zscored)
    Parameters:
        a: Array containing the weights of the condition
        b: 0 - unchanged, 1 - demeaned, 2 - zscored
    '''
    if b == 0:
        weights = a
    elif b == 1:
        weights = a - np.mean(a)
    elif b == 2:
        weights = (a - np.mean(a)) / np.std(a)
    return(weights)


# --- HRF computation ---
def twogamma_hrf(params=None, fs=100, ScaleToOne=False,DivideBySum=True):
    """
    Compute a canonical two-gamma HRF.

    Args:
        params (list or None): HRF parameters. Defaults to SPM [6, 16, 1, 1, 6, 0, 32].
        fs (int): Sampling frequency in Hz. Default is 100.
        ScaleToOne (bool): If True, scale HRF to unit peak.
        DivideBySum (bool): If True, normalize HRF by its sum.

    Returns:
        tuple: (hrf: np.ndarray, fs: int)
    """
    
    if params is None:
        params = [6, 16, 1, 1, 6, 0, 32]
    
    dt = 1 / fs # 1/Hz
    n_subs = int(params[6] / dt)
    u = np.arange(n_subs) - params[5] / dt

    posG = gamma.pdf(u * dt, params[0] / params[2], scale=params[2])
    negG = gamma.pdf(u * dt, params[1] / params[3], scale=params[3])

    hrf = posG - negG / params[4] if params[4] >= 0.001 else -negG

    max_hrf = np.max(hrf)
    min_hrf = np.min(hrf)
    sum_hrf = np.sum(hrf)
    abs_fmax = max(abs(min_hrf), abs(max_hrf))

    if ScaleToOne or sum_hrf <= 0.01:
        DivideBySum = False
    elif DivideBySum:
        ScaleToOne = False

    if abs_fmax < 0.0001:
        ScaleToOne = False

    if ScaleToOne:
        hrf /= abs_fmax
    if DivideBySum:
        hrf /= sum_hrf

    return hrf, fs




# ========================= MAIN PROCESS ======================================

# compute the HRF
hrf, fs = twogamma_hrf()

# determine PRT time resolution
prt_is_msec = prt_header['ResolutionOfTime'].lower() == 'msec'

### Delete rest condition from PRT, if included
if RestCondInPRT:
    del prt_data[RestCondNumber]
    prt_header['NrOfConditions'] = int(prt_header['NrOfConditions'])-1


sdm_data = [] # create empty design
  
for cond in range(len(prt_data)):
    ###  Check if parametric weights are defined in PRT and should be included
    add_parametric_predictor = (
        "ParametricWeights" in prt_header
        and prt_header["ParametricWeights"] == 1
        and define_para_sdm
        # if the parametric weights in this condition are not constant create a parametric predictor
        and len(np.unique(prt_data[cond]['Parametric weight'])) > 1
    )
    # Create weights if a parametric predictor is needed
    if add_parametric_predictor:
        weights = create_weights(prt_data[cond]['Parametric weight'], standardize_weights)
        
    ### Define the convolved predictors based on the PRT conditions
    # boxcar is created in second resolution and sampled to fs
    n_timepoints_sec = (n_vols * TR) / 1000
    n_timepoints_hires = int(n_timepoints_sec * fs)
    pred_block_hires = np.zeros(n_timepoints_hires, dtype=float)
    
    for start, stop in zip(prt_data[cond]['Time start'], prt_data[cond]['Time stop']):
        if prt_is_msec:
            start_sec = start / 1000
            stop_sec = stop  / 1000
        else:
            start_sec = start * TR / 1000 - TR/1000
            stop_sec = stop * TR / 1000
        start_idx = int(start_sec * fs)
        stop_idx = int(stop_sec * fs)
        pred_block_hires[start_idx:stop_idx+1] = 1.0
        
        
    pred_conv_hires = convolve(pred_block_hires, hrf, mode='full')[:len(pred_block_hires)]
    samples_per_TR = int(fs * TR / 1000)
    pred_conv = resample_poly(pred_conv_hires, up=1, down=samples_per_TR)[:n_vols]   
  
        
    if add_parametric_predictor:
        pred_weights_hires = np.zeros(n_timepoints_hires, dtype=float)
        for start, stop, weight in zip(
            prt_data[cond]['Time start'],
            prt_data[cond]['Time stop'],
            weights
        ):
            if prt_is_msec:
                start_sec = start / 1000
                stop_sec = stop  / 1000
            else:
                start_sec = start * TR / 1000 - TR/1000
                stop_sec = stop * TR / 1000
            start_idx = int(start_sec * fs)
            stop_idx = int(stop_sec * fs)
            pred_weights_hires[start_idx:stop_idx+1] = weight
         
        #create convolved weighted predictor
        pred_weights_conv_hires = convolve(pred_weights_hires, hrf, mode='full')[:len(pred_weights_hires)]
        pred_weights_conv = resample_poly(pred_weights_conv_hires, up=1, down=samples_per_TR)[:n_vols]
        
### Define Single Study Design Matrix based on the created predictors
    
    ## Add main predictor to model
    # scale predictor to 1  
    if ScalePred1:  
        pred_conv = pred_conv/max(pred_conv)
    # copy to design
    sdm_data.append({})
    sdm_data[-1]['NameOfPredictor'] = prt_data[cond]['NameOfCondition']
    sdm_data[-1]['ColorOfPredictor'] = prt_data[cond]['Color'].tolist()
    sdm_data[-1]['ValuesOfPredictor'] = pred_conv
    del(pred_block_hires, pred_conv_hires, pred_conv)
    
    ## Add parametric predictor to model
    if add_parametric_predictor:
        # copy to design
        sdm_data.append({})
        sdm_data[-2]['NameOfPredictor'] = prt_data[cond]['NameOfCondition'] + ' [Main]' # change name of main predictor
        sdm_data[-1]['NameOfPredictor'] = prt_data[cond]['NameOfCondition'] + ' [Parametric]'
        sdm_data[-1]['ColorOfPredictor'] = prt_data[cond]['Color'].tolist()
        sdm_data[-1]['ValuesOfPredictor'] = pred_weights_conv
        del(pred_weights_hires, pred_weights_conv_hires, pred_weights_conv, weights)

# Add constant to model
sdm_data.append({})
sdm_data[-1]['NameOfPredictor'] = 'Constant'
sdm_data[-1]['ColorOfPredictor'] = [255,255,255]
sdm_data[-1]['ValuesOfPredictor'] = np.ones(n_vols).astype(float)

# Define SDM header
sdm_header = {'FileVersion': 1, 'NrOfPredictors': len(sdm_data), 'NrOfDataPoints': int(n_vols),
              'IncludesConstant': 1, 'FirstConfoundPredictor': len(sdm_data)}

# Save SDM
bvbabel.sdm.write_sdm((prt_file[:-4] + '.sdm'), sdm_header, sdm_data)

