"""
This script replicates the SDM creation in BrainVoyager, using the canonic two-gamma HRF.
If the PRT contains parametric weights, a parametric SDM can be created.

"""

import numpy as np
from scipy.stats import gamma
import scipy.interpolate as interpolate
from scipy import signal
import bvbabel


# choose whether you would like to add parametric predictors to the SDM, if the PRT contains weighted condition intervals
define_para_sdm = False
# to be defined for parametric SDMs: use parametric weights as defined in PRT (0) or subtract mean of weights (1) or z-score weights (2)
standardize_weights = 1

# read all necessary files (PRT, FMR)
fmr_file = '/Users/judithe/Documents/BrainVoyager/SampleData/BV_GSG/derivatives/rawdata_bv/sub-01/ses-04/func/sub-01_ses-04_task-blocked_run-01_bold.fmr'
fmr_header, _  = bvbabel.fmr.read_fmr(fmr_file)
prt_file = '/Users/judithe/Documents/BrainVoyager/SampleData/BV_GSG/derivatives/rawdata_bv/sub-01/ses-04/func/sub-01_ses-04_task-blocked_run-01_bold.prt'
prt_header, prt_data = bvbabel.prt.read_prt(prt_file)

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


sdm_data = list() # create empty design


###  Define HRF as specified in BrainVoyager's Two-Gamma HRF Notebook

a1 =  6.0; sc1 =  6.0 
a2 = 16.0; sc2 = -1.0
if prt_header['ResolutionOfTime'] != 'Volumes':
    hz = 1000
else:
    hz = 1/int(fmr_header['TR'])*1000
    
x = np.linspace(0, 30, int(hz*30+1)) # input range: 0 - 30 seconds
y1 = sc1 * gamma.pdf(x, a1)
y2 = sc2 * gamma.pdf(x, a2)
hrf = y1 + y2
hrf = hrf/max(hrf)



### Check whether there are parametric weights in the current PRT condition         
   
for cond in range(len(prt_data)):
    
    add_parametric_predictor = False
    if "ParametricWeights" in prt_header:
        if prt_header['ParametricWeights'] == 1:
            if define_para_sdm:
                # if the parametric weights in this condition are not constant create a parametric predictor
                if len(np.unique(prt_data[cond]['Parametric weight'])) > 1:
                    add_parametric_predictor = True
                    weights = create_weights(prt_data[cond]['Parametric weight'], standardize_weights)
                
### Define the convolved predictors based on the PRT conditions

    # if PRT is defined in msec, the convolution is performed in msec resolution and the resulting predictor is downsampled to volume resolution
    if prt_header['ResolutionOfTime'] != 'Volumes':
    
        # create empty predictor in MSEC resolution
        pred_empty = np.zeros(fmr_header['NrOfVolumes']*int(fmr_header['TR'])).astype(float)                        
        # create block predictor
        pred_block = pred_empty
        for event in range(len(prt_data[cond]['Time start'])):
            pred_block[prt_data[cond]['Time start'][event]:prt_data[cond]['Time stop'][event]+1] = 1
    
        # create convolved predictor
        pred_block = np.concatenate((np.zeros(len(hrf)-1), pred_block)) # pad predictor with zeros
        pred_conv = signal.convolve(pred_block, hrf, 'valid', method='auto')
                
        # downsample predictor to TR
        # create time vectors
        time_msec = np.linspace(0, fmr_header['NrOfVolumes']*int(fmr_header['TR']), fmr_header['NrOfVolumes']*int(fmr_header['TR']))
        time_TR= np.arange(time_msec[0], time_msec[-1], int(fmr_header['TR']))
        f1 = interpolate.interp1d(time_msec, pred_conv)
        pred_conv = f1(time_TR)
        
        if add_parametric_predictor:
            pred_weights = np.zeros(fmr_header['NrOfVolumes']*int(fmr_header['TR'])).astype(float)
            for event in range(len(prt_data[cond]['Parametric weight'])):
                pred_weights[prt_data[cond]['Time start'][event]:prt_data[cond]['Time stop'][event]+1] = weights[event]
        
            # create convolved parametric predictor
            pred_weights = np.concatenate((np.zeros(len(hrf)-1), pred_weights)) # pad predictor with zeros
            pred_weights_conv = signal.convolve(pred_weights, hrf, 'valid', method='auto')
            # downsample parametric predictor to TR
            f2 = interpolate.interp1d(time_msec, pred_weights_conv)
            pred_weights_conv = f2(time_TR)
            
        
    # if PRT is defined in volumes, HRF is created in volume resolution and the convolution of the boxcar-function and the HRF is performed in volume resolution    
    else:
        # create empty predictor in VOLUME resolution
        pred_empty = np.zeros(fmr_header['NrOfVolumes']).astype(float)                        
        # create block predictor
        pred_block = pred_empty
        for event in range(len(prt_data[cond]['Time start'])):
            # in case of a PRT in volume resolution, the indexing has to change, as volume 1 represents index 0
            pred_block[prt_data[cond]['Time start'][event]-1:prt_data[cond]['Time stop'][event]] = 1
        
        # create convolved predictor
        pred_block = np.concatenate((np.zeros(len(hrf)-1), pred_block)) # pad predictor with zeros
        pred_conv = signal.convolve(pred_block, hrf, 'valid', method='auto')

        if add_parametric_predictor:
            pred_weights = np.zeros(fmr_header['NrOfVolumes']).astype(float)
            for event in range(len(prt_data[cond]['Parametric weight'])):
                pred_weights[prt_data[cond]['Time start'][event]-1:prt_data[cond]['Time stop'][event]] = weights[event]
            
            # create convolved parametric predictor
            pred_weights = np.concatenate((np.zeros(len(hrf)-1), pred_weights)) # pad predictor with zeros
            pred_weights_conv = signal.convolve(pred_weights, hrf, 'valid', method='auto')

            
### Define Single Study Design Matrix based on the created predictors
    
    ## Add main predictor to model
    # scale predictor to 1    
    pred_conv = pred_conv/max(pred_conv)
    # copy to design
    sdm_data.append({})
    sdm_data[-1]['NameOfPredictor'] = prt_data[cond]['NameOfCondition']
    sdm_data[-1]['ColorOfPredictor'] = prt_data[cond]['Color'].tolist()
    sdm_data[-1]['ValuesOfPredictor'] = pred_conv
    del(pred_empty, pred_block, pred_conv)
    
    ## Add parametric predictor to model
    if add_parametric_predictor:
        # scale predictor to 1    
        pred_weights_conv = pred_weights_conv/max(pred_weights_conv)
        # copy to design
        sdm_data.append({})
        sdm_data[-2]['NameOfPredictor'] = prt_data[cond]['NameOfCondition'] + ' [Main]' # change name of main predictor
        sdm_data[-1]['NameOfPredictor'] = prt_data[cond]['NameOfCondition'] + ' [Parametric]'
        sdm_data[-1]['ColorOfPredictor'] = prt_data[cond]['Color'].tolist()
        sdm_data[-1]['ValuesOfPredictor'] = pred_weights_conv
        del(pred_weights, pred_weights_conv, weights)

# Add constant to model
sdm_data.append({})
sdm_data[-1]['NameOfPredictor'] = 'Constant'
sdm_data[-1]['ColorOfPredictor'] = [255,255,255]
sdm_data[-1]['ValuesOfPredictor'] = np.ones(fmr_header['NrOfVolumes']).astype(float)

# Define SDM header
sdm_header = {'FileVersion': 1, 'NrOfPredictors': len(sdm_data), 'NrOfDataPoints': int(fmr_header['NrOfVolumes']),
              'IncludesConstant': 1, 'FirstConfoundPredictor': len(sdm_data)}

# Save SDM
bvbabel.sdm.write_sdm((prt_file[:-4] + '.sdm'), sdm_header, sdm_data)
      