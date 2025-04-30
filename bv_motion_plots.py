"""
Create motion plots similar to the ones created by BrainVoyager, using the *3DMC.sdm files
and prints the maximum motion
"""


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import bvbabel
import os.path
import sys
from PyQt5.QtWidgets import QFileDialog, QApplication

app = QApplication(sys.argv)

# change default plotting properties
mpl.rcParams['axes.labelsize'] = 'small'
mpl.rcParams['xtick.labelsize'] = 'small'
mpl.rcParams['ytick.labelsize'] = 'small'
mpl.rcParams['axes.titlesize'] = 'medium'

# set plotdisp to True if you would like to see the plots while running the script
plotdisp = False

motionfiles, _ = QFileDialog.getOpenFileNames(None, 'Select *3DMC.sdm files', os.getcwd(), 'Motion Files (*3DMC.sdm)')

for f in range(len(motionfiles)):
    sdm_header, sdm_data = bvbabel.sdm.read_sdm(motionfiles[f])
        
    motion = np.empty((0,sdm_header['NrOfDataPoints']))
    for pred in range(len(sdm_data)):
        motion = np.vstack([motion, sdm_data[pred]['ValuesOfPredictor']])
                                 
  
# Plot the Motion
    plt.style.use('dark_background')
    fig_motion = plt.figure('Motion', figsize=(10, 8))
    ax_motion = fig_motion.add_subplot(111)
    plt.subplots_adjust(right=0.75)
    ax_motion.set_ylim(np.floor(np.amin(motion)).astype(int), np.ceil(np.amax(motion)).astype(int))
    ax_motion.set_yticks(range(np.floor(np.amin(motion)).astype(int),np.ceil(np.amax(motion)).astype(int)+1))
    ax_motion.grid(axis = 'y')
    ax_motion.spines['right'].set_visible(False)
    ax_motion.spines['top'].set_visible(False)
    
    ax_motion.set_title('Rigid Body Motion Parameters - 3 Translations, 3 Rotations')
    
    # plot the motion
    plot_transx, = ax_motion.plot(motion[0,:], linewidth=1, color=np.array(sdm_data[0]['ColorOfPredictor'])/255, label=sdm_data[0]['NameOfPredictor'])
    plot_transy, = ax_motion.plot(motion[1,:], linewidth=1, color=np.array(sdm_data[1]['ColorOfPredictor'])/255, label=sdm_data[1]['NameOfPredictor'])
    plot_transz, = ax_motion.plot(motion[2,:], linewidth=1, color=np.array(sdm_data[2]['ColorOfPredictor'])/255, label=sdm_data[2]['NameOfPredictor'])
    plot_transx, = ax_motion.plot(motion[3,:], linewidth=1, color=np.array(sdm_data[3]['ColorOfPredictor'])/255, label=sdm_data[3]['NameOfPredictor'])
    plot_transy, = ax_motion.plot(motion[4,:], linewidth=1, color=np.array(sdm_data[4]['ColorOfPredictor'])/255, label=sdm_data[4]['NameOfPredictor'])
    plot_transz, = ax_motion.plot(motion[5,:], linewidth=1, color=np.array(sdm_data[5]['ColorOfPredictor'])/255, label=sdm_data[5]['NameOfPredictor'])
        
    handles, labels = ax_motion.get_legend_handles_labels()
    
    fig_motion.legend(handles= handles,loc="center right", frameon=False, framealpha = 1, fontsize = 'small')
    
    outpath = motionfiles[f].split('.')[0] + '_MotionPlot.png'
    fig_motion.savefig(outpath, dpi=600, format='png')
    
    print('')
    print('Input file:', motionfiles[f])
    print('Saved motion plot to:', outpath)
    
    if plotdisp == True:
        plt.show()
    plt.clf()
    
    #subtract first timepoint from all 6 timeseries, so that they start at 0
    motion_run = motion[:,0:] - motion[:,0,None]
    
    print('Maximum Motion with Respect to Reference Run:', np.max(np.absolute(motion)), ',', sdm_data[np.argmax(np.amax(np.absolute(motion) ,1))]['NameOfPredictor'])
    print('Maximum Motion within the Run:', np.max(np.absolute(motion_run)), ',', sdm_data[np.argmax(np.amax(np.absolute(motion_run) ,1))]['NameOfPredictor'])
    print('Maximum Motion Range within the Run:', np.max(np.amax(motion ,1) - np.amin(motion ,1)), ',', sdm_data[np.argmax(np.amax(motion ,1) - np.amin(motion ,1))]['NameOfPredictor'])
    

# Restore the default plotting parameters after script is finished
mpl.rcParams.update(mpl.rcParamsDefault)