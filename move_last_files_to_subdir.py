''' Simple script to move the last 5 dicom (noise volumes) files to a subfolder '''
 
import os
import shutil
from glob import glob

datapath = brainvoyager.choose_directory("Please select the top directory") # brainvoyager.sampledata_path + 'GSGData'
            
def last_8chars(x):
    return(x[-8:])
        
for currentdir, subdirs, filenames in os.walk(datapath):

    file_list2 = glob(os.path.join(currentdir, "*.dcm"))    
    sortedlist = sorted(file_list2, key = last_8chars)     
    listlen = len(sortedlist)
    brainvoyager.print_to_log('Length of list list in directory ' + currentdir + ': ' + str(listlen))
    if listlen > 10:
        tmpdir = os.path.join(currentdir, 'last5vols')
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
            brainvoyager.print_to_log('Created directory: ' + tmpdir)    
        for noise_volumes in range(listlen-5, listlen):    
            brainvoyager.print_to_log('Moving ' + sortedlist[noise_volumes])    
            shutil.move(sortedlist[noise_volumes], tmpdir) 