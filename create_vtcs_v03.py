''' create_vtcs_v03.py 

    purpose: script for creating vtc files in mni space. this script walks through subdirectories; if all 5 filetypes (vmr, fmr, ia, fa, mni) are found, a VTC file is created in that subdirectory.
    usage: click 'Python' icon in BrainVoyager 23, open this script. adapt manual settings if required. press the 'run' arrow.
    
    tested in BrainVoyager 23.0.9 on macOS 12.5 (Apple Silicon) with 1 nested directory level 
    not tested yet: 1. unicode paths # os.walk(u".") 2. other operating systems 3. filenames with spaces 
    
    240219hb
'''
 
import os

#------- manual settings -------
resolution      = 2
interpolation   = 2  # 0 -> nearest neighbour, 1 -> trilinear, 2 -> sinc
exttal          = False
datapath        = brainvoyager.choose_directory("Please select the top directory") # brainvoyager.sampledata_path + 'GSGData'
# ------------------------------

for currentdir, subdirs, filenames in os.walk(datapath):
    vtcfiles = dict() # {'VMR', 'FMR', 'IA', 'FA', 'MNI'} 
    brainvoyager.print_to_log('Reset files...')
    brainvoyager.print_to_log('Current directory: ' + currentdir)
    for filename in filenames:           
        print('Current directory: ' + currentdir) 
        if filename.endswith('.fmr') and not filename.endswith('_firstvol.fmr'): 
            vtcfiles['FMR'] = os.path.join(currentdir, filename)
            brainvoyager.print_to_log('FMR found: ' + vtcfiles['FMR'])
        if filename.endswith('.vmr'):     
            vtcfiles['VMR'] = os.path.join(currentdir, filename)
            brainvoyager.print_to_log('VMR found: ' + vtcfiles['VMR'])
        if filename.endswith('_IA.trf'):     
            vtcfiles['IA'] = os.path.join(currentdir, filename)     
            brainvoyager.print_to_log('Initial alignment file found: ' + vtcfiles['IA'])       
        if filename.endswith('_FA.trf'):     
            vtcfiles['FA'] = os.path.join(currentdir, filename)
            brainvoyager.print_to_log('Fine alignment file found: ' + vtcfiles['FA']) 
        if filename.endswith('MNI_a12.trf'):     
            vtcfiles['MNI'] = os.path.join(currentdir, filename)
            brainvoyager.print_to_log('MNI normalisation file found: ' + vtcfiles['MNI']) 
        if (len(vtcfiles.values()) == 5):    
            brainvoyager.print_to_log('Number of unique file types found: ' + str(len(vtcfiles)))    
            brainvoyager.print_to_log('All files present, proceed with creating VTC...')
            doc = brainvoyager.open_document(vtcfiles['VMR'])
            newname = os.path.join(currentdir, vtcfiles['FMR'].split(".fmr", 1)[0]+'_res'+str(resolution)+'_interp'+str(interpolation)+"_MNI.vtc")
            doc.vtc_creation_extended_tal_space = exttal
            brainvoyager.print_to_log('Creating ' + newname + ' from:' +'\n' +'\t'+vtcfiles['FMR'] + '\n' +'\t'+ vtcfiles['IA'] + '\n'+ '\t'+ vtcfiles['FA'] + '\n' +'\t'+ vtcfiles['MNI'])
            ok = doc.create_vtc_in_mni_space(vtcfiles['FMR'], vtcfiles['IA'], vtcfiles['FA'], vtcfiles['MNI'], newname, resolution, interpolation)
            doc.close() 
            vtcfiles = dict() # reset   
       
