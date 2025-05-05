# scripts_bids_derivatives
BrainVoyager Python Scripts
This set of simple Python scripts is intended for use with BrainVoyager 
and the Python library bvbabel (https://github.com/ofgulban/bvbabel). 
The scripts are optimized for the file naming and folder structure convention 
supported by the BrainVoyager Data Analysis Manager 
(https://download.brainvoyager.com/bv/doc/UsersGuide/DataManagementAndWorkflows/Overview.html), 
and they follow the recommendations of the BIDS standard (https://bids.neuroimaging.io/index.html).
While optimized for this structure, the scripts can be adapted for other folder organizations 
and naming conventions.

## 🧪 Tested Dataset
Most scripts were tested using the fMRI NeWbi dataset (https://www.newbi4fmri.com/):
Culham, J., Stubbs, K., Jackson, E., & Lagace Cusiac, R. (2020). newbi4fmri2020 Localizer. OpenNeuro.
DOI: 10.18112/openneuro.ds003433.v1.0.1, accessed 31 January 2023.


## 📁 Example Data Structure
```plaintext
/Users/user/Documents/BrainVoyager/Projects/\
└── NeWBI4fMRI2020/\
   ├── dataset_description.json
   ├── sourcedata/
   ├── rawdata/
   │    └── sub-01/
   │        └── ses-01/
   │            ├── anat/
   │            │   ├── sub-01_ses-01_T1w.json
   │            │   └── sub-01_ses-01_T1w.nii.gz
   │            └── func/
   │                ├── sub-01_ses-01_task-Localizer_run-01_bold.json
   │                ├── sub-01_ses-01_task-Localizer_run-01_bold.nii.gz
   │                └── sub-01_ses-01_task-Localizer_run-01_events.tsv
   └── derivatives/
       ├── rawdata_bv/
       │    └── sub-01/
       │        └── ses-01/
       │            ├── anat/
       │            │   ├── sub-01_ses-01_T1w.vmr
       │            │   └── sub-01_ses-01_T1w.v16
       │            └── func/
       │                ├── sub-01_ses-01_task-Localizer_run-01_bold.fmr
       │                ├── sub-01_ses-01_task-Localizer_run-01_bold.stc
       │                └── sub-01_ses-01_task-Localizer_run-01_bold.prt (or sub-01_ses-01_Localizer_run-1.prt in the example scripts)
       └── workflow_id-1_type-2_name-anat-preprocessing/
           └── sub-01/
               └── ses-01/
                   └── anat/
                       ├── sub-01_ses-01_T1w_IIHC.vmr
                       └── sub-01_ses-01_T1w_IIHC.v16
       
