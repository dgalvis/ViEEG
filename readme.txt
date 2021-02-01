ViEEG signals (ictal) reconstruction pipeline:
Descriptions: here we offer the processing pipeline for ViEEG (MEG) signal reconstruction for review. We used Freesurfer6 standard processing pipeline (recon-all) to process patient's individual MRI scans and CURRY 8 ® software to generate boundary element method (BEM) models used as head models and define ViEEG electrode coordinates. Next, MNE-Python pipeline (ViEEG_ictal_signal_reconstruction_pipeline.py) imports head models (in Freesurfer RAS coordinate system), ictal MEG sensor signals and ViEEG electrode coordinates (in Freesurfer RAS coordinate system) to perform ViEEG signal reconstruction using a beamformer algorithm (https://mne.tools/stable/generated/mne.beamformer.apply_lcmv_epochs.html#mne.beamformer.apply_lcmv_epochs) with source signal orientation chosen to maximise the source power. This step utilises the standard MNE-Python MEG source reconstruction pipeling using LCMV beamformer algorithm. Reconstructed ViEEG signals are then saved in matrix format with sufficient supplementory information to record and explain how the data (in both .mat data file and a seperate .json file).

Required sofware and versions:
Operating Systems: Ubuntu 16.04 LTS (Linux x86 64bits)
MNE-Python Toolbox, version: 0.19.0 (Follow the link for installation of MNE-Python Toolbox using Anaconda virtual environment for a complete installation with all required dependencies https://mne.tools/stable/install/mne_python.html)
Python, version: 3.7.4
CURRY 8 ®(Neuroscan Compumedics), version: 8.0.5
Freesurfer, version: 6.0.0 (Follow the link for installation of Freesurfer6 https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)


This code is found in the source reconstruction directory.

Compare results against existing clinical localisation:
The outputs of dynamical network models can be compared against clinical localisation at nodal level. The Excel file (clinical_localisations_for_comparison.xlsx) under folder sample_data/ contains data of existing clinical localisations for patient 1 and 4, as examples. The clinical localisations include resection margin (column resection_margin in Excel file), pathology identified by MRI (column MRI_pathology), seizure onset zone (SOZ) identified by iEEG if applicable (column iEEG_SOZ), abnormality identified by PET (column PET_abnormality), early-phase of MEG source localisation (column early-MSL), middle-phase of MEG source localisation (column mid-MSL), late-phase of MEG source localisation (column late-MSL) and the earliest solution of MEG source localisation and EEG source localisation (column earliest_solution). Value zero/0 means the ViEEG electrode is not in a clinical localisation and value one/1 means the ViEEG electrode is in a clinical localisation. N/A in MRI_pathology or PET_abnormality means a clinical localisation is reported as normal/negative. N/A in iEEG_SOZ means either iEEG recording was not conducted on this patient or clinical iEEG data is not available to us (patients from external surgical centres).


The code to develop functional networks and identify "node ictpogenicity" is found in the modelling directory. There is a separate, detailed README.txt which you can use to rerun the functional connectivity derivations and mathematical modelling. 
