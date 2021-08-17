#Copyright (C) 2019 Daniel Galvis
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.   

Packages:   

- MATLAB2019B
- Python (Anaconda)
- C (MILCA package)

Before running:  
- Install environment   
conda env create -f environment.yml  
pip install mne==0.22.1  

- Compile Mutual Information Code (make sure the file has executable permissions):   
./MI_initiator.sh

Notes:
- .fif and .pom files are not in the repository but available upon request (files too large)  
- You must be in the conda environment (mne) when you run the matlab code (as the code uses terminal commands to run the beamformer   


To Run (in matlab):  
- run_functional_connectivity.m (get time series, run beamformer, run Mutual information, get FC)
- run_ictogenicity.m (find seizure nodes with functional connectivity from previous)  
- run_figures.m to make the figures  

Note:
- results{idx}.mat and NI{idx}.mat contain all the results
- idx = 1 central example, idx = 2 upper left example, idx = 3 lower right example. See run_figures.m to load in the results (time series, functional connectivity matrices, node ictogenicities etc)   
