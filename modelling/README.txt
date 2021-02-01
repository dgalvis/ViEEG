REPOSITORY NAME: Virtual_Grid

CODE AUTHOR: Daniel Galvis (except where otherwise noted)

LICENSE:: GPL v3 
    Copyright (C) 2019 Daniel Galvis
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

PURPOSE : This code takes in beamformer reconstructed virtual sources from MEG seizure data, performs functional connectivity analysis on the time series, and uses the resulting networks to model node ictogenicity and infer an optimal resection site.

CITED WORKS : 
. Bonferroni - Holm multiple comparisons correction
    . Holm, S. 1979
    . Groppe et al., 2011
. Surrogates
    . Schreiber and Schmitz 1996 
    . Rummel et al., 2011 
. Mutual Information (MI5 and MI5_UNI, mutual information based least dependent component analysis)
    . Kraskov, H. et al., 2004 
. Phase-locking value
    . Jervis et al., 1983
    . Lachaux et al., 1999
    . Palva et al., 2012 (imaginary version)
. Amplitude-based correlation coefficient
    . Schoffelen et al., 2009 
    . Hipp et al., 2012 (orthogonalized version)
    . Brookes et al., 2012 (orthogonalized version)
. Functional Connectivity Generally
    . Palva et al. 2018
. Brain Network Ictogenicity and Node Ictogenicity / Theta Model 
    . Schmidt et al. 2014
    . Goodfellow et al., 2016 
    . Lopes et al., 2017 
    . Lopes et al., 2018 

. Codes: The following functions were written by someone other than Daniel Galvis
    . MEGpy/MutualInformationICA/*: written by Kraskov H et al., (2004) 
    . MEGpy/connectivity.py function out = MIhigherdim(d1, d2): written by Kraskov, H et al., (2004), edited by Lopes, M (2017, 2018), and written for Python by Daniel Galvis
    . MEGmat/BNI_single.m, MEGmat/BNI_find.m, MEGmat/NI_model.m, MEGmat/theta_model.m: were written by Lopes, M (2017, 2018) and edited for this pipeline by Daniel Galvis
    . MEGmat/bonf_holm.m: was written by Groppe, D (2011)
    * For more details, please refer to those functions.

PUBLICATIONS   : Submitted Nature Comms

DEPENDENCIES : This code uses both MATLAB and python. This code was mean to be run with python3.6.5 (see environmentVG_<os_name>.yml for suggested environment) and MATLAB2017a (see matlab_version_info.mat).

The python environment can be installed using anaconda with the following command: conda env create -f environmentVG_<os_name>.yml
(<os_name> is osx for MAC or linux for linux)

TO RUN: ./run_job.sh (must be executable chmod +x) after dependencies are correct

1) Compiles Kraskov H et al. (2004) mutual information code
./MI_initiator.sh

2) Edit the config files, these are found in subdirectories of net_results. This is where you add the frequency band of interest, the patient filenames (and seizure epoch times), the method used.
(example ./net_results/CC for the CC method with patient 1 and 4, number of surrogates has been reset to 19 to make it run faster in examples)

2b) There are params in ModelRun.m for the dynamical systems modelling you wish to change as well. This is T: how long to run each dynamical system, n_n: how many runs for noise (in manuscript I used 128, for example I use 8 to make it run in reasonable time without HPC).

3) Create Networks using python
ipython NetRun.py ./CC
ipython NetRun.py ./MI
(The config files in net_results/CC and net_results/MI are informative. They are the configurations and include the seizure ranges, bandpass filtering range, functional connectivity matrix)
(The networks are saves into net_results/CC and net_results/MI. They get used by the next step)
(NetRun.py contains a docstring with detailed steps of what is run)

4) Calculates Node Ictogenicity (Run modelling)
matlab -nodisplay -r "ModelRun('./CC');exit;"
matlab -nodisplay -r "ModelRun('./MI')exit;"
(the results of the python code must be completed.)
(ModelRun contains a docstring with detailed steps of what is run)
(Note: The dynamical system is run 128 times for each network, so that the results are not dependent on the noise in the model, but in this example, we run 8 times so it is fast. See below)

SAMPLE RESULTS:
--------
These have already been run. The sample results are included in the repository. Note that results may be slightly different from what was used in the manuscript for three reasons:
1) Random numbers are used in surrogate generation and mathematical modelling
2) 19 surrogates are used for CC instead of 199 in the manuscript (this is so running examples if faster for you)
3) The example data is trimmed down for file size. Running the band pass filtering and resampling may be slightly affected by this.
4) For mathematical modelling, fewer dynamical systems runs are used in the example (8 instead of 128 as in the manuscript).
However, the results have been checked and are very similar despite these things (The top ictogenic nodes predicted are the same, but using fewer noise runs (4) makes the solution more variable). 
The results used for publication are available upon request.

ipython NetRun.py ./CC: Creates .mat files in ./net_results/CC with the same filenames as the patient filenames
                        These contain structs (one for each frequency band in the config file). In example, this is a struct called "a"
                        a contains networks for each method used, in this case just CC, so a.CC is a network

                  ./MI: Same thing (except is is MI5_UNI)
                  Note: MI5_UNI is the Kraskov mutual information

                  In example, these is one of these for patient 1 and patient 4.
                  See the supplemental materials for a description of how the networks were derived (also NetRun.py docstring)

matlab -r "ModelRun('./CC')": Again Creates .mat files in ./analysis_results/CC with the same filenames as the patient filenames
 
This code performs the mathematical modelling. It uses the theta model using the functional network as connections. More details in supplemental materials and docstring

For ./CC, there are files in analysis_results/CC with the patient filenames: output.a.CC contains the useful information in the example.

- The global connectivity parameter is determined for the full network such that BNI = 0.5 (that is so that the average of all nodes is active 50% of the time)
The global coupling value is found in output.a.CC.out, this is the average across 8 runs (128 in manuscript) to get rid of noise.
The full set of 8 values is found in output.a.CC.out_full{1}.BNI_full_final. This is matrix of nodes x number of runs that has the fraction of active time for each one.
(Other contents of output_full{1} is explained in docstrings, but they relate to the iterations to calculate the global coupling value)

- The node ictogenicity value is in output.a.CC.NI_out{1}. For each node, it identifies the new BNI value when it is removed from the network. If BNI is decreased, this implies
the node contributes to the increased activity. Here you can see the rank order of most to least ictogenic nodes.

output.a.CC.NI_out{1}  matrix  number of nodes x 2, first column node ID, second column BNI without that node, sorted from smallest BNI (max NI) to largest BNI (min NI)
output.a.CC.NI_out_full{1} matrix of size (number of nodes - 1) x number of runs x number of nodes, The first two dimensions are the BNI_full of the reduced network. The
third dimension indicates there is one for each node (since every node gets removed one by one).

(The 8 runs are so that the BNI orderings are not dependent on the random variable, see supplemental materials)

output.a.CC.results{1} - This contains a second way of ordering the nodes (which is almost identical - indicating the method is robust). 
                         The docstrings outline what is contained in this struct but they were not used by M. Cao in his final analysis (He used output.a.CC.NI_out{1}). 
                         M.Cao did have access to output.a.CC.results{1}.sig and output.a.CC.results{1}.sigbh, which use a ranksum test to determine the significant nodes.
                         sig (ranksum) and sigbh (with Bonferroni Holms) compares the BNI of the 128 runs with node i removed to the 128 x (num_nodes -1) runs where i neq j is removed.
(same idea for MI)



REQUIRED CONTENT:
--------
INSTALL TIMES:
Installation of anaconda and MATLAB are required. Typical install times
Compiling mutual information code, almost immediate
--------
DEMO TIMES:
Notes: Running the CC measure is much faster than the MI measure! (for the network creation component)

MACBOOK PRO 2019 Catalina:
surrogate creation (19 surrogates) < 1 minute (for 75 - 100 nodes: the code makes number of time points divisible by 10 which helps)
CC (75 - 100 nodes, 19 surrogates) < 1 minute
MI (75 - 100 nodes, 19 surrogates) < 1 hour
(this is for two patients)

Modelling (75 - 100 nodes) 3-6 hours (with 8 repeats for each network, see sample results above) (for two patients)
(This code was originally run on an HPC with 16-32 cores and ~15 jobs running in parallel)

---------
EXPECTED OUTPUTS:
expected network output located in (./net_results/CC and ./net_results/MI)
expected modelling output located in (./analysis_results/CC and ./analysis_results/MI)
---------
DATA SAMPLE:
Located in ./patients directory)
---------
SOFTWARE DEPENDENCIES
MATLAB2017A (see matlab_version_info.mat, load into matlab)
Python (environmentVG_osx.yml or environmentVG_linux.yml - install instructions above)
---------
RUN INSTRUCTIONS:
Listed above
NetRun.py - docstring with instructions for data preprocessing and network creation
ModelRun.mat - docstring with instructions for node ictogenicity quantification
---------

CONFIGURATION FILES (in a subdirectory of net_results):

pats:
    Patient_1_ictal_ViEEG_signals.mat:
        low: [0]
        high: [17]
    Patient_4_ictal_ViEEG_signals.mat:
        low: [0]
        high: [20] # Low and high should have same number of elements. You can add different epochs
    # Add more patients if desired

methods:
    a:
        method: ['CC'] # Other methods can be added here (There are several unused ones in the code from the literature). 
        lf: 1
        hf: 25
    #b: (if you want other frequency bands)

details: 
    freq_new: 512 # resampling frequency 512 Hz
    max_workers: 8 # This gets reset based on your computer anyways
    nsurr: 19 #surrogate number
    nseg: 10 # number of segments to divide the data into
    sub_size: 0.25 # size of the subsegments (0.25*full_segment size). (nseg and sub_size define sliding window)
    TH: True # Use a ranksum test to compare the 10 data subsegment networks to the 19*10 surrogate subsegment networks
    BH: False # Bonferroni-Holms correction (not used)

---------
The supplemental materials explains the sliding windows, surrogate correction, and ranksum test. 
---------
METHODS (see config):

'CC' and 'MI5_UNI' were used in this project

'CC' - amplitude based correlation coefficient
'MI5_UNI' - Kraskov et al. 2004 mutual information
---
Others:

'PLV' - phase locking value (suitable when using a tight frequency band)
'iPLV' - imaginary phase locking value (suitable when using a tight frequency band)
'oCC' - orthogonalised CC (suitable when there are strongly correlated signals)
'MI'  - This is an unsuitable but fast mutual information (do not use this!)
'MI5' - Kraskov et al. 2004 mutual information but with multivariate surrogates (throws out linear correlations)

BH: True #in config - would use a Bonferroni-Holms correction for the ranksum test to compare true networks to surrogate networks
         #this is only suitable when looking for very strong connections

Univariate surrogates: We use these (don't throw out linear correlations) - Schreiber and Shmitz (1996)
Multivariate surrogates: Throw out linear corrlations (only used with 'MI5' the nonlinear method)
