function ModelRun(indir, outdir)
%All code was created by Daniel Galvis except where otherwise noted (MEGmat/BNI_find.m, MEGmat/BNI_single.m MEGmat/NI_model.m, MEGmat/theta_model.m, MEGmat/bonf_holm.m, MEGpy/MutualInformationICA/*, and MEGpy/connectivity.py->out=MIhigherdim(d1,d2))
%
%Copyright (C) 2019 Daniel Galvis
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <https://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------%
% function ModelRun(indir, outdir)
% note: ModelRun(indir) works as ModelRun(indir, indir)
% ------------------------------------------------------------------------%
% inputs:
% indir  : a subdirectory of ./net_results
% outdir : a subdirectory of ./analysis_results
% indir contains the config file and the functional connectivity results
% outdir will contain files with the same names as in ./net_results/indir
% but with the results of the BNI and NI analysis with the theta model
% ------------------------------------------------------------------------%
% outputs:
% For each functional connectivity network in the directory
% ./net_results/indir, there is an output 
% For each patient/seizure file, there is a file in ./net_results/indir
% with a network for all epochs, frequency bands, connectivity method
% For each patient/seizure file, this will create a file in
% ./net_results/outdir with the analysis of each network
%
% input is of the form a.CC(:,:,1:N) 'a','b'... -freq band, 'CC','MI' ...
% - conn method, N - number of epochs for that 'pat' (see config.ini)
% output is of the from output.a.CC{1} ... a.CC{N} with the analysis for
% that network
% ------------------------------------------------------------------------%
% For each epoch:
% BNI_find.m - a list of attempted weights w and the resulting BNI for it
%          - it continues until it find a w such that BNI ~ 0.5
% BNI_single.m - finds the final w, which interpolates the result of BNI_find
% to try and optimize w such that BNI = 0.5
% NI_model.m - given that w, it removes each node one at a time and
% recalculates BNI. Ictogenic nodes lower BNI when removed
% analysis.m - see ./MEGmat/analysis.mat, but this rank orders the nodes in a
% couple ways and saves a few other interesting results
% ------------------------------------------------------------------------%

    % All extra functions are in ./MEGmat
    addpath('MEGmat');

    % Activate parallel pool (with max number of workers on the machine)
    myCluster = parcluster('local');
    parpool('local', myCluster.NumWorkers);

    % seed the random number generator (seed is saved in the output .mat
    % file)
    rng('shuffle');
    seed = randi(2^32-1);
    rng(seed);
    
    % Inputs will be in ./net_results/indir
    % Outputs will be in ./analysis_results/outdir (but usually we just
    % want this to be in ./analysis_results/indir so allow only one input
    % to ModelRun(indir) = ModelRun(indir, indir)
    try
        outdir = fullfile('analysis_results',outdir);
    catch
        outdir = fullfile('analysis_results',indir);
    end
    indir = fullfile('net_results',indir);
    
    % Iterate over all Patient/seizure files in ./net_results/indir
    files = dir(fullfile(indir,'Patient*.mat'));

    % We will need this temp file to create the output struct
    temp = {};temp.out = {};temp.out_full = {};
    temp.NI_out = {}; temp.NI_out_full = {};temp.results = {};
    

    % parameters for running the theta model (BNI, NI)
    params.T=4*10^6;         % # time steps
    params.n_n=8;          % # runs for noise
    params.I_0=-1.2;         % distance to SNIC
    params.I_sig=5*1.2*0.1;  % noise level

    % Iterate over patients
    for file = files'
        disp(file.name);
        % Create a struct in with the functional networks
        % The functional networks are in the dictionary format from the
        % config file for python (for Patient X)
        % in.a.CC - frequency band a, method CC
        % in.b.MCC - frequency band b, method MCC
        % These arrays are of size (source,source,time_segment_list)
        % Where time_segment_list is the number of segments analysed (see
        % config, low/high in 'pats') for patient X seizure Y file
        in = load(fullfile(indir, file.name));
        % fieldnames a,b,c (frequency bands) etc.
        fn = fieldnames(in);
        % Result will have the same form as variable in
        output = in;
        % Iterate over frequency bands
        for k =1:numel(fn)
            % 'a', 'b' ...
            fr_bd   = in.(fn{k});
            % fieldnames CC, MCC, PLV etc
            fn2 = fieldnames(fr_bd);
            % Iterate over functional connectivity methods
            for k2= 1:numel(fn2)
                % functional connectivity methods CC, MCC, PLV etc.
                f_c_m = fr_bd.(fn2{k2});
                
                % Format output of the BNI, NI modelling
                output.(fn{k}).(fn2{k2}) = temp;
                out = {}; out_full = {}; NI_out = {}; NI_out_full = {}; results = {};               
                % Iterate over analysed time segments
                for k3= 1:size(f_c_m,3)
                    % connectivity matrix for the given time segment,
                    % functional connectivity method, frequency band,
                    % patient/seizure file
                    network = f_c_m(:,:,k3);
                    
                    % Find the weight w, such the the connectivity matrix
                    % gives a BNI = 0.5. This finds a close w value on
                    % either side
                    [out, out_full] = BNI_find(network, out, out_full, k3, params);
                    % Linear interpolation of the results from BNI_find to
                    % find a weight w that is hopefully very close to
                    % giving BNI = 0.5
                    [out, out_full] = BNI_single(network, out, out_full, k3, params);
                    
                    % Calculate NI. Remove each source and recalculate BNI
                    % to see if it increases or decreases
                    [NI_out, NI_out_full] = NI_model(network, out, NI_out, NI_out_full, k3, params);
                    % Rank order the nodes in terms of ictogeniticy and
                    % find the significantly ictogenic nodes (see analysis
                    % for all the fun stuff it does)
                    [results] = analysis(results, out_full, NI_out_full, k3);
                end
                % Put the result in output, which has a nearly identical
                % form to the input in ./net_results/indir
                % output.a.CC{1} is modelling results for in.a.CC(:,:,1)
                % output.b.CC{2} is modelling results for in.b.CC(:,:,2)
                % etc.
                output.(fn{k}).(fn2{k2}).out = out;
                output.(fn{k}).(fn2{k2}).out_full = out_full;
                output.(fn{k}).(fn2{k2}).NI_out = NI_out;
                output.(fn{k}).(fn2{k2}).NI_out_full = NI_out_full;
                output.(fn{k}).(fn2{k2}).results = results;
                % save the seed as well for posterity
                % as with the network results, each Patient/Seizure in
                % 'pats' has it's own file
                save(fullfile(outdir,file.name),'output','seed','params');
                
            end
        end
    end
end
