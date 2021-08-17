%=========================================================================%
%All code was created by Daniel Galvis except where otherwise noted 

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

% This code is for the supplemental materials ground truth check on
% beamforming. See each section for details

%=========================================================================%

%% 1) Initialise everything
clear; clc; close all; restoredefaultpath;

addpath('mat_funs');

idx = 1;
load(['results', num2str(idx),'.mat'], 'FC_beam_full', 'FC_full', 'struct_conn_true', 'node_conns');

try
    parpool('local', 18);
catch
    try
        parpool('local')
    catch
        disp('parpool active');
    end
end

% parameters for running the theta model (BNI, NI)
params.T=4*10^6;         % # time steps
params.n_n=8;          % # runs for noise
params.I_0=-1.2;         % distance to SNIC
params.I_sig=5*1.2*0.1;  % noise level

figure(); imagesc(FC_beam_full);title('beamformer');colorbar;
figure(); imagesc(FC_full);title('original');colorbar;
figure(); imagesc(struct_conn_true);colorbar;
pause(0.1);

%% 2) Identify global weight for  BNI = 0.5 
w_beam = BNI_find(FC_beam_full, params);
w_orig = BNI_find(FC_full, params);

%% 3) Identify NI (node ictogenicity)
NI_beam = NI_model(FC_beam_full, w_beam, params);
NI_orig = NI_model(FC_full, w_orig, params);
save(['NI', num2str(idx),'.mat'], 'NI_beam', 'w_beam', 'NI_orig', 'w_orig');



