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

%% 1) Initialise everything (clear stuff, seeds, parpool)
clear; close all;clc;restoredefaultpath;
idx = 3;
rng(idx);
seeds = randi(2^32-1, [4,1]);

addpath('mat_funs');

try
    parpool('local', 18);
catch
    try
        parpool('local');
    catch
        disp('parpool active');
    end
end

%% 2) Generate Structural Connectivity for Theta Model Connectivity
num_nodes = 64;

% Connections between non-seizure nodes
conn11 = 0.0;
% Connections between seizure and non-seizure nodes
conn12 = 0.0;
% Connections between seizure nodes
conn22 = 1.0;

% Select seizure nodes
if idx == 1
    node_conns = [28,29,36,37]; % Middle
elseif idx == 2
    node_conns = [37,38,45,46]; % Top Left
elseif idx == 3
    node_conns = [19,20,27,28]; % Bottom Right
else
    error('This idx has no node choice. See section (2)');
end

not_node_conns = setdiff(1:num_nodes, node_conns);

% Structural Connectivity (All connection values)
struct_conn = zeros(num_nodes);
% Structural Connectivity bool (connections between seizure nodes)
struct_conn_true = zeros(num_nodes);

% Assign Structural Connectivity values
for i = 1:(num_nodes-1)
    for j = (i+1):num_nodes
        if sum(i==node_conns) && sum(j==node_conns)
            struct_conn(i,j) = conn22;
            struct_conn_true(i,j) = 1;
        elseif sum(i==node_conns) || sum(j==node_conns)
            struct_conn(i,j) = conn12;
        else
            struct_conn(i,j) = conn11;
        end
    end
end

% Symmetric! Struct Conn
struct_conn = struct_conn + struct_conn';
struct_conn_true = struct_conn_true + struct_conn_true';

% Number of connections (between seizure nodes)
num_conns = sum(struct_conn_true(:)) / 2;


figure();
imagesc(struct_conn);
title('Structural Connectivity');
xlabel('node i');
ylabel('node j');
pause(0.1);

%% 3) Parameters for the Theta Model
params.T=4*10^6;         % # time steps
params.n_n=8;            % # runs for noise
%params.I_0=-1.2;         % distance to SNIC
params.I_sig=5*1.2*0.1;  % noise level


% Non-seizure nodes are less excitable
params.I_0(not_node_conns) = -3;
% Seizure nodes are more excitable
params.I_0(node_conns) = -1.2;

% Global connectivity weight
w = 200;

%% 4) Run the theta model to find time series (for doing functional connectivity)

time_secs = 90;
T = 1000 * time_secs; % 90 seconds sampled at 1000 Hz
I_0 = params.I_0;
I_sig = params.I_sig;
flag = 'BNI';

% Get the time series
out = theta_model_time_series(struct_conn, T,w,I_0,I_sig,flag, seeds(2));

% The time points associated with the time series
times = 0:0.001:time_secs;
times = times(1:end-1);

%% 5) Save the information for the beamformer
sz_times = times;
sz_ts = out';
seizure_range = [sz_times(1), sz_times(end)];
seizure_id = 1;
sfreq = 1000;
save(fullfile('ViEEG_forward_inverse_simulations',['model_output', num2str(idx),'.mat']), ...
    'sz_times', 'sz_ts', 'seizure_range', 'seizure_id', 'sfreq');

%% 6) Terminal command for .py and import the beamformed time series
system(['cd ViEEG_forward_inverse_simulations;ipython ViEEG_manuscript_forward_inverse_simulation_Beamforming.py ', num2str(idx), ';cd ..']);
load(fullfile('ViEEG_forward_inverse_simulations',['model_beamform', num2str(idx),'.mat']), 'data_beam');

data_beam = data_beam(10001:end,:);
%% 7) The original time series


data = out(10000:end,:);

times = times(10000:end);
times = times - times(1);

figure();
subplot(221);
plot(times, data(:, node_conns));
axis tight;
xlabel('times (au)');
ylabel('time series seizure nodes');
title('original time series (seizure nodes)');
subplot(222);
plot(times, data(:, not_node_conns));
axis tight;
xlabel('times (au)');
ylabel('time series non-seizure nodes');
title('original time series (non-seizure nodes)');
subplot(223);
plot(times, data_beam(:, node_conns));
axis tight;
xlabel('times (au)');
ylabel('time series seizure nodes');
title('beamformed time series (seizure nodes)');
subplot(224);
plot(times, data_beam(:, not_node_conns));
axis tight;
xlabel('times (au)');
ylabel('time series non-seizure nodes');
title('beamformed time series (non-seizure nodes)');
pause(0.1);


%% 8) Normalise data for Mutual Information
data_beam = (data_beam - mean(data_beam, 1)) ./ std(data_beam, 1);
data = (data-mean(data,1)) ./ std(data, 1);
%% 9) Functional Connectivity (original time series)

% Setup for parpool
ct = 0;
vals = zeros(num_nodes*(num_nodes-1) / 2, 2);
for i = 1:(num_nodes-1)
    for j = (i+1):num_nodes
        ct = ct + 1;
        vals(ct,:) = [i,j];
    end
end

% Run Mutual Information for each node pair
rng(seeds(3));
rand_idx = randi(2^32,num_nodes*(num_nodes-1) / 2,1);
FC_aux = zeros(num_nodes*(num_nodes-1)/2,1);
parfor i = 1:ct
    FC_aux(i) = MIhigherdim([data(:, vals(i,1)), data(:,vals(i,2))], 3,1,1, rand_idx(i));
end

% Turn back from a list to a matrix
FC = zeros(num_nodes, num_nodes);
for i = 1:ct
    FC(vals(i,1), vals(i,2)) = FC_aux(i);
end

% Symmetric
FC = FC + FC';
FC_full = real(FC);

figure();
imagesc(FC_full);
pause(0.1);



%% 10) Functional Connectivity Beamformed

% Setup for parpool
ct = 0;
vals = zeros(num_nodes*(num_nodes-1) / 2, 2);
for i = 1:(num_nodes-1)
    for j = (i+1):num_nodes
        ct = ct + 1;
        vals(ct,:) = [i,j];
    end
end

% Run Mutual Information for each node pair
rng(seeds(4));
rand_idx = randi(2^32,num_nodes*(num_nodes-1) / 2,1);
FC_aux = zeros(num_nodes*(num_nodes-1)/2,1);
parfor i = 1:ct
    FC_aux(i) = MIhigherdim([data_beam(:, vals(i,1)), data_beam(:,vals(i,2))], 3,1,1, rand_idx(i));
end

% Turn back from a list to a matrix
FC_beam = zeros(num_nodes, num_nodes);
for i = 1:ct
    FC_beam(vals(i,1), vals(i,2)) = FC_aux(i);
end

% Symmetric
FC_beam = FC_beam + FC_beam';
FC_beam_full = real(FC_beam);

figure();
imagesc(FC_beam_full)
pause(0.1);

%% 11) Save Results
save(['results', num2str(idx),'.mat'], ...
    'FC_full', 'FC_beam_full', 'struct_conn', 'struct_conn_true', ...
    'seeds', 'times', 'data', 'data_beam', 'node_conns','w');



        

