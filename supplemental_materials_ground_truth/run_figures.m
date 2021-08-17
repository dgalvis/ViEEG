%% Choose the idx: (results{idx}.mat and NI{idx}.mat must exist already)
clc;clear;close all; restoredefaultpath;
idx = 3;

%% Results from run_functional_connectivity.m
load(['results', num2str(idx),'.mat'], ...
    'FC_full', 'FC_beam_full', 'struct_conn', 'struct_conn_true', ...
    'seeds', 'times', 'data', 'data_beam', 'node_conns','w');
num_nodes = length(FC_full);
not_node_conns = setdiff(1:num_nodes, node_conns);

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



figure('Position',[10,10,900,1500]);
subplot(311);
imagesc(struct_conn);
title('structural connectivity');
subplot(312);
imagesc(FC_full);
title('original functional connectivity');
subplot(313);
imagesc(FC_beam_full);
title('beamformed functional connectivity');

%% Results from run_ictogenicity

load(['NI', num2str(idx),'.mat'], 'NI_beam', 'NI_orig');
NI_beam = (0.5 - NI_beam) ./ 0.5;
NI_orig = (0.5 - NI_orig) ./ 0.5;

NI_beam = squeeze(median(median(NI_beam, 2), 1));
NI_orig = squeeze(median(median(NI_orig, 2), 1));
hotspots_beam = zeros(sqrt(length(FC_beam_full)));
hotspots_orig = zeros(sqrt(length(FC_full)));
ct = 1;
for i = 1:8
    for j = 8:-1:1
        hotspots_beam(i,j) = NI_beam(ct);
        hotspots_orig(i,j) = NI_orig(ct);
        %pause(0.1)
        ct = ct + 1;
    end
end

figure();
subplot(223);hold all;
title('beamformed grid 8 x 8');
imagesc(hotspots_beam);colorbar;
axis tight;

subplot(221);hold all;
title('original grid 8 x 8');
imagesc(hotspots_orig);colorbar;
axis tight;

subplot(224);hold all;
title('beamformed histogram 8 x 8');
histogram(NI_beam(not_node_conns), 10,'normalization', 'probability');
histogram(NI_beam(node_conns), 10,'normalization', 'probability');
ylabel('fraction of nodes')
xlabel('value of NI^{i}');
legend('NI not ict', 'NI ict');

subplot(222);hold all;
title('original histogram 8 x 8');
histogram(NI_orig(not_node_conns), 10,'normalization', 'probability');
histogram(NI_orig(node_conns), 10,'normalization', 'probability');
legend('NI not ict', 'NI ict');
ylabel('fraction of nodes')
xlabel('value of NI^{i}');



