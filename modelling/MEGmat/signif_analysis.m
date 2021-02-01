function [sig, sigbh] = signif_analysis(NI_out_full, count)
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
% inputs:
% NI_out_full{count} = array (N-1, n_n, N)
%                      dim 1 BNI of each leftover node
%                      dim 2 noise run
%                      dim 3 removed node
%                      mean over first 2 dims for BNI after node removal
%                      mean over first 1 dims for distribution per node
%                      removal
% count: current count
% ----------------------------------------------------------------------- %
% output:
% sig: significant nodes. using a ranksum test on the distribution for node i
% BNI being smaller than for randomly removing a node
% sigbh: bonferroni holms corrected significant nodes
% ----------------------------------------------------------------------- %
    % x(i,j) - BNI calculated on noise run i when node j is removed
    x = NI_out_full{count};
    x = squeeze(mean(x,1));
    
    % Calculate p values of BNI when node j is removed vs. a random node
    % Only keep if the BNI is significantly smaller for removal of j
    for j = 1:size(x,2)
    [p(j),h(j)] = ranksum(x(:,j),x(:),'tail','left');
    end
    % Bonferroni holdms correction
    [pbh, hbh] = bonf_holm(p,0.05);

    % Indices of the significant nodes (which sources are significant)
    sig = 1:size(x,2);
    sig = sig(h==1);
    p_new = p(sig);
    [~,idx] = sort(p_new);
    sig = sig(idx);
    %disp('significant nodes');
    %disp(sig)

    % Indices of the Bonferroni-holms significant nodes (sources)
    sigbh = 1:size(x,2);
    sigbh = sigbh(hbh==1);
    p_newbh = pbh(sigbh);
    [~,idxbh] = sort(p_newbh);
    sigbh = sigbh(idxbh);
    %disp('BH significant nodes');
    %disp(sigbh)
end
