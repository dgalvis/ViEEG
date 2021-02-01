function [results] = analysis(results,out_full,NI_out_full,count)
%All code was created by Daniel Galvis except where otherwise noted (MEGmat/BNI_model.m, MEGmat/BNI_single.m MEGmat/NI_find.m, MEGmat/theta_model.m, MEGmat/bonf_holm.m, MEGpy/MutualInformationICA/*, and MEGpy/connectivity.py->out=MIhigherdim(d1,d2))
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
% function [results, count] = analysis(results, out_full, NI_out_full, count)
% inputs:
% results: results{1} ... results{count - 1} 
% out_full{count}
% NI_out_full{count} = array (N-1, n_n, N)
%                      dim 1 BNI of each leftover node
%                      dim 2 noise run
%                      dim 3 removed node
%                      mean over first 2 dims for BNI after node removal
%                      mean over first 1 dims for distribution per node
%                      removal
% out_full{count}.w_save - all attempted w values
% out_full{count}.BNI = array ( sources x noise runs x iterations)
%                             (mean over first 2 dimensions is BNI)
% out_full{count}.w_final = w = out{count} from BNI_find
% out_full{count}.BNI_final = actual BNI for the above w
% out_full{count}.BNI_final_full = array ( sources x noise runs)
%                                        (for the best w)
%                                        (mean over dim 1 to get a
%                                        distribution)
% count: current count
% ----------------------------------------------------------------------- %
% outputs:
% results{count}.full_order = NI_out{count} = rank order from most to least
% ictogenic [idx, BNI(idx)]
% results{count}.sig, results{count}.sigbh - see signif_analysis.m
% results{count}.BNI_pd, results{count}.BNI_ci - see how close to BNI = 0.5
% our choice of best w actually got (in BNI_find.m/BNI_single.m)
% results{count}.new_order - a different way of rank ordering most to least
% ictogenic nodes based on area under the curve of fraction of the noise
% runs that node i in in the top 1 ... num sources when each is sorted for
% the original ictogenicity (see below)
% ----------------------------------------------------------------------- %
% Code Author: Daniel Galvis, 2019

    % Calculate the most to least ictogenic nodes
    % and the BNI values after they are removed
    x = NI_out_full{count};
    x = mean(mean(x,2),1);
    x = squeeze(x);
    [~,idx] = sort(x);
    results{count}.full_order = [idx,x(idx)];

    % use a wilcoxon ranksum test of BNI with node i removed vs. random BNI
    % If BNI with node i removed is significantly less, then ictogenic
    % with and without bonferroni holms correction
    [results{count}.sig, results{count}.sigbh] = signif_analysis(NI_out_full, count);
    
    % Calculate the distribution of BNI values and confidence intervals
    % for the original BNI calculation of the whole network
    % (longer theta runs T>4*10^6 could reduce this deviation)
    x = squeeze(mean(out_full{count}.BNI_full_final,1))';
    pd = fitdist(x,'Normal');
    ci = paramci(pd);  
    results{count}.BNI_pd = pd;
    results{count}.BNI_ci = ci;

    % A second way of calculated most to least ictogenic nodes
    x = NI_out_full{count};
    % This contains BNI for each noise run (n_n) for each removed node i
    % source number vs n_n array after transpose taken
    x = squeeze(mean(x,1))';
    % sort each column and keep idxs 
    % for each column most ictogenic (smallest BNI) is at the top
    % idx is the source numbers in that order
    % The final product is an NI ordering for each noise run (n_n), they
    % should be mostly the same but some changes in order due to noise will
    % occur
    [~,idx] = sort(x);
    
    % rank_mat(i,j) - array source_number x source_number
    % rank_mat(i,j) - for source j, how many of the n_n noise runs was it ranked
    % in the top i of ictogenic nodes
    % rank_mat(:,j) - a curve that always goes to n_n as j goes to source number
    rank_mat = zeros(size(x,1),size(x,1));
    % iterate over ranked in the top i of ictogenic nodes
    for i = 1:size(x,1)
            % keep nodes ranked in top i over all n_n runs
            aux = idx(1:i,:);aux = aux(:);
            % for each node find how many times it was ranked in top i
            % nodes over all n_n noise runs
            rank_mat(i,:) = sum(((1:size(x,1)) == aux),1);
    end
    results{count}.rank_mat = rank_mat;
    results{count}.idx = idx;
    
    % A simple integral overeach of the curves (curves normalized so that x
    % and y axes go between 0 and 1)
    % Basically this is right hand rule integral if you assume that for
    % each node the number of times it is in to top 0 of nodes is 0 and
    % include that as the leftmost point
    areas = sum(rank_mat,1)/(size(x,1)*size(x,2));
    
    % In this case, the nodes with the largest area are the most ictogenic
    % because this just means they were in the top i of ictogenicity for
    % most values of i and for most of the noise runs
    [~,idx] = sort(areas,'descend');
    % Save the ordering and the area values
    results{count}.new_order=[idx',areas(idx)'];
    
end

