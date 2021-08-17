function NI = NI_model(net_aux, w, params)
% Find BNI for a given w = out{count} after removing 1 node at a time
% M.A.Lopes, 2017
% Adjustments: D Galvis 2019
% inputs:
% net_aux: where net_aux is a sources x sources array (connectivity matrix)
% out{count} - the interpolated best w value from BNI_find
% NI_out - previous results NI_out{1} ... NI_out{count - 1 }
% NI_out_full - previous results NI_out_full{1} ... NI_out_full{count - 1}
% count: current count
% params.n_n - number of noise runs
% params.T - time steps
% params.I_0 - initial condition
% params.I_sig - amount of noise
% ----------------------------------------------------------------------- %
% outputs:
% NI_out{count} = rank order of most to least ictogenic [idx, BNI(idx]
% NI_out_full{count} = array (N-1, n_n, N)
%                      dim 1 BNI of each leftover node
%                      dim 2 noise run
%                      dim 3 removed node
%                      mean over first 2 dims for BNI after node removal
%                      mean over first 1 dims for distribution per node
%                      removal
% ----------------------------------------------------------------------- %
    % inport parameters
    n_n   = params.n_n;       % number of noise runs
    T     = params.T;         % time steps
    I_0   = params.I_0;       % distance to SNIC (initial conditions)
    I_sig = params.I_sig;     % noise level
    flag= 'BNI';              % proper normalization (some use 'NI')
    
    
    
    % We only consider undirected networks here
    if issymmetric(net_aux)
        disp('symmetric');
        net = net_aux;
    else
        disp('symmetric now!!');
        net = net_aux + net_aux';
    end

    N=length(net);    % # nodes
    NI=zeros(N-1,n_n,N);
    
    % Seeds so that different noise runs will give different outputs
    % with the parallel processing
    seeds = randi(2^32-1, [n_n, N]);
    
    % Remove each node one at a time
    for node_R=1:N
        disp(['node #: ',num2str(node_R)]);
        net_auxB = net;
        net_auxB(node_R,:)=[];
        net_auxB(:,node_R)=[];
        % Run n_n times to account for noise and/or get a distribution
        parfor noise=1:n_n
            NI(:,noise,node_R)=theta_model(net_auxB,T,w,I_0,I_sig,flag,seeds(noise,node_R));
            
        end
        
    end


end
