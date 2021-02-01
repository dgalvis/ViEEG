function [out,out_full] = BNI_single(net_aux,out,out_full,count,params)
% Find BNI for a given w = out{count}
% M.A.Lopes, 2017
% Adjustments: D Galvis 2019
% inputs:
% net_aux: where net_aux is a sources x sources array (connectivity matrix)
% out{count} - the interpolated best w value from BNI_find
% out_full{count} - all BNI for attempted w in BNI_find
% count: current count
% params.n_n - number of noise runs
% params.T - time steps
% params.I_0 - initial condition
% params.I_sig - amount of noise
% ----------------------------------------------------------------------- %
% outputs:
% out_full{count}.w_final = w = out{count} from BNI_find
% out_full{count}.BNI_final = actual BNI for the above w
% out_full{count}.BNI_final_full = array ( sources x noise runs)
%                                        (for the best w)
%                                        (mean over dim 1 to get a
%                                        distribution)
% ----------------------------------------------------------------------- %
    % inport parameters
    n_n   = params.n_n;       % number of noise runs
    T     = params.T;         % time steps
    I_0   = params.I_0;       % distance to SNIC (initial conditions)
    I_sig = params.I_sig;     % noise level
    flag= 'BNI';              % proper normalization

    % We only consider undirected networks here
    if issymmetric(net_aux)
        disp('symmetric');
        net = net_aux;
    else
        disp('symmetric now!!');
        net = net_aux + net_aux';
    end

    % w is the best w from BNI_find
    w   = out{count};

    % Calculate the actual exact BNI for that interpolated w value
    % Hopefully this is very close to 0.5
    N = length(net);
    BNI = zeros(N,n_n);
    
    % Seeds so that different noise runs will give different outputs
    % with the parallel processing
    seeds = randi(2^32-1, [n_n, 1]);
    parfor noise= 1:n_n
        BNI(:,noise) = theta_model(net,T,w,I_0,I_sig,flag,seeds(noise));
    end
    
    % Results
    out_full{count}.w_final = w;
    out_full{count}.BNI_final = squeeze(mean(mean(BNI)));
    out_full{count}.BNI_full_final = BNI;
end
