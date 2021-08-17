function [out, seed] =theta_model_time_series(network,T,w,I_0,I_sig,flag,seed)

    if nargin < 7
        rng('shuffle');
        seed = rng;
    else
        rng(seed);
        seed = rng;
    end
    
    if nargin < 6
        disp('Missing arguments. Please read the instructions.');    
        out=[];
        return;
    end
    
    dt=10^-2;  % time step
    threshold=0.9; % threshold for BNI
    window_epochs=6*4/dt; % window for BNI

    I_sig=I_sig/sqrt(dt); 
    N=length(network); % number of nodes
    
    % normalisation of coupling
    if strcmpi(flag,'BNI')
        wnet=w*network/N;
    elseif strcmpi(flag,'NI')
        wnet=w*network/(N+1);
    else
        disp('The flag argument should be either BNI or NI.');  
        BNI=[];
        return;
    end
    wnet=wnet';
    
    if length(I_0)==1
        I_0=I_0*ones(N,1);
    end
    if length(I_0)~=N
        disp('Please ensure I_0 is a vector Nx1.');   
        out=[];

        return;
    end
    D=size(I_0);
    if D(1)<D(2)
        I_0=I_0';
    end
    signal=zeros(T,N);
    theta_s=-real(acos((1+I_0)./(1-I_0))); % stable point if I_0 < 0
    theta_old=theta_s; % initial condition  

    % Compute time series
    for time=1:T-1
        I=I_0+I_sig*randn(N,1)+wnet*(1-cos(theta_old-theta_s));
        theta_new=theta_old+dt*(1-cos(theta_old)+(1+cos(theta_old)).*I);
        signal(time+1,:)= 0.5*(1-cos(theta_old-theta_s));
        theta_old=theta_new;
    end
    
    out = signal;
end