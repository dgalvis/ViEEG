function w = BNI_find(net_aux, params)
% Find coupling at which BNI~0.5 using the theta model 
% M.A.Lopes, 2017
% Adjustments: D Galvis 2019
% inputs:
% net_aux: where net_aux is a sources x sources array (connectivity matrix)
% out: the previous output out{1} ... out{count - 1} filled in 
% out_full: the previous output out_full{1} ... out_full{count - 1} filled
% in
% count: current count
% params.n_n - number of noise runs
% params.T - time steps
% params.I_0 - initial condition
% params.I_sig - amount of noise
% ----------------------------------------------------------------------- %
% outputs
% out{count} = w such that BNI = 0.5
% out_full{count}.w_save - all attempted w values
% out_full{count}.BNI = array ( sources x noise runs x iterations)
%                             (mean over first 2 dimensions is BNI)
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
           

    w=25;             % initial guess of coupling 
    n_max=10;         % max number of attempts
    crit=0.05;        % criteria for the tuning
    BNI_ref=0.5;      % reference of BNI we are looking for
    displac=0.01;     % displacement to help find BNI_ref

    N=length(net);    % # nodes

    % Auxilliary variables
    it=0;
    z=1;
    x1=0;
    x2=0;
    
    % These will hold the results for all attempted w values
    BNI=zeros(N,n_n,n_max);
    w_save=zeros(n_max,1);

    
    while(z)
        it=it+1;
        % Seeds so that different noise runs will give different outputs
        % with the parallel processing
        seeds = randi(2^32-1, [n_n, 1]);
        parfor noise=1:n_n
            BNI(:,noise,it)=theta_model(net,T,w,I_0,I_sig,flag,seeds(noise));
        end
        w_save(it)=w;
        bni_aux1=squeeze(mean(mean(BNI(:,:,it)),2));
        bni_aux2=squeeze(mean(mean(BNI(:,:,1:it)),2));
        disp(['Iteration: ' num2str(it) ' | BNI = ' num2str(bni_aux1) ' | w = ' num2str(w)]);
        if it==1
            % Lucky guess:        
            if bni_aux1-BNI_ref<crit && bni_aux1-BNI_ref>0
                x1=1;
            end
            if BNI_ref-bni_aux1<crit && BNI_ref-bni_aux1>0 
                x2=1; 
            end       

            if bni_aux1<BNI_ref
                w=w*2;
            else
                w=w/2;
            end

        else
            % 1st: Find a point above and below BNI_ref
            L1=sum(bni_aux2>BNI_ref);
            L2=sum(bni_aux2<BNI_ref);        
            if L1*L2==0
                if bni_aux1<BNI_ref
                    w=w*2;
                else
                    w=w/2;
                end
                
                if it==n_max
                    z=0;
                end
                continue;
            end
            % 2nd: Fine tuning
            if bni_aux1-BNI_ref<crit && bni_aux1-BNI_ref>0
                x1=1;            
            end
            if BNI_ref-bni_aux1<crit && BNI_ref-bni_aux1>0
                x2=1;            
            end        
            [bni_aux3,index]=sort(bni_aux2);       
            ind1=find(bni_aux3<BNI_ref,1,'last');
            ind2=find(bni_aux3>BNI_ref,1);

            slope=(bni_aux3(ind2)-bni_aux3(ind1))/(w_save(index(ind2))-w_save(index(ind1)));
            yy0=bni_aux3(ind2)-slope*w_save(index(ind2));
            if x1==1
                w=(BNI_ref-displac-yy0)/slope;
            elseif x2==1
                w=(BNI_ref+displac-yy0)/slope;
            else
                w=(BNI_ref-yy0)/slope;
            end
            %w=(w_save(index(ind1))+w_save(index(ind2)))/2; 

        end    
        if (x1+x2==2 || it==n_max)
            z=0;
        end
    end
    w_save(it+1:end)=[];
    BNI(:,:,it+1:end)=[];


    aux = [mean(squeeze(mean(BNI,1)),1)',w_save];
    % sort min to max network BNI
    [~,idx] = sort(aux(:,1));
    % set [BNI, w] values according to sorted BNI (for interpolation below)
    aux = aux(idx,:);



    % the final w value 
    try
        % interpolate, BNI_single.m can make sure this value is good
        w = interp1(aux(:,1),aux(:,2),0.5);
        disp(['result= ', num2str(w)]);
    catch
        disp('nan');
        w = nan;
    end

end
