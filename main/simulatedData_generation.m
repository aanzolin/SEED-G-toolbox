%% Function "simulatedData_generation" 
%  for the generation of simulated datasets (EEG-like), fitting a known connectivity pattern (ground-truth network)
%  
%  Created on  March 2 2019
%  Modified on November 4 2020
%       - Add inputs "val_Range", "ModelDel","AR_choice"
%       - Add output "flag_out", "ModelDel"
%       - Add parameter "max_att"
%       
%  Modified on April 8 2021
%% @author: alessandra anzolin (aanzolin@mgh.harvard.edu)
%
%% Inputs:
%       - DataLength:  number of samples of the generated time series
%       - Trials : number of trials composing the simulated dataset
%       - ModelDel: model distributed on the lags including the real sources
%               on the diagonal 
%               IF the model is not avalable, the function will generate
%               one according to the following inputs.            
%               In this case, impose ModelDel=[]; 
%       - SNR: Signal to Noise Ratio (default = 3)
%       - real_samp: real EEG/cortical data organized as samples x channels/subjects x trials
%       - sig_num: connectivity model size, which correspond to the number of genrated time series 
%       - density: percentage of non-null connections. 
%               Range [0.01 1]. We reccomend choosing density<=0.3
%       - val_Range: range of possible connections values specified as 
%               [-a a], meaning possible value:|a|<1
%       - AR_perc: number of real sources included in the model 
%               expressed as percentage of signals to generate
%               Range [0.01 1]. 
%       - AR_choice: 
%               0 --> the final number of real sources is actually "AR_perc"
%               1 --> after the model generation, an extra real sorce is assigned 
%               to every isolated node (the reason fo this choice is related 
%               with the fact that isolated nodes, very common in big models, do show
%               the spectral properties of white noise)
%       - popt: MVAR model order
%
%% Outputs
%
%       - Y_gen: simulated generated data 
%                samples x signals x trials
%       - E_gen: noise used for the generation
%       - Model Del: model distributed on the lags (including AR real components on the diagonal)
%       - flag_out:
%               1: data-set successfully generated
%               0: impossible to generate the data-set with selected combination of parameters

function [Y_gen, E_gen, ModelDel, flag_out]=simulatedData_generation(DataLength,Trials,ModelDel,SNR,real_samp,sig_num,density,val_Range,AR_perc,AR_choice,popt)

if nargin<4
    SNR = 10;
end

%%%Input for MVGC-Toolbox function 'var_to_tsdata' 
mtrunc      = 0;       %default
decayfac    = 100;     %default
Singtr      = 1;

%%% Maximum number of attempts to generate the dataset given a connectivity model
max_att     = Trials*100; 

%%% Threshold for the signal amplitude
SigLim = 80;      

%%% Connectivity model generation if not provided by the user
if isempty(ModelDel)
    if density>0.3
        disp('Warning: The simulated dataset might not be sucessfully generated for network density higher that 30%. Check output flag_out.')
    elseif sig_num>30 && density>0.1
        disp('Warning: The simulated dataset might not be sucessfully generated for the selected combination of density and number of time series. Check output flag_out.')
    end

    Sw =            eye(sig_num);                 % residual covariance matrix - White Noise 
    AR_num =        round(sig_num*AR_perc);       % number of Auto-Regressive components
    MinDelta =      0.1;                          % minimum difference between two different connections
    Del_Range =     1:popt;
    K =             (sig_num*(sig_num-1))*density;    % number of non-null connections in the imposed model
    
    %%% Inputs for function "arfit_v3"
    selector=      'AIK';
    no_const=      'zero';

    %%% AR parameters evaluation from real data
    if AR_choice
        nAR = sig_num;
    else
        nAR = AR_num;
    end
    
    for r=1:nAR
        for tt=1:size(real_samp,3)
            indperm = randperm(size(real_samp,2));
            EEG_ch = real_samp(:,indperm(r),tt);
            orgEEG{tt}=EEG_ch;
        end %cycle on real trial
        clear tt
        
        ar_coef(:,r) = arfit_v3(orgEEG, popt, popt, selector, no_const);
    end
    clear r
else
    Nod = size(ModelDel,1);
    Sw = eye(Nod); % residual covariance matrix - White Noise 
end %if isempty(ModelDel)

flag=1;
cont=0;

while flag==1
    cont=cont+1;
    flag=0;

    if isempty(ModelDel) 
        %%% Model Structure generation (ch x ch)
        if sig_num == 2
            A = makerandCIJ_dir(sig_num,K);
            inDiag = 1;
        else
            [A, inDiag]=makerandCIJ_dir_realARcomponents(sig_num,K,AR_num);
        end
        
        %%% Model Generation (ch x ch)
        [Model, DelayMatrix]=get_ConnectivityModel(A,val_Range, MinDelta, Del_Range);
        %%% Lag separation (ch x ch x lag)
        ModelDel=rearrangeModel(Del_Range,Model,DelayMatrix);
        
        switch AR_choice
            case {0}
                ARpos = inDiag;
            case {1}
                % Check for isolated nodes
                Nod = size(Model,1);
                Dbin = distance_bin(A);
                Dbin(find(eye(Nod)))=Inf;
                
                for j=1:Nod
                    if length(unique(Dbin(j,:)))==1
                        indInf(j) = 1;
                    else
                        indInf(j) = 0;
                    end
                end 
                
                ARpos = find(indInf); 
                AR_num = length(ARpos);
        end %swich AR_choice

        if AR_num > size(real_samp,2)
            Error('More real data needed in input.')
        end
        
        %%% Add AR coefficients on the main diagonal
        for ii=1:AR_num
            ModelDel(ARpos(ii),ARpos(ii),:)=ar_coef(1:popt,ii)';
        end
        clear ii
    end %if isempty(ModelDel)

    Tr=1;
    num_iter=0;

    while Tr~=Trials+1
        num_iter=num_iter+1;
        
        if num_iter>max_att %Max number of attempts
            Tr=1;
            flag=1;
            num_iter=0;
            break
        end
        
        [Y,E,mtrunc]=var_to_tsdata(ModelDel,Sw,DataLength,Singtr,mtrunc,decayfac);
        
        Ctr=find(abs(Y)>SigLim);
        
        if ~isempty(Ctr)
            continue
        else
            for ch=1:Nod
                nc=randn(1,size(Y,2));
                Ynorm=norm(Y(ch,:));
                ncnorm=norm(nc);
                noise=(Ynorm/(ncnorm*sqrt(SNR)))*nc;
                Ynoise(ch,:)=Y(ch,:)+noise;
            end
            
            Y_gen(:,:,Tr)=Ynoise';
            E_gen(:,:,Tr)=E';
            Tr=Tr+1;
        end
    end 
    
    clear Tr Ctr num_iter
    flag_out=1;
    
    if cont>10
        flag_out=0;
        Y_gen=[];
        E_gen=[];
        break
    end
end    
