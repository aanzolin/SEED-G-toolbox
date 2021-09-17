%% Function "simulatedData_generation_time_varying_connectivity"
%  for the generation of simulated datasets (EEG-like), fitting a time-varying connectivity pattern (ground-truth network)
%
%  Created on  November 5 2019
%  Modified on April 8 2021
%% @authors: alessandra anzolin (aanzolin@mgh.harvard.edu)
%%           jlenia toppi
%% Inputs:
% DataLength :  number of samples
% Trials :      number of trials
% var_perc :    percentage of possible connections to be varied in time
% var_val :     amplitude variation (ES. var_val = 2 means that the final value of the connection will be will be
%               double or half of its initial value)
% ModelDel:     model distributed on the lags + real sources --> if you
%               want to generate the model into the function assign
%               "ModelDel=[]" and insert the following parameters
% SNR :         Signal to Noise Ratio. (default value: 3)
% EEGsamp :     real EEG data (samples x signals x trials)
% sig_num :     model size -> number of signals to be generated
% den :         connections density selected for the model (0.01 0.1 0.3...)
% val_Range:    range of possible connections values
%               format:[-a a] --> possible value:|a|<1
% AR_perc :     number of real sources expressed in term of PERCENTAGE of signals be to generate (0.1 0.2...)
% ARpos
% popt :        optimum MVAR model order
% trans_samp :  number of sample involved in the transition phase
%               0 = 'step' variation
%               1 = 'ramp' variation
%% Outputs:
% Y_gen :       generated data (samples x signals x trials)
% E_gen :       "noise" used for the generation
% Model Del:    model distributed on the lags + real sources
% Model :       basic connectivity model
% var_pos :     indices of the time-varying connections
% flag_out:     "1": data-set successfully generated
%               "0": impossible to generate the data-set with selected parameters

function [Y_gen, E_gen, ModelDel, Model, flag_out]=simulatedData_generation_time_varying_connectivity(DataLength,Trials,Model,SNR,EEGsamp,AR_perc,ARpos,DelayMatrix,popt,trans_samp)

%%% Time-Varying Implementation Parameters
CampIntPerc=[0.5 0.5];
CampInt=CampIntPerc*DataLength;
NumIntCost=length(CampInt);
Del_Range = 1:popt;

%%%Input for MVGC-Toolbox function 'var_to_tsdata'
mtrunc=     0;       %default
decayfac=   100;     %default
Singtr=     1;

%%% Threshold for the signal amplitude
SigLim = 80;

if nargin<4
    SNR = 10;
end

sig_num =       size(Model,1);
Sw =            eye(sig_num);                  % residual covariance matrix - White Noise (input Barnett Toolbox)
AR_num =        round(sig_num*AR_perc);        % number of Auto-Regressive components
max_att =       Trials*100;                    % maximum number of attempts to generate the chosen num. of trials with the same model

%%% Inputs for function "arfit_v3"
selector=      'AIK';
no_const=      'zero';

%%% AR parameters evaluation from real data
for r=1:AR_num
    for tt=1:size(EEGsamp,3)
        EEG_ch = EEGsamp(:,r,tt);
        orgEEG{tt}=EEG_ch;
    end %cycle on real trial
    clear tt
    ar_coef(:,r) = arfit_v3(orgEEG, popt, popt, selector, no_const);
end

flag=1;
cont=0;

while flag==1
    cont=cont+1;
    flag=0;
    
    %%% 1.TV model generation
    %   1.1 Model generation in the constant intervals
    for i=1:NumIntCost
        startInt(i)=sum(CampInt(1:i-1))+sum(trans_samp(1:i-1))+1;
        endInt(i)=startInt(i)-1+CampInt(i);
        aux=reshape(repmat(Model(:,:,i),1,CampInt(i)),[sig_num sig_num CampInt(i)]);
        TVmod(:,:,startInt(i):endInt(i))=aux;
    end
    clear i
    %   1.2 Model generation in the transient
    if trans_samp~=0
        for i=1:NumIntCost-1
            slope=(TVmod(:,:,startInt(i+1))-TVmod(:,:,endInt(i)))/(trans_samp(i)+1);
            for j=1:trans_samp
                TVmod(:,:,endInt(i)+j)=TVmod(:,:,endInt(i))+j*slope;
            end
        end
    end
    clear i
    
    for mo=1:size(TVmod,3)
        ModelDel(:,:,:,mo)=rearrangeModel(Del_Range,TVmod(:,:,mo),DelayMatrix);
        
        %%% Add Auto-Regressive coefficients on the main diagonal
        for ii=1:AR_num
            ModelDel(ARpos(ii),ARpos(ii),:,mo)=ar_coef(1:popt,ii)';
        end
        
        clear ii
    end %cycle on the number of models
    
    clear mo
    
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
        
        [Y,E]=var_to_tsdata_complete(ModelDel,Sw,DataLength,Singtr,mtrunc,decayfac);
        Ctr=find(abs(Y)>SigLim);
        
        if ~isempty(Ctr)
            continue
        else
            for ch=1:sig_num
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
