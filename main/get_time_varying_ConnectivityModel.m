function [Model, ARpos, DelayMatrix, var_pos] = get_time_varying_ConnectivityModel(real_samp,sig_num,density,val_Range,popt,var_val,var_perc,AR_perc,AR_choice)

%% Input:
%       - real_samp: real EEG/cortical data organized as samples x channels/subjects x trials     
%       - sig_num: connectivity model size, which correspond to the number of genrated time series 
%       - density: percentage of non-null connections. 
%               Range [0.01 1]. We reccomend choosing density<=0.3
%       - val_Range: range of possible connections values specified as 
%               [-a a], meaning possible value:|a|<1
%       - popt: MVAR model order
%       - var_val: factor for amplitude variation
%               2 means double the original value of the connection
%       - var_perc: percentage of connections to modify
%               0.5 means that half of the connections are variable over
%               time
%       - AR_perc: number of real sources included in the model 
%               expressed as percentage of signals to generate
%               Range [0.01 1]. 
%       - AR_choice: 
%               0 --> the final number of real sources is actually "AR_perc"
%               1 --> after the model generation, an extra real sorce is assigned 
%               to every isolated node (the reason fo this choice is related 
%               with the fact that isolated nodes, very common in big models, do show
%               the spectral properties of white noise)

%% Output:
%       - Model:          weighted connectivity model.
%                         NxN matrix where the strength of non-null
%                         connections is randomly assigned from the range
%                         'ValuesRange'
%       - ARpos:          position of the AR coefficeints on the main diagonal
%       - DelayMatrix:    NxN matrix in which the entry ij (i,j=1,...N) 
%                         corresponds to the lag associated with the connection
%                         ij. Delay values are randomly selected within the range
%                         DelayRange.
%       - var_pos:        position of the time-varying connections
%
% Created on July 20 2018
%% @author: alessandra anzolin


%%% Fixed parameters for time-varying model generation
CampIntPerc =   [var_perc var_perc];
NumIntCost =    length(CampIntPerc);
Sw =            eye(sig_num);                 % residual covariance matrix - White Noise (input Barnett Toolbox)
AR_num =        round(sig_num*AR_perc);       % number of Auto-Regressive components
MinDelta =      0.1;                          % minimum difference between two different connections
Del_Range =     [1:popt];
K =             (sig_num*(sig_num-1))*density;    % number of connections in the model

% extracting AR parameters from real EEG data to be used in the model for
% reproducing in the simulated signals the same spectral properties of real
% EEG

% input for function "arfit_v2"
selector=      'AIK';
no_const=      'zero';

for r=1:AR_num
    
    for tt=1:size(real_samp,3)
        EEG_ch = real_samp(:,r,tt);
        orgEEG{tt}=EEG_ch;
    end %cycle on real trial
    clear tt
    
    ar_coef(:,r) = arfit_v3(orgEEG, popt, popt, selector, no_const);
end
clear r

% Time-varying Model generation (ch x ch)
[A, inDiag]=makerandCIJ_dir_realARcomponents(sig_num,K,AR_num);
[Model, DelayMatrix]=get_ConnectivityModel(A,val_Range, MinDelta, Del_Range);
%Model = repmat(Model,[1 1 length(var_val)]);

for i=2:NumIntCost
    [Model(:,:,i), var_pos]=applyTransitionTVModel(Model,var_perc,var_val);
end

% to add real sources
switch AR_choice
    case {0}
        ARpos = inDiag;
    case {1} % no isolated nodes
        Dbin = distance_bin(A);
        Dbin(find(eye(sig_num)))=Inf;
        for j=1:sig_num
            if length(unique(Dbin(j,:)))==1
                indInf(j) = 1;
            else
                indInf(j) = 0;
            end
        end %Cycle on Dbin rows
        ARpos=find(indInf); 
end %switch AR_choice