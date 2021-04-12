%% "Script for running the funtion "simulatedData_generation_intertrial_var"
% Created on  Apr 27 2019
% Modified on Apr  8 2021
% 
%% @author: alessandra anzolin (aanzolin@mgh.harvard.edu)

clear all; clc;
%1 Load real data for AR coefficients
%       The real sources provided in the structure MixSub_EEG come from
%       resting state EEG recorde on healthy subjects. Each channel was
%       randomly extracted from one subjects in the dataset 
%
%       Instruction for the user:
%       - download EEG_real_sources.mat (or sLOR_cortical_sources) 
%         from the github page 
%       - specify the directory
%       - if a different read dataset is preferred, be sure to organize it
%       as 'samples x channels x trials' (trials can be equal to 1)
datadir     = 'C:\Users\aless\Dropbox\work\SEED-G_002\GITHUB\req_structures';
name        = 'EEG_real_sources.mat'; %200samples x 45channels x 40realizations
% name      = 'sLOR_cortical_sources.mat'; %125samples x 65dipoles x 130realizations
Struct      = loadname(fullfile(datadir,name));

%2 Set parameters
DataLength =    250;
Trials =        30;
ModelDel =      [];
SNR =           10;
EEGsamp =       Struct.samp;
sig_num =       6;
den =           0.3;
val_Range =     [-0.5 0.5];
AR_perc =       0.3;
AR_choice =     1;
popt =          10;

%% 3 Set parameters for inter-trial variability simulation
%       - SimType: The SEED-G toolbox models two kinds of non-idealities across trials: 
%               'Modified Connections':  variability of the intensity of existing connections hypothesizing that connections values are not consistent across repetitions  
%               'Spurious Links': variability in network density in order to account for the presence of spurious connections in some trials. 
%       - percentage of modified trials
SimType = 'Modified Connections'; 
perc_modTrial = 0.3;

if isequal(SimType,'Modified Connections')
    perc_adSp = 0;
    perc_modCon = 0.5;
    perc_modValue = 0.5;
    mod_Dir = 1;

elseif isequal(SimType,'Spurious Effect')
    perc_adSp = 0.3;
    perc_modCon = 0;
    perc_modValue = 0;
    mod_Dir = 0;
end

[Y, E, Model, Delay_mat, flag_out] = simulatedData_generation_intertrial_var(DataLength,...
    Trials,ModelDel,SNR,EEGsamp,sig_num,den,val_Range,AR_perc,AR_choice,popt,...
    perc_modTrial, perc_adSp, perc_modCon,perc_modValue, mod_Dir);
