%% "Script for running the funtion "simulatedData_generation"
% Created on  Apr 19 2019
% Modified on Apr  8 2021
% 
%% @author: alessandra anzolin (aanzolin@mgh.harvard.edu)

%1 Load real data for AR coefficients
%       The real sources provided in the structure EEG_real_sources.mat come from
%       resting state EEG recorded on healthy subjects. Each channel was
%       randomly extracted from one subjects in the dataset 
%
%       Instruction for the user:
%       - download EEG_real_sources.mat (or sLOR_cortical_sources) 
%         from the github page 
%       - specify the directory
%       - if a different read dataset is preferred, be sure to organize it
%       as 'samples x channels x trials' (trials can be equal to 1)

clear all; clc;
datadir     = 'C:\Users\aless\Dropbox\work\SEED-G_002\GITHUB\req_structures';
name        = 'EEG_real_sources.mat'; %200samples x 45channels x 40realizations
% name      = 'sLOR_cortical_sources.mat'; %125samples x 65dipoles x 130realizations
Struct      = loadname(fullfile(datadir,name));

%2 Set parameters
DataLength  = 250;
Trials      = 30;
DataType    = 'scalp'; %'cortex'
ModelDel    = []; 
SNR         = 10;
samp        = Struct.samp;
sig_num     = 10;
density     = 0.3;
val_Range   = [-0.5 0.5];
AR_perc     = 0.3;
AR_choice   = 1;
popt        = 10;

[pseudoEEG, E_gen, ModelDel, flag_out]=simulatedData_generation(DataLength,...
    Trials,ModelDel,SNR,samp,sig_num,density,val_Range,AR_perc,AR_choice,popt);

if flag_out
    disp('simulated dataset succesfully generated!')
end