%% "Script for running the funtion "simulatedData_generation_time_varying_connectivity"
% Created on  Apr 19 2019
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

%2 Set general parameters
Nodes=          10;           % channels number
Conn_Range=     [-0.5 0.5];   % range for connections values
Density=        0.2;          % network density
AR_perc =       0.5;          % percentage of real sources included in the model
AR_choice=      1;
MinDelta=       0.1;          
popt=           16;           
Trials =        100;
DataLength=     516;
SNR =           10;

%3 Set Time-Varying Parameters
var_perc=       0.5;          % percentage of connections to modify
var_val=        2;            % factor for amplitude variation
trans=          100;          % 0='step' 
                              % N>0 -->'ramp' The variation of the
                              % connection value happens in N samples

%4 Ground-truth connectivity model generation
[Model, inDiag, DelayMatrix, modCon] = get_time_varying_ConnectivityModel(Struct.samp,Nodes,Density,Conn_Range,popt,var_val,var_perc,AR_perc,AR_choice);

%5 Dataset generation
[Y_gen, E_gen, ModelDel, Model, flag_out]=simulatedData_generation_time_varying_connectivity(DataLength,Trials,Model,SNR,Struct.samp,AR_perc,inDiag,DelayMatrix,popt,trans);     
disp(flag_out)
