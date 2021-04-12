%% "Script for running the funtion "simulatedData_generation_with_volume_conduction"
% Created on  Apr 19 2019
% Modified on Apr  8 2021

%1 Load real data for AR coefficients
%       The real sources provided in the structure sLOR_cortical_sources come from
%       resting state EEG recorded on healthy subjects, recontructed in the
%       source domain emploing the algorithm sLORETA[1]. 
%
%       Instruction for the user:
%       - sLOR_cortical_sources (or ECoG signals)
%       - specify the directory
%       - if a different read dataset is preferred, be sure to organize it
%       as 'samples x channels x trials' (trials can be equal to 1)
clear all; clc;
datadir   = 'C:\Users\aless\Dropbox\work\SEED-G toolbox\real data';
name      = 'sLOR_cortical_sources.mat'; %125samples x 65dipoles x 130realizations
Struct    = loadname(fullfile(datadir,name));

%2 Set parameters
DataLength =    250;
Trials =        30;
ModelDel =      [];
SNR =           10;
sig_num =       6;
density =       0.3;
val_Range =     [-0.5 0.5];
AR_perc =       0.3;
AR_choice =     1;
popt =          10;

%3. Specify if the data should be projected on the scapl after the
%        generation
forward = 1;
%       1 for forward problem solution
%       0 otherwise (no prjection in the sensor space)

% If you select forwar = 1 --> please specify channel labels and ROI labels
% 3.1 Specify list of electrodes for forward model
    sensors_labels = {'C3' 'C4' 'CP1' 'CP2' 'CP5' 'CP6' 'Cz' 'F3' 'F4'...
        'F7' 'F8' 'Fz' 'AFz' 'FC1' 'FC2' 'FC5' 'FC6' 'Fp1' 'Fp2' 'Fpz' 'O1'...
        'O2' 'Oz' 'P3' 'P4' 'P7' 'P8' 'POz' 'Pz' 'T7' 'T8'}; % Example with 31 channels
%3.2 ROI selection
    %   Each generated signal will be located on the brain cortex according to
    %   the provided list of ROIs labels
    ROI_labels = {'Brodmann area 19_L' 'Brodmann area 19_R'...
        'Brodmann area 4_L' 'Brodmann area 4_R' 'Brodmann area 5_L' 'Brodmann area 5_R'};
    if length(ROI_labels)~=sig_num
        error('The number of ROI lables must be equal to the number of generated signals.')
    end
%3.3 path for dependencies NYH and Fieldtrip
nyh_path = 'C:\Users\aless\Dropbox\work\SEED-G toolbox\dependencies\NYH'; 
ft_path = 'C:\Users\aless\Dropbox\work\SEED-G toolbox\dependencies\Fieldtrip'; 

[Y_gen, E_gen, ModelDel, flag_out]=simulatedData_generation_with_volume_conduction(DataLength,...
    Trials,ModelDel,SNR,Struct.samp,sig_num,density,val_Range,AR_perc,AR_choice,popt,forward,sensors_labels,ROI_labels,nyh_path,ft_path);
