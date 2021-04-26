%% SEED-G toolbox application:
%% Partial Directed Coherence (PDC) estimation on pseudo-EEG signals
%
% Created on Apr 22 2021
%% @author: alessandra anzolin (aanzolin@mgh.harvard.edu)

clear all;close all;clc

%1. Change directory
datadir = 'C:\Users\aless\Dropbox\work\SEED-G toolbox\demo\simulated data';
savedir = datadir;
numIter = 3;

%3. Setting parameters for PDC estimation and asymptotic statistic
popt=10;                         % optimal order
EndFreq=30;
a=0.05;                          % significance level
corrtype='no_corr';              %'FDR', 'Bonf', 'no_corr'

for it=1:numIter %iterations
    disp(it)
    
    %1. Load pseudo-EEG dataset simulated using SEED-G
    name = sprintf('Sim_dataset_%s',num2str(it));
    EEG = loadname(fullfile(datadir,name));
    %%% Multiple comparison correction
    nodes = size(EEG.samp,2);
    N=nodes*(nodes-1)*EndFreq;
    
    if isequal('Bonf',corrtype)
        aSignif = a/N;
    elseif isequal('FDR',corrtype)
        aSignif = a*((N+1)/(2*N));
    elseif isequal('no_corr',corrtype)
        aSignif = a;
    end
    
    %2. Estimate VAR model parameters
    [~,pf,A] = mvar(EEG.samp',popt,1,1);
    
    %3. Compute PDC
    PDC_str = asymp_pdc(EEG.samp',A,pf,EndFreq,'euc',aSignif);
    PDC = PDC_str.pdc;      %pdc values
    PDCpat = PDC_str.th;    %threshold values
    
    %4. Mean across frequencies
    PDC_m=mean(PDC(:,:,1:EndFreq),3);
    PDCpat_m=mean(PDCpat(:,:,1:EndFreq),3);
    PDCfilt = PDC_m;
    PDCfilt(PDC_m<PDCpat_m)=NaN;
    PDCfilt=PDCfilt-triu(tril(NaN(size(PDCfilt))));
    PDCfilt(isnan(PDCfilt))=0;
    
    %5.Subplot 'Model vs Estimates'
    model = sum(EEG.model,3);
    model = model - triu(tril(model));
        %%% Remove AR components
    
    %6. Evaluate false positive and false negative
    model(model~=0) = 1;
    PDCfilt(PDCfilt~=0) = 1;
    ef = PDCfilt-model;
    indpos_f = find(ef>0);
    indneg_f = find(ef<0);
    fpGlob_f(it) = length(indpos_f)/length(find(model==0))*100;
    fnGlob_f(it) = length(indneg_f)/length(find(model==1))*100;
    %%% Print the percentage of false positive and false negative
    fprintf('Iteration %s \n', num2str(round(it)))
    fprintf('False Positive Rate = %s \n',num2str(round(fpGlob_f(it))))
    fprintf('False Negative Rate = %s \n',num2str(round(fnGlob_f(it))))

    %7. Subplot 'Model VS PDC'
    j = figure;
    scrsz = get(0,'ScreenSize');
    set(j, 'Position', [scrsz(1)+50 scrsz(2)+200 scrsz(3)-100 scrsz(4)-400]);
    set(j, 'Visible', 'on');

    %%% Model 
    M = digraph(model);
    subplot(1,2,1)
    h = plot(M,'Layout',"circle");
    h.EdgeColor = [1.00,0.41,0.16];
    h.NodeColor = [0.49,0.18,0.56];
    h.MarkerSize = 12;
    h.NodeFontSize = 20;
    h.LineWidth = 2;
    h.ArrowSize = 10;
    ax = gca;
    title('Model','FontSize',20);
    ax.TitleHorizontalAlignment = 'left';

    %%% Estimated network
    G = digraph(PDCfilt);
    subplot(1,2,2)
    h = plot(G,'Layout',"circle");
    h.EdgeColor = [1.00,0.41,0.16];
    h.NodeColor = [0.49,0.18,0.56];
    h.MarkerSize = 12;
    h.NodeFontSize = 20;
    h.LineWidth = 2;
    h.ArrowSize = 10;
    ax = gca;
    title('PDC','FontSize',20);
    ax.TitleHorizontalAlignment = 'left';

    suptitle(sprintf('Simulated dataset %s',num2str(it)))
end %if calcolo pdc - caso reale
    