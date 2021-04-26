function c=dtf_alg_special(u,nFreqs,metric,alg,criterion,maxIP,alpha)
%Compute directed transfer function measure given by "option" from series j-->i.
%
%function c=dtf_alg_special(u,alg,criterion,nFreqs,metric,maxIP,alpha)
%
% input: x         - data
%        alg       - algorithm (1: Nutall-Strand),(2: mlsm) ,
%                              (3: Vieira Morf),  (4: QR artfit)
%        criterion - AR order selection criteria =>
%                                   1: AIC; 2: Hanna-Quinn; 3: Schwarz;
%                                   4: FPE, 5: fixed order in MaxIP
%        nFreqs - number of point in [0,fs/2] frequency scale
%        metric   euc  - Euclidean ==> DTF
%                 diag - diagonal ==> DC
%                 info - information ==> iDTF
%        maxIP - externally defined maximum IP
%        alpha - Significance level for null hypothesis testing
%                alpha = .05 is default;
%                if alpha = zero, no asymptotic statistics is computed.
%
% output:  c.A         - AR coefficient matrix estimate
%          c.pf        - covariance matrix estimate
%          c.p         - VAR model order
%          c.metric    - metric for DTF calculation
%          c.nfreqs    - number of point in [0,fs/2] frequency scale
%          c.dtf       - |DTF(f)|^2 estimate
%          c.th        - Threshold DTF(f) value with (1-avalue) significance level.
%          c.dtf_th    - above threshold DTF(f) otherwise DTF(f) = NaN.
%          c.ic1,c.ic2 - confidence interval
%          c.criterion - AR order selection criterion
%          c.alg       - AR estimation algorithm
%          c.alpha     - significance level
%          c.Pass      - 
%          c.Portmanteasu - 
%          c.SS        - Power spectra
%          c.coh       - Coherence function
%          c.chLabels  - Default channel label = [].
%
%% Example of use:
%                 u=sunmeladat([4 3]);  % Andrews & Herzberg 1936-1972
%                                       % sunspot-melanoma series
%                 u=detrend(u);         % Detrend the series 
%                 c=dtf_alg_special(u,64,'diag',1,1,30,0.01); 
%                 figure; xplot(c); %pretty plot     

% If the number of input parameters is smaller than seven, following default
% values are assumed for DTF calculation.
if nargin < 7, alpha = 0; end;      % do not calculate asymptotic statistics
if nargin < 6, maxIP = 30; end;     % defaults value 
if nargin < 5, criterion =  1; end; % AIC order choice
if nargin < 4, alg = 1;  end;       % Nutall-Strand is default AR estimator
if nargin < 3, metric = 'diag'; end;% DC estimation
if nargin < 2, nFreqs = 128; end;   % 128 points on freq scale.
 
[m,n]=size(u);
if m > n,
   u=u.';
end;

nSegLength = length(u);

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv]=mvar(u,maxIP,alg,criterion);

%==========================================================================
%    Testing for adequacy of MAR model fitting through Portmanteau test
%==========================================================================
h = 20; % testing lag
VARadequacy_signif = 0.05;
aValueVAR = 1 - VARadequacy_signif;
flgPrintResults = 1;
[Pass,Portmanteau,st,ths]=mvarresidue(ef,nSegLength,IP,aValueVAR,h,...
                                                          flgPrintResults);
%==========================================================================
%            DTF, threshold and confidence interval calculation.
%==========================================================================

% if alpha == 0, no asymptotic statistics is performed. ASYMP_DTF returns
% only the DTF. This option is much faster!!
c=asymp_dtf_special(u,A,pf,nFreqs,metric,alpha);

c.A = A; c.pf = pf; c.nfreqs = nFreqs;
c.criterion = criterion; c.alg = alg;
c.Pass = Pass; c.Portmanteau = Portmanteau; 

% Power spectra and coherence calculation
c.SS = ss_alg(A, pf, nFreqs);
c.coh = coh_alg(c.SS);

% Statistically significant DTF on frequency scale
if alpha ~= 0,
   dtf_temp = ((abs(c.dtf)-c.th) > 0).*c.dtf + ((abs(c.dtf)-c.th) <= 0)*(-1);
   dtf_temp(ind2sub(size(dtf_temp),find(dtf_temp == -1))) = NaN;
   c.dtf_th = dtf_temp;
else
   c.dtf_th = [];
end;

% Experimental:
c.chLabels=[];