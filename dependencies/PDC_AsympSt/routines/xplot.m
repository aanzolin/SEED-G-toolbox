%% Connectivity plot in matrix layout with power spectra in the main
%  diagonal.flg
%
%function [hxlabel,hylabel] = xplot(c,flgPrinting,fs,w_max,chLabels, ...
%                                     flgColor,flgScale,flgMax,flgSignifColor)
%
% input: 
%    c.{SS,Coh,L,Lpatnaik,LTra,L2vinf,L2vsup,metric}- data structure
%    flgPrinting= [1 1 1 1 1 0 1]; %Plot everything, except coherence.
%     blue-line    | | | | | | 7 Spectra (0: w/o; 1: Linear; 2: Log)
%     gray         | | | | | 6 Coherence
%     dashed-blue  | | | | 5 Plot lower confidence limit (***legacy)
%     dashed-blue  | | | 4 Plot confidence interval
%     red          | | 3 Significant PDC in red line (***legacy)
%     dashed-black | 2 Patnaik threshold level in black dashed-line
%     green        1 PDC/DTF in green line or black w/o statistics,
%                                    see flgSignifColor bellow.
%    fs       - sampling frequency 
%    w_max    - frequency scale upper-limit
%    chLabels - channel identification labels
%    flgColor - 0: white background;
%               1: white [.1 1], light-blue [.01 .1], purple [0 .01]
%               2: white [.1 1], light-gray [.01 .1], gray [0 .01]
%    flgScale - 1: [0 1] / {if max(PDC/DTF) > 1}:[0 max(PDC/DTF)]
%               2: [0 {if max(PDC/DTF) > 1}:max(PDC/DTF)]/[0 1]/[0 .1]/[0 .01]
%                                 based on flgMax(PDC/DTF/Threshold/CI/all)
%               3: [0 max(PDC/DTF/Thr/CI/all)]
%                                 based on flgMax(PDC/DTF/Threshold/CI/all)
%    flgMax   - {'PDC'; 'DTF'; 'Thr'; 'CI'; 'TCI'; 'all'} measure used
%                 as upper limit for Scale. PDC/DTF max value; Thr =
%                 Patnaik threshold maximum value; CI = maximum upper 
%                 confidence interval value; TCI = threshold or CI max
%                 value; all = plotted either PDC/DTF, Thr, CI maximum
%                 value; See also flgScale.
%    flgSignifColor -  0: black line
%                      1: black / gray  => significant / not PDC/DTF
%                      2: red   / gray  =>      "
%                      3: red   / green =>      "
%                      4: red   / black =>      "
%                      5: black / green =>      "
%                    line plotting color of PDC/DTF on the frequency scale. 
% output: graphics
%         hxlabel, hylabel = graphic label's handles
%
% Examples: Calculate c using alg_pdc or alg_dtf.
%
%           xplot(c); % Defaults flgPrinting, fs, w_max and flgColor
%
%           xplot(c,[1 1 1 0 0 1 1],1,100,[],0);
%                        % PDC/DTF, threshold, power spectra, coherence
%                        % plotS. fs=1 Hz; default channel label;
%                        % flgColor=0 => no color or rescalling is used.

function [hxlabel,hylabel] = xplot(c,flgPrinting,fs,w_max,chLabels, ...
                                   flgColor,flgScale,flgMax,flgSignifColor)

warning off

flgPrint = 'Print';

if nargin < 7
   flgScale = 1; % 1: [0 1] / {if max(PDC/DTF) > 1}:[0 max(PDC/DTF)]
   %             % 2: [0 {if max(PDC/DTF) > 1}:max(PDC/DTF)]/[0 1]/[0 .1]/[0 .01]
   %             %                  based on flgMax(PDC/DTF/Threshold/CI/all)
   %             % 3: [0 max(PDC/DTF/Thr/CI/all)]
   %             %                  based on flgMax(PDC/DTF/Threshold/CI/all)
end;

if nargin < 8,
   flgMax = 'TCI'; % {'PDC'; 'DTF'; 'Thr'; 'CI'; 'TCI'; 'all'}
end;

if nargin < 6,
   flgColor = 0; % 0: white background;
   %             % 1: white [.1 1], light-blue [.01 .1], purple [0 .01]
   %             % 2: white [.1 1], light-gray [.01 .1], gray [0 .01]
end;

if nargin < 9,
   flgSignifColor = 3; % 0: black line
   %                   % 1: black / gray -> significant /not signif PDC/DTF
   %                   % 2: red  /  gray ->      "
   %                   % 3: red  / green ->      "
   %                   % 4: red  / black ->      "
end;

flgYTickLabel = 0; % TO BE DEFINED.

knargin = 6; % Number of input arguments control <<legacy code - review>>>!

if isfield(c,'dtf')
   L = c.dtf;       % |DTF|^2 (N x N x freq)
   if isfield(c,'dtf_th')
      LTra   = c.dtf_th; % Significant PDC^2 on freq range, otherwise NaN.
   end;
   flgType = 'DTF';
elseif isfield(c,'pdc')
   L = c.pdc;       % |PDC|^2 (N x N x freq)
   if isfield(c,'pdc_th')
      LTra   = c.pdc_th; % Significant |PDC|^2 on freq range, otherwise NaN.
   end;
   flgType = 'PDC';
else
   error('Variable c does not hold PDC/DTF analysis results.')
end

SS = c.SS;       % Spectra
Coh = c.coh;     % |Coh|^2
Lpatnaik = c.th; % Patnaik threshold values for alpha (level of significance)
L2vinf = c.ic1;
L2vsup = c.ic2;
metric = c.metric; % Valid options: "euc", "diag" or "info"
[N,q,nFreqs]=size(L); % N = q = Number of channels/time series;
%                     % nFreqs = number of points on frequency scale

nodesett = 1:N;

if nargin <  (knargin-4),  flgPrinting = [1 1 1 0 0 0 1]; end;
if nargin <= (knargin-4),  fs = 1; end;
if nargin <  (knargin-2),  w_max = fs/2; end;
if nargin <   knargin,     flgColor = 0; end;

if w_max > fs/2 + eps,
   error(['The parameter w_max should be =< Nyquist frequency,' ...
      'i.e, w_max <= fs/2.'])
end;

if c.alpha == 0, % No asymptotic statistics calculation performed.
   flgPrinting = flgPrinting .* [1 0 0 0 0 1 1];
end;

if exist('shadedplot.m','file') ~= 2, flgColor = 0; end;

w = 0:fs/(2*nFreqs):w_max-fs/(2*nFreqs);
nPlotPoints = length(w);
w_min = w(1);

if fs == 1,
   str_w_max = '.5';
elseif fs < 1,
   str_w_max = sprintf('%0.5g', w_max);
elseif fs < 20,
   str_w_max = sprintf('%1.2f', w_max);
   
elseif fs < 200,
   if mod(fs,2),
   str_w_max = sprintf('%3.1f', w_max);
   else
   str_w_max = sprintf('%2.0f', w_max);
   end      
elseif fs <= 10000,
   str_w_max = sprintf('%5.0f', w_max);
else
   str_w_max = 'fs/2';
   
end

if nargin < (knargin-1) || isempty(chLabels),
   if isfield(c,'chLabels'),
      chLabels = c.chLabels;
   else
      chLabels = [];
   end;
   if ~isempty(chLabels) && max(size(chLabels)) < N,
      error('1 NOT ENOUGH CHANNEL LABELS.');
   end;
elseif max(size(chLabels)) < N,
   if isfield(c,'chLabels'),
      chLabels = c.chLabels;
   else
      chLabels=[];
   end;
   if ~isempty(chLabels) && max(size(chLabels)) < N,
      error('2 NOT ENOUGH CHANNEL LABELS 2.');
   else
      disp('3 NOT ENOUGH CHANNEL LABELS. Default labels assumed.');
   end;
end;

hxlabel=0; % x-axis labels' handles
hylabel=0; % y-axis labels' handles

%==========================================================================
%==========================================================================

for j = 1:N,
   s = nodesett(j);
   for i = 1:N,
      r = nodesett(i);
      if j ~= i || ( j == i && flgPrinting(7) ~= 0)
%         h=subplot2(N,N,(i-1)*N+j);
         h=subplot2(N,N,(i-1)*N+j);
      end;
%=========================================================================
%                Main Diagonal Power spectrum plotting
%=========================================================================
      if (j == i) && flgPrinting(7) ~= 3, %Power spectrum
         if flgPrinting(7) ~= 0,

            SStmp = abs(getCij(SS,r,s,nPlotPoints));
            Ltmp  = abs(getCij(L, r,s,nPlotPoints));

            switch flgPrinting(7), % Main diagonal plotting SS and/or PDC

               case 1,
                  %Normalization of linear power scale
                  SStmp = SStmp/max(SStmp);
                  h12   = plot(w,SStmp);
                  ax(1) = gca;
                  ax(2) = ax(1);
                  if j == 2,
                     hh=ylabel('Spectra [a.u.]');
                     pos = get(hh,'Position');
                     pos(1) = 0.005;
                     set(hh,'Position',pos, ...
                        'FontSize',[8]);
                  end;
                  set(h12,'LineWidth',2.5,'Color',[0 0 0.7]);
               case 2,
                  SStmp=log(SStmp);
                  SStmp=(SStmp-min(SStmp))/(max(SStmp)-min(SStmp));
                  h12 = plot(w,SStmp);
                  ax(1)=gca;
                  ax(2)=ax(1);
                  if j == 2,
                     hh=ylabel('log Spectra [a.u.]');
                     pos = get(hh,'Position');
                     pos(1) = 0.005;
                     set(hh,'Position',pos, ...
                        'FontSize',[8]);
                  end;
                  set(h12,'LineWidth',2.5,'Color',[0 0 0.7]);

%                case 3, % PDCii or DTFii plot on main diagoinal. ****
%                   atrib = 'g-';
%                   atrib = [0.41,0.41,0.41];
%                   plot(w,Ltmp,'Color',atrib,'LineWidth',2.5);
%                   hold on
%                   ax(1)=gca;
%                   ax(2)=ax(1);
% 
%                   if flgPrinting(2)==1
%                      atrib='k--'; % Patnaik signif. level in black line
%                      plot(w,abs(getCij(Lpatnaik,r,s,nPlotPoints)), ...
%                         atrib,'LineWidth',1.5);
%                   end;
%                   if flgPrinting(3)
%                      atrib='k-'; % Significant PDCn in black line
%                      plot(w,abs(getCij(LTra,r,s,nPlotPoints)), ...
%                         atrib,'LineWidth',2.5);
%                   end;

               case 4, %Normalization of linear power scale
                  SStmp=SStmp/max(SStmp);
                  [ax,h1,h2]=plotyy(h,w,SStmp,w,Ltmp);
                  if j == 2,
                     hh=ylabel('Spectra [a.u.]');
                     pos = get(hh,'Position');
                     pos(1) = 0.005;
                     set(hh,'Position',pos, ...
                        'FontSize',[8]);
                  end;

                  set(h1,'LineWidth',2.5,'Color',[0 0 0.7]);
                  set(h2,'LineWidth',2.5,'Color',[0.41 0.41 0.41]);
                  hold(ax(1), 'all'); hold(ax(2), 'all');
                  if flgPrinting(3)
                     atrib='r-'; % Significant PDCn in black (or red) line
                     plot(w,abs(getCij(LTra,r,s,nPlotPoints)), ...
                        atrib,'LineWidth',2.5);
                  end;
                  if flgPrinting(2)==1
                     atrib='k--'; % Patnaik signif. level in black line
                     plot(w,abs(getCij(Lpatnaik,r,s,nPlotPoints)), ...
                        atrib,'LineWidth',1.5);
                  end;


               case 5, % log Spectra
                  SStmp = log(SStmp);
                  SStmp = (SStmp-min(SStmp))/(max(SStmp)-min(SStmp));
                  [ax,h1,h2]=plotyy(h,w,SStmp,w,Ltmp);
                  if j == 2,
                     hh=ylabel('log Spectra [a.u.]');
                     pos = get(hh,'Position');
                     pos(1) = 0.005;
                     set(hh,'Position',pos, ...
                        'FontSize',[8]);
                  end;

                  set(h1,'LineWidth',2.5,'Color',[0 0 0.7]);
                  set(h2,'LineWidth',2.5,'Color',[0.41 0.41 0.41]);
                  hold(ax(1), 'all'); hold(ax(2), 'all');
                  if flgPrinting(3)
                     atrib='k-'; % Significant PDCn in black (or red) line
                     plot(w,abs(getCij(LTra,r,s,nPlotPoints)), ...
                        atrib,'LineWidth',2.5);
                  end;
                  if flgPrinting(2)==1
                     atrib='k--'; % Patnaik signif. level in black line
                     plot(w,abs(getCij(Lpatnaik,r,s,nPlotPoints)), ...
                        atrib,'LineWidth',1.5);
                  end;

            end;
            hold on
            grid on

            ylim =[-0.050 1.10];

            set(ax(1),'XLim',[w_min w_max],'XTick',[],'XTickLabel',[' '],...
               'YLim',ylim); %,'YTick',[]);
            if j == N && (flgPrinting(7) == 4 || flgPrinting(7) == 5)
               set(ax(2),'XLim',[w_min w_max],'XTick',[], ...
                  'XTickLabel',[' '],'YLim', ylim); %, ... % [-0.02 1.05], ...
               % 'YTick',[0 .5 1]); %, 'YTickLabel',[' 0'; '.5'; ' 1']);
            else
               set(ax(2),'XLim',[w_min w_max],'XTick',[], ...
                  'XTickLabel',[' '],'YLim', ylim);%, ... % [-0.02 1.05], ...
               %                  'YTick',[0 .5 1]);%, 'YTickLabel',[' ']);
            end;
            set(h, 'Color', 0.8*[1 1 1]); % Background color
            if j == N,
               hxlabel(j)=labelitx(j,chLabels);
               set(h,'XTickLabel',[' ';' ';' ';' ';' ';' '])
            else
               set(h,'XTickLabel',[' ';' ';' ';' ';' ';' '])
            end;
            if j == 1,
               hylabel(i)=labelity(i,chLabels);
            end;
            if j == N,
               hxlabel(j)=labelitx(j,chLabels);
            end;
         end;
         if j == 1,
            hylabel(i)=labelity(i,chLabels);
         end;
         if j == N,
            hxlabel(j)=labelitx(j,chLabels);
         end;


%% ========================================================================
%                      PDC and Coh2 plotting
%% ========================================================================
      else % PDC and coherence
         if flgPrinting(1),
            %atrib='g-'; % PDCn in green lines
            Ltmp          = abs(getCij(L,r,s,nPlotPoints));
            Cohtmp        = abs(getCij(Coh,r,s,nPlotPoints));
            if c.alpha ~= 0,
               LTratmp        = abs(getCij(LTra,r,s,nPlotPoints));
               Lpatmp         = abs(getCij(Lpatnaik,r,s,nPlotPoints));
               L2vsuptmp      = abs(getCij(L2vsup,r,s,nPlotPoints));
               %kSignif        = abs(getCij(LTra,r,s,nPlotPoints));
               flgSignif      = sum(~isnan(LTratmp));
               indexSignif    = ~isnan(LTratmp);
               indexNotSignif = isnan(LTratmp);
               atrib2 = [0.41,0.41,0.41];
            else
               if c.alpha == 0,     % With no statistics, PDC is
                  atrib1 = [0 0 0]; % plotted in black lines.
                  atrib2 = [0 0 0]; 
               end;
            end;
            grayline = 0.71;
            switch flgSignifColor
               case 0
                  atrib1 = [0 0 0];  %black line
                  atrib2 = [0 0 0];  %black line
               case 1
                  atrib1 = [0 0 0];  % 1: black / gray -> significant /not signif PDC/DTF
                  atrib2 = grayline*[1 1 1];  %gray line
               case 2
                  atrib1 = [1 0 0];  % 2: red  /  gray ->      "
                  atrib2 = grayline*[1 1 1];  %gray line
               case 3
                  atrib1 = [1 0 0];  % 3: red  / green ->      "
                  atrib2 = [0 1 0];  %    green ->      "
               case 4
                  atrib1 = [1 0 0];  % 4: red        "
                  atrib2 = [0 0 0];  %    black ->      "
               otherwise
                  atrib1 = [0 0 0];  % 5: black        "
                  atrib2 = [0 1 0];  %    green ->      "
            end;
            h12 = plot(h,w,Ltmp,'Color',atrib2,'LineWidth',[2.5]);            
            ax(1)=gca;
            ax(2)=ax(1);
            grid off
            hold on

            if c.alpha ~= 0,
               maxCI  = max(L2vsuptmp(indexSignif));
               maxPDC = max(Ltmp);
               maxThr = max(Lpatmp);
               maxCoh = max(Cohtmp);
            else
               maxCI  = max(Ltmp);
               maxPDC = maxCI;
               maxThr = maxCI;
               maxCoh = max(Cohtmp);
            end;


            switch flgScale
               case 1 % [0 1] / [0 max(PDC)]
                  if maxPDC <= 1,
                     ylim = [-0.050 1.10];
                     maxValue = 1.0;
                  else
                     ylim = maxPDC*[-0.05 1.1];
                     maxValue = maxPDC;
                  end

               case {0,2,3} % [0 max(PDC)]/[0 1]/[0 .1]/[0 .01]
                  switch upper(flgMax)
                     case 'PDC'
                        maxValue = maxPDC;
                     case 'Thr'
                        if flgPrinting(2),
                           maxValue = max([maxThr maxPDC]);
                        else
                           maxValue = maxPDC;
                        end;
                     case 'CI'
                        if flgSignif && flgPrinting(4),
                           maxValue = max([maxCI maxPDC]);
                        else
                           maxValue = maxPDC;
                        end
                     case 'COH'
                        if flgPrinting(6),
                           maxValue = max([maxCoh maxPDC]);
                        else
                           maxValue = maxPDC;
                        end
                     case 'TCI'
                        if (flgSignif && flgPrinting(4)) && flgPrinting(2),
                           maxValue = max([maxPDC maxThr maxCI]);
                        elseif flgSignif && flgPrinting(4),
                           maxValue = max([maxPDC maxCI]);
                        elseif flgPrinting(2),
                           maxValue = max([maxPDC maxThr]);
                        else
                           maxValue = maxPDC;
                        end

                     case 'ALL'
                        if (flgSignif && flgPrinting(4)) && flgPrinting(2) ...
                                                         && flgPrinting(6),
                           maxValue = max([maxPDC maxThr maxCI maxCoh]);
                        elseif (flgSignif && flgPrinting(4)) && flgPrinting(2),
                           maxValue = max([maxPDC maxThr maxCI]);
                        elseif (flgSignif && flgPrinting(4)) && flgPrinting(6),
                           maxValue = max([maxPDC maxCoh maxCI]);
                        elseif flgPrinting(2) && flgPrinting(6),
                           maxValue = max([maxPDC maxThr maxCoh]);
                        elseif flgPrinting(2),
                           maxValue = max([maxPDC maxThr]);
                        elseif (flgSignif && flgPrinting(4)),
                           maxValue = max([maxPDC maxCI]);
                        elseif flgPrinting(6),
                           maxValue = max([maxPDC maxCoh]);
                        else
                           maxValue = maxPDC;
                        end
                     otherwise
                        maxValue = 1.0;
                  end
            end;
            if flgScale == 2
                  if maxValue < 0.0001,
                     ylim = [-0.00004 0.00210];
                  elseif maxValue < 0.0105,
                     ylim = [-0.0002 0.0105];
                  elseif maxValue < 0.105,
                     ylim = [-0.002 0.105];
                  elseif maxValue < 1.05
                     ylim = [-0.02 1.05];
                  else
                     ylim = maxValue*[-0.02 1.05];
                  end;
            elseif flgScale == 3
                  ylim = maxValue * [-0.02 1.05];
            end;
%% ========================================================================
            switch flgColor
               case 0,
                  % White background.
               case 1,
                  if maxValue < 0.001,
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        [0.7 0.7 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.8*[1 0.7 1]);
                  elseif maxValue < 0.01,
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        [0.7 0.7 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.8*[1 0.7 1]);
                  else
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        [0.7 0.7 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.8*[1 0.7 1]);
                  end;
               case 2,
                  if maxValue < 0.001,
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        0.9*[1 1 1]); %light gray
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.75*[1 1 1]);  % darker gray
                  elseif maxValue < 0.01,
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        0.9*[1 1 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.75*[1 1 1]);
                  else
                     shadedplot([w_min w_max],[0.01 0.01],[0.1 0.1], ...
                        0.9*[1 1 1]);
                     hold on
                     shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                        0.75*[1 1 1]);
                  end;
            end;
         end;
         grid off

%======================================================================
%                         Labeling axis
%======================================================================
         if i == N,
            if j == 1,
               hylabel(i)=labelity(i,chLabels);
               hxlabel(j)=labelitx(j,chLabels);
               if flgColor == 4,
                  if max(Ltmp) < 0.001,
                     set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' 0';'  ';'  ';'  ';'  ';str_w_max], ...
                        'FontSize',10,'FontWeight','bold')
                  else
                     set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                        'XTick',[0 .1 .2 .3 .4 .5], ...
                        'XTickLabel',[' 0';'  ';'  ';'  ';'  ';str_w_max], ...
                        'FontSize',10)
                  end;
               else
                  set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                     'XTick',2*w_max*[0 .1 .2 .3 .4 .5], ...
                     'XTickLabel',{' 0';'  ';'  ';'  ';'  ';str_w_max}, ...
                     'FontSize',10)
               end;
            else
               hxlabel(j)=labelitx(j,chLabels);
               set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                  'XTick',[0 .1 .2 .3 .4 .5], ...
                  'XTickLabel',[' ';' ';' ';' ';' ';' '])
            end;
         elseif i==1 && j == 2 && flgPrinting(7)==0,
            hylabel(i)=labelity(i,chLabels);
            set(h,'XLim', [w_min w_max],'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '])
         elseif j == 1,
            hylabel(i)=labelity(i,chLabels);
            set(h,'XLim', [w_min w_max],'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ... 
               'FontSize',10)
            if i == N,
               set(h,'FontSize',10,'FontWeight','bold');
            end;
         else
            set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'FontSize',10)
         end;

         hold on

%==========================================================================
%==========================================================================


         % Plot lower and upper bound of confidence interval
         if flgPrinting(4),
            if flgPrinting(4) >= 2, % Plot error bar
               if flgPrinting(4)== 2,
                  kstep=1; % Errorbar spacing parameter
               else        % flgPrinting(4) == 3.
                  kstep=length(w)/32; % Limit number of error bar on x-axis
                  %                   % to 32.
               end;
               L2vsuptmp = abs(getCij(L2vsup,r,s,nPlotPoints));
               L2vsuptmp(indexNotSignif) = NaN;
               Ltmp2 = Ltmp; %= = abs(getCij(L,r,s,nPlotPoints)); 
               Ltmp2(indexNotSignif)=NaN;
               varPDC2 = L2vsuptmp-Ltmp2;
               % Sparse error bar plotting for confidence interval 
               web=kstep:kstep:length(w); % for sparse plotting interval.
               b=zeros(length(web),2,1);
               b(:,1,1) = varPDC2(web); 
               b(:,2,1) = varPDC2(web);
               
               Nweb = length(web);
               

               if flgPrinting(4)== 2,
                  indexTemp=ones(1,Nweb);
                  %web(webb) index with NaN values.
                  isnanIndex=find(isnan(varPDC2(web)));
                  indexTemp(isnanIndex) = 0;
                  webb=web(indexTemp>0); % webb contains not NaN web indices.
                  webbProbe = [-1 webb Nweb+2]; % Temporary variable for probing not NaN elements.
                  webbDiff = diff(webbProbe);
                  indexNotNaN = find(webbDiff>1);
                  %webbDiffIndex = webb(webbDiff(1:end-1) > 1);
                  nDiff2 = sum(webbDiff > 1);

                  if nDiff2 ~= 0,
                     for kDiff = 1:nDiff2-1,
                        indexPlotVar0   = webb(indexNotNaN(kDiff));
                        indexPlotVarEnd = webb(indexNotNaN(kDiff+1)-1);

                        webbIndexTemp = web(indexPlotVar0:indexPlotVarEnd);
                        wtmp = w(webbIndexTemp); wtmp=[wtmp(1) wtmp];
                        ytmp = Ltmp(webbIndexTemp); ytmp=[ytmp(1); ytmp];
                        btmp = b(webbIndexTemp,:,:);  btmp=[[0 0]; btmp];
                        if is_octave(),
                            hE = boundedline(wtmp, ytmp,btmp,'-k');
                        else
                            hE = boundedline(wtmp, ytmp,btmp,'-k', 'alpha');
                        end;
                     end;
                  else
                     wtmp = w(webb); wtmp=[wtmp(1) wtmp];
                     ytmp = Ltmp(webb); ytmp=[ytmp(1); ytmp];
                     btmp = b(webb,:,:);  btmp=[[0 0]; btmp];
                     if is_octave()
                        hE = boundedline(wtmp, ytmp,btmp(:,1), '-k');
                     else
                        hE = boundedline(wtmp,ytmp,btmp(:,1),'-k', 'alpha');
                     end;
                  end;
               else
                  hE = errorbar(w(web),Ltmp(web),varPDC2(web));
               end;
               
               if exist('hE','var'),
                  set(hE,'LineWidth',0.75, ...
                     'color',[0 0 0]); %[0.4 0.4 0.4]);
               end;

            else
               if flgPrinting(4) && flgSignif
                  atrib='k--'; % Lower-bound significant level
                  tmp= getCij(L2vinf,r,s,nPlotPoints);
                  % lower-bound can be negative.
                  tmp(indexNotSignif)=NaN;
                  plot(w,tmp, atrib,'LineWidth',1.5,'Color',[.7 .7 .7]);
               end;
               if flgPrinting(4) && flgSignif
                  atrib='k--'; % Upper-bound significant level
                  tmp = abs(getCij(L2vsup,r,s,nPlotPoints));
                  tmp(indexNotSignif)=NaN;
                  plot(w,tmp, atrib,'LineWidth',1.5,'Color',[.7 .7 .7]);
               end;
            end;
         end;
                  
         set(gca, ...
            'Box'         , 'on'     , ...
            'TickDir'     , 'in'     , ...
            'TickLength'  , [.02 .02] , ...
            'XMinorTick'  , 'off'      , ...
            'YMinorTick'  , 'off'      );
         if flgPrinting(6), %|Coh|^2 plot in gray-line
            plot(w,Cohtmp,'-','LineWidth',2.5, ...
               'Color',[.7 .7 .7]);
         end;

         plot(h,w,Ltmp,'Color',atrib2,'LineWidth',2.5);
         if c.alpha ~= 0,
            plot(h,w,LTratmp,'Color',atrib1,'LineWidth',2.5);
         end;
         
         if flgPrinting(2),
            atrib='k--'; % Patnaik significance level in black line
            plot(w,abs(getCij(Lpatnaik,r,s,nPlotPoints)),atrib, ...
                  'LineWidth',1.5);
         end;
      end; % PDC and coherence
   end;
end;
%

supAxes = [.08 .08 .84 .84];

[ax,h1] = suplabel('Frequency','x',supAxes); set(h1, 'FontSize',14)

if strcmp(flgType,'DTF'),
   gType = 'gamma';
   switch metric
      case 'euc'
         pType = 'DTF';
      case 'diag'
         pType = 'DC';
      case 'info'
         pType = '_{i}DTF';
      case 'ratio'
         pType = 'Ratio';
      otherwise
         error('Unknown metric.')
   end;
else
   gType = 'pi';
   switch metric
      case 'euc'
         pType = 'PDC';
      case 'diag'
         pType = '_{g}PDC';
      case 'info'
         pType = '_{i}PDC';
      case 'ratio'
         pType = 'iPDC/gPDC Ratio';
      otherwise
         error('Unknown metric.')
   end;
end;

switch lower(flgPrint)
   case 'screen'
      switch metric
         case 'euc'
            [ax,h1]=suplabel(['{\mid\' gType '_{\it{{i}{j}}}{(\lambda)\mid}^{2}}'],...
               'y',supAxes);
         case 'diag'
            [ax,h1]=suplabel(['{\mid{_{\it{g}}}\' gType '_{\it{{i}{j}}}{(\lambda)\mid}^{2}}'],...
               'y',supAxes);
         case 'info'
            [ax,h1]=suplabel(['{\mid{_i}\' gType '_{\it{{i}{j}}}{(\lambda)\mid}^{2}}'],...
               'y',supAxes);
         case 'ratio'
            [ax,h1]=suplabel(['i' flgType ' to g' flgType ' Ratio'],...
               'y',supAxes);
         otherwise
            error('Unknown metric.')
      end;
   otherwise %Print
      if is_octave(),
         xplot_title(c); % In octave suplabel('x' and 'y') does not work.
      else
         [ax,h1] = suplabel(['{|' pType '(\lambda)|}^{2}'],...
            'y',supAxes);
         set(h1, 'FontSize',14)
         pos = get(h1,'Position');
         pos(1) = pos(1)+0.05; %0.0545
         set(h1,'Position',pos) % Adjust ylabel position
         
      end;
end;

% Adjusting axis labels positions.
for k=1:N,
   set(hxlabel(k),'Units','normalized');
   set(hylabel(k),'Units','normalized');
   pos = get(hylabel(k),'Position');
   pos(1) = -0.135*N/3;
   set(hylabel(k),'Position',pos);
   pos = get(hxlabel(k),'Position');
   pos(2) = -0.145*N/3;
   set(hxlabel(k),'Position',pos);
end;

%==========================================================================
function [hxlabel]=labelitx(j,chLabels) % Labels x-axis plottings
if isempty(chLabels)
   hxlabel = xlabel(['j = ' int2str(j)]);
   set(hxlabel,'FontSize',12); % , ... %'FontWeight','bold', ...
%      'FontName','Arial') % 'FontName','Arial'
else
   hxlabel=xlabel([chLabels{j}]);
   set(hxlabel,'FontSize',12) %'FontWeight','bold')
end;

%% ========================================================================

function [hylabel]=labelity(i,chLabels) % Labels y-axis plottings
if isempty(chLabels)
   hylabel=ylabel(['i = ' int2str(i)],...
      'Rotation',90);
   set(hylabel,'FontSize',12); %, ... %'FontWeight','bold', ...
%      'FontName','Arial')  % 'FontName','Arial', 'Times'
else
   hylabel=ylabel([chLabels{i}]);
   set(hylabel,'FontSize',12); %'FontWeight','bold','Color',[0 0 0])
end;

%==========================================================================
