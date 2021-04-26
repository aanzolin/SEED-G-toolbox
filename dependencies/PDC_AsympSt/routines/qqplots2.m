function [h]= qqplots2(res, resa, pairs, connectivity, measure, kFreq,strFreq)
%     Quantile-quantile plots with PDC/DTF measure for two pairs
%     at lambda = 0.2 rad/s versus Gaussian distribution for significant 
%     connection and Chi2 distribution for H0 (null hypothesis)

if nargin < 7
   strFreq = '.1';
end;

if nargin < 6
   f = 6; % Get third frequency, e.g. f = 0.2 as freq = [0 .1 .2 .3 .4].
else
   f = kFreq;
end;

if nargin < 3,
   pairs = [[2,1]; [2,3]];
   connectivity = [1 0]; % 0 = not connected; 1 = connected.
   measure = 'iDTF';
end;

h=figure;
set(h,'NumberTitle','off','MenuBar','none', ...
   'Name', 'Q-Q plots')
for j=1:2, 
   subplot(1,2,j);
   title([int2str(pairs(j,2)) '-->' int2str(pairs(j,1))]);
   hold on
   r1  = res(:,pairs(j,1),pairs(j,2),f);
   ra1 = resa(:,pairs(j,1),pairs(j,2),f);
   if connectivity(j),
      qqnorm(r1,ra1,measure,strFreq);
      axis square
      disp('qqnorm')
   else
      qqchi2(r1,ra1,measure,strFreq);
      axis square
      disp('qqChi2')
   end;
end;

%==========================================================================
function []= qqchi2(res, resa, measure,strFreq)
%     Quantile-quantile plot for weighted Chi2 distribution
%     res contains values estimated from sample, and
%     resa asymptotic values (mean, var1, patden, patdf)
%
[m, du]=size(res);
[res_s,inx] = sort(res);
x = linspace(1.0/m, 1-1.0/m, m);

%y = icdf('chi2',x, resa(4));
y = chi2inv(x, resa(4));

perc_inf = 0.25;
perc_sup = 0.75;

tick_label = {'.001', '.500', '.750', '.950', '.990', '.999'};
tick_value = [0.001 0.5 0.75 0.95 0.99 0.999];
ytick = y(round(tick_value*m));

xmi = res_s(round(m/100));   % 1% from the lower end.
xma = res_s(round(m - 4));   % Five samples from upper end.

xp25 = res_s(round(perc_inf*m));
xp75 = res_s(round(perc_sup*m));

% yp25=icdf('chi2',perc_inf,resa(4));
% yp75=icdf('chi2',perc_sup,resa(4));

yp25 = chi2inv(perc_inf,resa(4));
yp75 = chi2inv(perc_sup,resa(4));

% Passing a line thru both coordinates
a = (yp75 - yp25)/(xp75 - xp25);
b = yp25 - a*xp25;

ymi = a*xmi + b;
if ymi < 0,
   xmi = -b/a;
   ymi = 0;
end;

yma = a*xma + yp25 - a*xp25;

h1=plot(res_s(2:end-1), y(2:end-1), 'k+','MarkerSize',[8], ...
                                         'MarkerEdgeColor',[.6 .6 .6],...
   'MarkerFaceColor',[.6 .6 .6],'LineWidth',[1.5]);
hold on

if is_octave,
   xlabel(['{|' measure '(\lambda {=' strFreq ')}|}^{2}'])
else
   xlabel(['{|' measure '(\lambda {=' strFreq ')}|}^{2}'], ...
      'Interpreter', 'TeX', 'fontsize', 14)
end;
ylabel('Weighted Chi-square quantile')

h2 = plot([xmi,xma],[ymi,yma], 'k-','linewidth',[2]);

plot([xp25 xp75], [yp25 yp75], 'go')

set(gca,'YTick', ytick, 'YTickLabel',tick_label);
axis square

grid
%
if min(res_s(3))/max(res_s(m-2)) > 0.1
   axis([min(res_s(3)) max(res_s(m-2)) y(1) y(end)]);
else
   axis([0 max(res) y(1) y(end)]);
end;

%==========================================================================
function []=qqnorm(res, resa,measure,strFreq)
%     Quantile-quntile plot for normal distribution
%     res contains sample values 
%     resa asymp values (mean, var1, patden,patdf)
%

m = size(res,1);
[res_s,inx] = sort(res);

x = linspace(1.0/m, 1-1.0/m, m);

% resa(1) % mean
% resa(2) % standard-deviation



y = norminv(x,resa(1),resa(2));  % y = icdf('norm',x,resa(1),resa(2));


perc_inf = 0.25; % usually 25%
perc_sup = 0.75; % and 75%

tick_label = {'.001','.010','.050','.250','.500','.750','.950','.990','.999'};
tick_value = [0.001 0.01 0.05 0.25 0.5 0.75 0.95 0.99 0.999];
ytick = y(round(tick_value*m));

xmi = min(res);
xma = max(res);
%
xp25 = res_s(round(perc_inf*m));
xp75 = res_s(round(perc_sup*m));

% yp25=icdf('norm', perc_inf, resa(1), resa(2));
% yp75=icdf('norm', perc_sup, resa(1), resa(2));

yp25 = norminv(perc_inf, resa(1), resa(2));
yp75 = norminv(perc_sup, resa(1), resa(2));

a = (yp75 - yp25)/(xp75 - xp25);

ymi = a*xmi + yp25 - a*xp25;
yma = a*xma + yp25 - a*xp25;

h1 = plot(res_s(2:end-1), y(2:end-1), 'k+','MarkerSize',[8], ...
              'MarkerEdgeColor',[.6 .6 .6], 'MarkerFaceColor',[.6 .6 .6], ...
              'LineWidth',[1.5]);

hold on

if is_octave,
   xlabel(['{|' measure '(\lambda {=' strFreq ')}|}^{2}'])
else
   xlabel(['{|' measure '(\lambda {=' strFreq ')}|}^{2}'], ...
      'Interpreter', 'TeX', 'fontsize', 14)
end;

ylabel('Normal quantile')
h2 = plot([xmi,xma],[ymi,yma], 'k-','linewidth',[2]);

plot([xp25 xp75], [yp25 yp75], 'go')

set(gca,'YTick',ytick, 'YTickLabel',tick_label);
grid

if min(res)/max(res) > 0.1
   axis([min(res) max(res) y(1) y(end)]);
else
   axis([0 max(res) y(1) y(end)]);
end;
