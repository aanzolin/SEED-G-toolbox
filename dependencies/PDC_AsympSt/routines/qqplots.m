function []= qqplots(res, resa,pairs,connectivity,measure)
%     Plot figure with quantiles for PDC measure for two pairs
%     at lambda = 0.2 rad/s.
%     Gaussian for significant connection, Chi2 for H0 (null hypothesis)

f = 3; % Get third frequency, i.e. f = 0.2 as freq = [0 .1 .2 .3 .4].
% See asymp_pdc.m for more detail.

if nargin < 3,
   pairs = [[1,2]; [2,1]];
   connectivity = [0 1]; % 0: not connected; 1: connected.
   measure = 'DTF';
end;

h=figure;
set(h,'NumberTitle','off','MenuBar','none', ...
      'Name', 'AsympDTF Simple Examples 1 - Quantile-quantile plots')
for j=1:2,
   subplot(1,2,j);
   r1 = res(:,pairs(j,1),pairs(j,2),f);
   ra1= resa(:,pairs(j,1),pairs(j,2),f);
   if connectivity(j),
      qqnorm(r1,ra1,measure);

      axis square
      disp('qqnorm')
   else
      qqchi2(r1,ra1,measure);

      axis square
      disp('qqchi2')
   end;
end;
%==========================================================================
function []= qqchi2(res, resa, measure)
%     quantile plot for weighted Chi2 distribution
%     res contains values estimated from sample, and
%     resa asymp values (mean,var1, patden, patdf)
%
[m, du]=size(res);
[res_s,inx] = sort(res);
x = linspace(1.0/m, 1-1.0/m, m);

y = icdf('chi2',x, resa(4));

tick_label = {'.001', '.500', '.750', '.950', '.990', '.999'};
tick_value = [0.001 0.5 0.75 0.95 0.99 0.999];
ytick = y(round(tick_value*m));

xmi = min(res);
xma = max(res);

xp25 = res_s(round(0.25*m));
xp75 = res_s(round(0.75*m));
yp25=icdf('chi2',0.25,resa(4));
yp75=icdf('chi2',0.75,resa(4));

a = (yp75 - yp25)/(xp75 - xp25);

ymi = a*xmi + yp25 - a*xp25;
yma = a*xma + yp25 - a*xp25;

h1=plot(res_s(2:end-1), y(2:end-1), 'k+','MarkerSize',[8], 'MarkerEdgeColor',[.6 .6 .6],...
   'MarkerFaceColor',[.6 .6 .6],'LineWidth',[1.5]);
hold on
xlabel(['$\left|' measure '(\lambda {=.1)} \right| ^{2}$'],'Interpreter', 'Latex', 'fontsize', 14)
ylabel('Weighted Chi-square Quantile')
%h2 = plot([xmi,xma],[ymi,yma], 'r:','linewidth',[1.5]);
h2 = plot([xmi,xma],[ymi,yma], 'k-','linewidth',[2]);

set(gca,'YTick',ytick, 'YTickLabel',tick_label);
axis square
v = axis;

grid

if min(res)/max(res) > 0.1
   axis([min(res) max(res) y(1) y(end)]);
else
   axis([0 max(res) y(1) y(end)]);
end;

%==========================================================================
function []=qqnorm(res, resa,measure)
%     qqplot for normal distribution

%     res contains sample values, resa asymp values (mean, var1, patden,
%     patdf)
%
m = size(res,1);
[res_s,inx] = sort(res);

x = linspace(1.0/m, 1-1.0/m, m);

% resa(1) % mean
% resa(2) % standard-deviation

y = icdf('norm',x,resa(1),resa(2));

tick_label = {'.001', '.050', '.250', '.750', '.950', '.999'};
tick_value = [0.001 0.05 0.25 0.75 0.95 0.999];
ytick = y(round(tick_value*m));

xmi = min(res);
xma = max(res);
%
xp25 = res_s(round(0.25*m));
xp75 = res_s(round(0.75*m));

yp25=icdf('norm',0.25,resa(1),resa(2));
yp75=icdf('norm',0.75,resa(1),resa(2));
%
a = (yp75 - yp25)/(xp75 - xp25);

ymi = a*xmi + yp25 - a*xp25;
yma = a*xma + yp25 - a*xp25;

h1=plot(res_s(2:end-1), y(2:end-1), 'k+','MarkerSize',[8], 'MarkerEdgeColor',[.6 .6 .6],...
   'MarkerFaceColor',[.6 .6 .6],'LineWidth',[1.5]);

hold on
xlabel(['$\left|' measure '(\lambda {=.1)} \right| ^{2}$'],'Interpreter', 'Latex', 'fontsize', 14)
ylabel('Normal Quantile')
%h2 = plot([xmi,xma],[ymi,yma], 'b--','linewidth',[1.5]);
h2 = plot([xmi,xma],[ymi,yma], 'k-','linewidth',[2]);

set(gca,'YTick',ytick, 'YTickLabel',tick_label);

grid

if min(res)/max(res) > 0.1
   axis([min(res) max(res) y(1) y(end)]);
else
   axis([0 max(res) y(1) y(end)]);
end;
