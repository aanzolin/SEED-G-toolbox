function []= qqplots3(res, resa,pairs,connectivity,measure,kFreq)
%     Plot figure with quantiles for PDC/DTF measure for two pairs
%     at lambda = 0.2 rad/s.
%     Gaussian for significant connection, Chi2 for H0 (null hypothesis)

if nargin < 6
   f = 6; % Get third frequency, e.g. f = 0.2 as freq = [0 .1 .2 .3 .4].
else
   f = kFreq;
end;

if nargin < 3,
   pairs = [[2,1]; [1,2]];
   connectivity = [0 1]; % 0 = not connected; 1 = connected.
   measure = ' ';
end;

h=figure;
set(h,'NumberTitle','off','MenuBar','none', ...
   'Name', 'Quantile-quantile plots (square-root plot for chi2')
for j=1:2,
   subplot(1,2,j);
   title([int2str(pairs(j,2)) '-->' int2str(pairs(j,1))]);
   hold on
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
%     resa asymp values (mean, var1, patden, patdf)
%

[m, du]=size(res);
[res_s,inx] = sort(res);
x = linspace(1.0/m, 1-1.0/m, m);
y = icdf('chi2', x, resa(4));

perc_inf = 0.25; % usually 25%
perc_sup = 0.75; % and 75%

tick_label = {'.001', '.250', '.500', '.750', '.950', '.990', '.999'};
tick_value = [0.001 0.25 0.5 0.75 0.95 0.99 0.999];


ytick = y(round(tick_value*m));

xmi = res_s(round(m/100)); % 1% from the lower end.
xma = res_s(round(m - 4));   % Five samples from upper end.

xp25 = res_s(round(perc_inf*m));
xp75 = res_s(round(perc_sup*m));
yp25=icdf('chi2',perc_inf,resa(4));
yp75=icdf('chi2',perc_sup,resa(4));
% Passing a line thru both coordinates
a = (yp75 - yp25)/(xp75 - xp25);
b = yp25 - a*xp25;


ymi = a*xmi + yp25 - a*xp25;
if ymi < 0,
   xmi = (a*xp25 - yp25)/a;
   ymi = 0;
end;

yma = a*xma + yp25 - a*xp25;


% h1=plot(res_s(2:end-1), y(2:end-1), 'k+','MarkerSize',[8], 'MarkerEdgeColor',[.6 .6 .6],...
%    'MarkerFaceColor',[.6 .6 .6],'LineWidth',[1.5]);

% % % % ytick = y(round(tick_value*m));
% % % % 
% % % % % Linear  in new squared-root scale:
% % % % % Discard extreme values from both ends.
% % % % xmi = sqrt(res_s(round(m/100))); % 1% from the lower end.
% % % % xma = sqrt(res_s(round(m-5)));   % Five samples from upper end.
% % % % 
% % % % xp25 = sqrt(res_s(round(perc_inf * m)))
% % % % xp75 = sqrt(res_s(round(perc_sup * m)))
% % % % yp25 = sqrt(icdf('chi2', perc_inf, resa(4)))
% % % % yp75 = sqrt(icdf('chi2', perc_sup, resa(4)))
% % % % 
% % % % a = (yp75 - yp25)/(xp75 - xp25);
% % % % 
% % % % ymi = a*xmi + yp25 - a*xp25
% % % % yma = a*xma + yp25 - a*xp25
% % % % 
% % % % h1=plot(sqrt(res_s(2:end-1)), sqrt(y(2:end-1)), 'k+','MarkerSize',[8], ...
% % % %    'MarkerEdgeColor',[.6 .6 .6], 'MarkerFaceColor',[.6 .6 .6], ...
% % % %    'LineWidth',[1.5]);


h1=plot(sqrt(res_s(2:end-1)), sqrt(y(2:end-1)), 'k+','MarkerSize',[8], ...
   'MarkerEdgeColor',[.6 .6 .6], 'MarkerFaceColor',[.6 .6 .6], ...
   'LineWidth',[1.5]);
hold on
xlabel(['$\left|' measure '(\lambda{=.1)} \right|$'],'Interpreter', 'Latex', ...
   'fontsize', 14)
ylabel('(Weighted chi-square quantile)$^{1/2}$','Interpreter','Latex')

% % h2 = plot([sqrt(abs(xmi))-(-b/a),sqrt(abs(xma))-(-b/a)],[sqrt(abs(ymi)), ...
% %            sqrt(abs(yma))], 'r-','linewidth',[3]);
h2 = plot([sqrt(abs(xmi)),sqrt(abs(xma))],[sqrt(abs(ymi)), ...
           sqrt(abs(yma))], 'r-','linewidth',[3]);

plot(sqrt([xp25 xp75]), sqrt([yp25 yp75]), 'go-')
% plot([xp25 xp75], [yp25 yp75], 'b--')
%plot(sqrt([xmi xma]), sqrt([ymi yma]), 'bo-')

set(gca,'YTick',sqrt(ytick), 'YTickLabel',tick_label);
axis square
axis
save debugQqplotChi2_2

grid

% if min(res)/max(res) > 0.1
%    axis([sqrt(min(res)) sqrt(max(res)) sqrt(y(1)) sqrt(y(end))]);
% else
%    axis([0 sqrt(max(res)) sqrt(y(1)) sqrt(y(end))]);
% end;

%==========================================================================
function []=qqnorm(res, resa,measure)
%     qqplot for normal distribution

%     res contains sample values, resa asymp values (mean, var1, patden,
%     patdf)
%
m = size(res,1);
[res_s,inx] = sort(res);

x = linspace(1.0/m, 1-1.0/m, m);

% resa(1) == mean;  resa(2) := standard-deviation
y = icdf('norm',x,resa(1),resa(2));

tick_label = {'.001', '.010', '.050', '.250', '.750', '.950', '.990', '.999'};
tick_value = [0.001 0.01 0.05 0.25 0.75 0.95 0.99 0.999];
ytick = y(round(tick_value*m));

perc_inf = 0.25; % usually 25%
perc_sup = 0.75; % and 75%

% xmi = res_s(2);
% xma = res_s(m-1);
xmi = res_s(round(m/500));
xma = res_s(round(m-5));
%
xp25 = res_s(round(perc_inf*m));
xp75 = res_s(round(perc_sup*m));

yp25=icdf('norm', perc_inf, resa(1), resa(2));
yp75=icdf('norm', perc_sup, resa(1), resa(2));
%
a = (yp75 - yp25)/(xp75 - xp25);

ymi = a*xmi + yp25 - a*xp25;
yma = a*xma + yp25 - a*xp25;

h1=plot(res_s(2:end-1), y(2:end-1), 'k+','MarkerSize',[8], 'MarkerEdgeColor',[.6 .6 .6],...
   'MarkerFaceColor',[.6 .6 .6],'LineWidth',[1.5]);

xlabel(['$\left|' measure '(\lambda{=.1)} \right| ^{2}$'],'Interpreter', 'Latex', 'fontsize', 14)
ylabel('Normal quantile')

h2 = plot([xmi,xma],[ymi,yma], 'k-','linewidth',[3]);
%plot([xp25 xp75], [yp25 yp75], 'ro')

plot([xp25 xp75], [yp25 yp75], 'go')
plot([xp25 xp75], [yp25 yp75], 'b--')
plot([xmi xma], [ymi yma], 'bo')

set(gca,'YTick',ytick, 'YTickLabel',tick_label);
grid
if min(res)/max(res) > 0.1
   axis([min(res) max(res) y(1) y(end)]);
else
   axis([0 max(res) y(1) y(end)]);
end;
