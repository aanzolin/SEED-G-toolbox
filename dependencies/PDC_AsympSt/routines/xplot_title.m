function [] = xplot_title(alpha,metric,measure)
%function [] = xplot_title(alpha,metric,measure)
%
% measure should be either 'pdc'| 'dtf'

if nargin == 1,
    c=alpha;
    if isfield(c,'dtf')
        measure = 'DTF';
        alpha = c.alpha;
        metric =  c.metric;
    elseif isfield(c,'pdc')
        measure = 'PDC';
        alpha = c.alpha;
        metric =  c.metric;
    end
elseif nargin < 3,
    measure = 'PDC';
end;
alphastr = int2str(100*alpha);
measure = upper(measure);
metric = lower(metric);
switch metric
    case 'euc'
        
    case 'diag'
        if measure == 'DTF',
            measure = 'DC';
        else
            measure = '_gPDC';
        end;
    case 'info'
        measure = ['_i' measure];
    otherwise
        error('Unknown metric.')
end;

suptitle(['{|' measure '|}^2 (' '{\alpha = ' alphastr '%}' ')'])