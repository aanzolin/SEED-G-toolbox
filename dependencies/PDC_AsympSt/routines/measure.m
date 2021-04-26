function c = measure(c)
% Returns the type of result contained in c structure, either DTF or PDC.
%function c = measure(c)
%       Returns the type of result contained in c structure, either DTF or
%       PDC.

if isfield(c,'dtf')
   c = 'DTF';
elseif isfield(c,'pdc')
   c = 'PDC';
else
   error('Structure is not a PDC or DTF analysis results.')
end;
   