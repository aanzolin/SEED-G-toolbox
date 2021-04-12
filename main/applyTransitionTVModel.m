function [TV_mod, pos]=applyTransitionTVModel(model,var_perc,var_val)

% Positive transition
indpos=find(abs(model)<0.25 & model~=0); %Vincolo: la connessione incrementata non deve superare 0.7 (circa) di valore

NtransPos=ceil(var_perc*length(indpos));
indpos1=randperm(length(indpos));
if ~isempty(indpos1)
    indpos=indpos(indpos1(1:NtransPos));
else
    indpos=[]; % no connections to increase 
end

% Negative transition
indneg = find(abs(model)>0.25);
NtransNeg=ceil(var_perc*length(indneg));
indneg1=randperm(length(indneg));

if ~isempty(indneg1)
    indneg=indneg(indneg1(1:NtransNeg));
else
    indneg=[]; 
end

pos = [indpos;indneg];

TV_mod = model;
TV_mod(indpos)=roundn((var_val+0.1*rand)*model(indpos),-2);
TV_mod(indneg)=roundn((1/(var_val+0.1*rand))*model(indneg),-2);

clear mod

