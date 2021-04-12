function [Model, DelayMatrix, MaxDelay] = get_ConnectivityModel(ModelStructure, ValuesRange, MinDelta, DelayRange)
%%% Input:
%       - ModelStructure: binary connectivity model.  
%                         NxN matrix (N: number of nodes) where elements equal to 1 
%                         correstpond to non-null connections. Every other element is
%                         equal to 0.
%       - ValuesRange:    range of possible connections values
%       - MinDelta:       minimum difference between two different values of connections
%       - DelayRange:     range of possible delay values
%%% Output:
%       - Model:          weighted connectivity model.
%                         NxN matrix where the strength of non-null
%                         connections is randomly assigned from the range
%                         'ValuesRange'
%       - DelayMatrix:    NxN matrix in which the entry ij (i,j=1,...N) 
%                         corresponds to the lag associated with the connection
%                         ij. Delay values are randomly selected within the range
%                         DelayRange.
%       - MaxDelay:       maximum possible delay
%
% Created on June 5 2018
% Modified on May 29 2019
%       to allow negative connection values in ValuesRange
%% @author: manuela petti (manuela.petti@uniroma1.it)

Structure = ModelStructure;
N_conn = length(find(ModelStructure));
nNodes = size(Structure,1);

%%% assignment of values to the connections
ModelValues = min(ValuesRange):MinDelta:max(ValuesRange);

if any(ModelValues==0)
    x=find(ModelValues==0);
    ModelValues(x)=[];
end

ModelValues_pos = ceil((rand(1,N_conn))*length(ModelValues));
ModelStructure(find(ModelStructure)) = ModelValues(ModelValues_pos);
Model = ModelStructure;

%%% Delay values assignment
DelayValues = min(DelayRange):1:max(DelayRange);
DelayValues_pos = ceil((rand(1,N_conn))*length(DelayValues));
Structure(find(Structure)) = DelayValues(DelayValues_pos);
DelayMatrix = Structure;

%%% Max delay 
diagSecondarie = zeros([nNodes-1,1]);
for rr=1:nNodes
    for cc=(rr+1):nNodes
        diagSecondarie(cc-rr) = diagSecondarie(cc-rr)+DelayMatrix(rr,cc);
    end
end
MaxDelay = max(diagSecondarie);
