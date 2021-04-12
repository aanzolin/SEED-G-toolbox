function [ModelDel] = rearrangeModel(DelayRange,Model,DelayMatrix)
% Function that returns a 3D Model with the connections distributed across 
% the delays
% Created on June 5 2018
%% @author: manuela petti (manuela.petti@uniroma1.it)

minDel=min(DelayRange);
maxDel=max(DelayRange);
i=0;

for v=minDel:maxDel
    i=i+1;
    ind=find(DelayMatrix==v);
    MD=zeros(size(Model,1),size(Model,2));
    MD(ind)=Model(ind);
    ModelDel(:,:,i)=MD;
    
    clear MD ind
end
