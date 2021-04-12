%%% Input
%       mod :           connectivity model (binary, 2D)
%       num_tr :        tot trials number
%       perc_tr: percentage of trial with variable/spurious connections
%               Range [0 1]
%       perc_sp: percentage of spurious connection
%               Range [0 1]
%       perc_con: percentage of connections variable across the trails
%               Range [0 1]
%       perc_modValue: the entity of the variation as percentage-change with 
%               respect to original weight of the connection
%               Range [0 1]
%       mod_sig: direction (increase or decrease with respect to the 
%               original weight) of the variation       
%               1 for positive changes 
%               -1 for negative changes
%%% Output
%       mat_rep: altered model (in a specified percentage of trials)
%       flag : 1 if the total number of spurious connections to inster in
%               higher than the number of null elements in the model.
%
%  Created  on June 15 2018
%% @author: alessandra anzolin (aanzolin@mgh.harvard.edu)
%
%  Modified on June 25 2018 by alessandra anzolin
%       - output "delayMat_rep" added
%
% Modified on March 23 2019 by manuela petti
%       - in the spurious model, the value of perc_sp can be negative (some of the
%         real connection disappear for some of the trials)
%       - max delay is computed trial by trial. Output "max_maxDel" added.  
%       - flag added 
%
% Modified on November 23 2019 by alessandra anzolin
%       - input "mod_Sig" added


function [mat_rep, delayMat_rep, max_maxDel, flag] = get_ConnectivityModel_withITV(mod, num_tr, perc_tr, ...
    perc_sp, perc_con, mod_con, mod_sig, ValuesConn_Range, delay_mat, DelayRange,inDiag)

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
mat_rep =       [];
delayMat_rep =  [];
max_maxDel =    [];
tr_mod =   round((num_tr*perc_tr)); 
tr_rep =   num_tr - tr_mod;            
n_con =    length(find(mod~=0));        
n_sp =     round((n_con*perc_sp))+1;  
n_mod =    round((n_con*perc_con))+1; 

%%% Replicate model for each trial
mod_rep =  repmat(mod,[1 1 tr_rep]);
%%% Replicate delay matrix for each trial
del_rep =  repmat(delay_mat,[1 1 tr_rep]);

nch = size(mod,1);
diagSecondarie = zeros([nch-1,1]);

for rr=1:nch
    for cc=(rr+1):nch
        diagSecondarie(cc-rr) = diagSecondarie(cc-rr)+delay_mat(rr,cc);
    end
end
maxDel_rep = repmat(max(diagSecondarie),[1 tr_rep]);

%%% Loop on trial to alter
for i = 1:tr_mod
    diagSecondarie = zeros([nch-1,1]);
    mat = mod;
    del = delay_mat;
    
    %%% Case1: 'Modified Connections'
    if n_mod>0
        ind = find(mat~=0);
        ind_rand = randperm(length(ind));
        ind_n = ind(ind_rand(1:n_mod)); % estrarre i primi "n_mod" indici (ormai disposti casualmente)
        num= mod_sig;
        mat(ind_n) = mat(ind_n) + num.*((mat(ind_n)*mod_con));
        clear ind ind_n ind_rand num
        flag = [];
    end
    
    %%% Case2: 'Spurious Connections'
    if n_sp~=0
        
        if n_sp>0
            mat_d=mat;
            mat_d(find(eye(nch)))=NaN;
            mat_d(inDiag,:)=NaN;
            NewPossible = find(mat_d==0); % selezionare tutti gli indici 
            
            if length(NewPossible)<n_sp
                flag(i) = 1;
                return
            else
                flag(i) = 0;
            end
            
            New_mixed = randperm(length(NewPossible));
            newselected = NewPossible(New_mixed(1:n_sp));
            
            ModelValues = min(ValuesConn_Range):0.1:max(ValuesConn_Range); % 0.1 is the minimum difference between two different values of connections
            mat(newselected) = ModelValues(randi([1,length(ModelValues)],n_sp,1));
            
            DelayValues = min(DelayRange):1:max(DelayRange);
            DelayValues_pos = ceil((rand(1,n_sp))*length(DelayValues));
            del(newselected) = DelayValues(DelayValues_pos);
            
        elseif n_sp<0
            ind_conn = find(mat~=0);
            tic
            while ~exist('modified_mat')
                modified_mat = mat;
                selected = ind_conn(randi(length(ind_conn),1,abs(n_sp)));
                modified_mat(selected) = 0;
                
            end
            mat = modified_mat;
            del(selected) = 0;
            flag = [];
        end
    end
    
    mat_mod(:,:,i)=mat;
    mat_del(:,:,i)=del;
    clear mat del
    
    for rr=1:nch
        for cc=(rr+1):nch
            diagSecondarie(cc-rr) = diagSecondarie(cc-rr)+mat_del(rr,cc,i);
        end
    end
    maxDel_rep(tr_rep+i) = max(diagSecondarie);
end

max_maxDel = max(maxDel_rep);
mat_rep = cat(3,mod_rep,mat_mod);
delayMat_rep = cat(3,del_rep,mat_del);
