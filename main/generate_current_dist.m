%% Function "generate_current_dist" 
%  Compute J: Generates cortical current density distribution.
%  
%  Created on  March 2 2019
%% @author: jlenia toppi 

function [J]= generate_current_dist(Y_gen,sa,type,source_ind)

% INPUT:
%       Y_gen: samp x ROI (simulated time series)
%       sa: struct containing head model (lead field)
%       type: 'dist' - gaussian distribution around the centroide
%             'one' - centroides
%       source_ind: dipoles of interest

%OUTPUT:
%       J: current density distribution 

[samp, source_num, trial] = size(Y_gen);
sigma_range =[10 40];
sigma = sigma_range(1) + diff(sigma_range)*rand(source_num, 1);
J=zeros(samp,size(sa.pos,1),trial);

switch type
    case 'one'
        J(:,source_ind,:)=Y_gen;
    case 'dist'
        cortex75K.vc = sa.pos;
        cortex75K.tri = sa.tri;
        for ss=1:source_num
            %%% generate Gaussian distributions on cortical manifold
            [~, source_amp(:, ss)] = graphrbf(cortex75K, sigma(ss), source_ind(ss));
            J(:,source_ind)=Y_gen;
            for i=1:size(J,1)
                J(i,:,ss)=source_amp'*Y_gen(i,ss);
            end
        end
end

J=multitransp(J,1);
end