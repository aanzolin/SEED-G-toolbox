%%% Modified verison of the function makerandCIJ_dir from the
% Brain Connectivity Toolbox [1] available at https://sites.google.com/site/bctnet/
% 
% References:
% [1] Complex network measures of brain connectivity: Uses and interpretations.
% Rubinov M, Sporns O (2010) NeuroImage 52:1059-69.
%
%%% Modified by Alessandra Anzolin in 'makerandCIJ_dir_realARcomponents'
%       1. control on Auto-Regressive Coefficients 
%       2. new output added 
%          inPos --> vector containinf the position of the AR coefficients
%          in the MVAR connectivity model

function [CIJ, inPos] = makerandCIJ_dir_realARcomponents(N,K,ARn)
% inputs:
%           N = number of vertices
%           K = number of edges
%           ARn = number of non-zero connection on the main diagonal
% output:
%           CIJ = directed random connection matrix
%
% Generates a random directed binary connection matrix, with size (N,K) 
% with ARn non-zero connection on the main diagonal (autoregressive part)
% Olaf Sporns, Indiana University, 2007/2008

ind = ~eye(N);
nPos = randperm(N);
inPos = nPos(1:ARn);

for ar=1:ARn
    posD = inPos(ar);

    ind(posD,:)=0;
    clear posD
end

i = find(ind==1);

rp = randperm(length(i));
irp = i(rp);

CIJ = zeros(N);
CIJ(irp(1:K)) = 1;
