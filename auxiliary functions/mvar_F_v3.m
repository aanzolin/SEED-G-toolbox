function [ARF,RCF,PE,DC] = mvar_F_v2(Y, Pmax, Mode)
% Estimates Multi-Variate AutoRegressive model parameters 
% function  [AR,RC,PE] = mvar(Y, Pmax);
%
% INPUT:
%   Y	Multivariate data series 
%   Pmax 	Model order
%
% OUTPUT
%   AR    multivariate autoregressive model parameter (same format as in [4]	
%   RC    reflection coefficients (= -PARCOR coefficients)
%   PE    remaining error variance
%
% All input and output parameters are organized in columns, one column 
% corresponds to the parameters of one channel.
%
% A multivariate inverse filter can be realized with 
%       [AR,RC,PE] = mvar(Y,P);
%	e = mvfilter([eye(size(AR,1)),-AR],eye(size(AR(1))),Y);
%
% see also: MVFILTER, COVM, SUMSKIPNAN, ARFIT2
%
% REFERENCES:
%  [1] M.S. Kay "Modern Spectral Estimation" Prentice Hall, 1988. 
%  [2] S.L. Marple "Digital Spectral Analysis with Applications" Prentice Hall, 1987.
%  [3] M. Kaminski, M. Ding, W. Truccolo, S.L. Bressler, Evaluating causal realations in neural systems:
%	Granger causality, directed transfer functions and statistical assessment of significance.
%	Biol. Cybern., 85,145-157 (2001)
%  [4] T. Schneider and A. Neumaier, A. 2001. 
%	Algorithm 808: ARFIT-a Matlab package for the estimation of parameters and eigenmodes 
%	of multivariate autoregressive models. ACM-Transactions on Mathematical Software. 27, (Mar.), 58-65.
%  [5] A. Schlogl 2002. 
%	Validation of MVAR estimators or Remark on Algorithm 808: ARFIT, 
%	ACM-Transactions on Mathematical Software. submitted.

%	Copyright (C) 1996-2002 by Alois Schloegl <a.schloegl@ieee.org>	

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

%% Modified by alessandra anzolin for multi-trial input
% Initialization original 
%%% [N,M] = size(Y);
% Initialization for SEED-G: specify trial 1 because the dataset can be
% multi-trial
[N,M] = size(Y{1}); 

if nargin<2
        Pmax=max([N,M])-1;
end

if iscell(Y)
        Pmax = min(max(N ,M ),Pmax);
        clear C;
else
    %%%%% Estimate Autocorrelation funtion 	
    if 0
        [tmp,LAG]=xcorr(Y,Pmax,'biased');
        for K=0:Pmax
            %C{K+1}=reshape(tmp(find(LAG==K)),M ,M );	
            C(:,K*M+(1:M))=reshape(tmp(find(LAG==K)),M ,M );	
        end
    else
        for K =0:Pmax
            %C{K+1}=Y(1:N-K,:)'*Y(K+1:N ,:)/N ;
            %C{K+1}=Y(K+1:N,:)'*Y(1:N-K,:)/N; % =Rxx(-k)=conj(Rxx(k)) in [2] with K=k+1; 
        end
    end
end
if nargin<3
    % tested with a bootstrap validation, Mode 2 or 5 are recommended
    %Mode=5;  % M*P << N
    %Mode=5;  % 5*6 << 100, test=900, permutations 1000
    Mode=2;   % M*P ~~ N
end

TotTrial=length(Y);
            [C(:,1:M),n] = covm(Y{1},'M');
            temp=zeros*C;
            for NumTrial=2:TotTrial
                temp=temp+C;
                [C(:,1:M),n] = covm(Y{NumTrial},'M');
            end
            temp=temp+C;    
            C=temp./TotTrial;
            PE(:,1:M)  = C(:,1:M)./n;
            
            if Mode==2
                %%%%% multi-channel Levinsion algorithm 
                %%%%% using Nutall-Strand Method [2]
                %%%%% Covariance matrix is normalized by N=length(X)-p 
                C(:,1:M) = C(:,1:M)/(N);
                F = Y;
                B = Y;
                PEF = C(:,1:M);
                PEB = C(:,1:M);
                
                for K=1:Pmax
                    D=covm(F{1}(K+1:N,:),B{1}(1:N-K,:),'M');

                    for NumTrial=2:TotTrial
                        D=D+covm(F{NumTrial}(K+1:N,:),B{NumTrial}(1:N-K,:),'M');
                    end 
                    D=D./TotTrial;
                    ARF(:,K*M+(1-M:0)) = D / PEB;	
                    ARB(:,K*M+(1-M:0)) = D'/ PEF;	
                    
                    for NumTrial=1:TotTrial
                        tmp = F{NumTrial}(K+1:N,:) - B{NumTrial}(1:N-K,:)*ARF(:,K*M+(1-M:0)).';
                        B{NumTrial}(1:N-K,:) = B{NumTrial}(1:N-K,:) - F{NumTrial}(K+1:N,:)*ARB(:,K*M+(1-M:0)).';
                        F{NumTrial}(K+1:N,:) = tmp;
                    end
                    
                    for L = 1:K-1
                        tmp      = ARF(:,L*M+(1-M:0))   - ARF(:,K*M+(1-M:0))*ARB(:,(K-L)*M+(1-M:0));
                        ARB(:,(K-L)*M+(1-M:0)) = ARB(:,(K-L)*M+(1-M:0)) - ARB(:,K*M+(1-M:0))*ARF(:,L*M+(1-M:0));
                        ARF(:,L*M+(1-M:0))   = tmp;
                    end
                    
                    RCF(:,K*M+(1-M:0)) = ARF(:,K*M+(1-M:0));
                    RCB(:,K*M+(1-M:0)) = ARB(:,K*M+(1-M:0));
                    
                    PEF=zeros*PEF;
                    PEB=zeros*PEB;
                    
                    for NumTrial=1:TotTrial
                        PEF = PEF+covm(F{NumTrial}(K+1:N,:),F{NumTrial}(K+1:N,:),'M');
                        PEB = PEB+covm(B{NumTrial}(1:N-K,:),B{NumTrial}(1:N-K,:),'M');
                    end

                    PEF=PEF./(TotTrial);
                    PEB=PEB./(TotTrial);
                    PE(:,K*M+(1:M)) = PEF;        
                end
                
            end

MAR    = zeros(M,M*Pmax);
DC     = zeros(M);

for K  = 1:Pmax
        DC = DC + ARF(:,K*M+(1-M:0))'.^2; %DC meausure [3]
end

