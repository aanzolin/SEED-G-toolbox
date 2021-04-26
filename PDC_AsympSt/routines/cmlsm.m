function [npf,na,nef]=cmlsm(u,IP)
% MVAR least squares estimator
%
%function [npf,na,nef]=cmlsm(u,IP);
%
%input:   u   - vector of rows
%         IP  - order
% output: npf - error covariance
%         na  - model
%         nef - residue

[m,n]=size(u);
[b,SU,e]=mlsmx(u,IP);
na=reshape(b,m,m,IP);
npf=SU*n; % see normalization
nef=e;

%==========================================================================
%
% 01/30/1998 - L.A.B. 11/4/2000 (LAB reviewed)
%
function [b,SU,e] = mlsmx(Y,p)
[K,T] = size(Y);
Z = zmatrm(Y,p);
Gamma = Z*Z';
U1 = Gamma\Z; % It is equivalent to inv(Gamma)*Z;
% Gamma = Gamma/T;  % not used (???)
SU = (Y*Y'-Y*Z'*U1*Y');
SU = SU/(T-K*p-1);
b = kron(U1,eye(K))*reshape(Y,K*T,1);
e = reshape(reshape(Y,K*T,1)-kron(Z',eye(K))*b,K,T); 

%==========================================================================
% Computation of Z - data structure (no estimation of the mean)
%
% function Z=zmatr(Y,p);
%
% input:  Y - data in row vectors 
%         p - model covariance order
%
% output: Z
%
% 01/30/1998 - L.A.B.
%
function Z=zmatrm(Y,p)
[K,T] = size(Y);
y1 = [zeros(K*p,1);reshape(flipud(Y),K*T,1)];
Z =  zeros(K*p,T);
for i=0:T-1
   Z(:,i+1)=flipud(y1(1+K*i:K*i+K*p));
end
