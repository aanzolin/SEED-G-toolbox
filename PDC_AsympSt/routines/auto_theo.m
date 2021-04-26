function [ R ] = auto_theo(A,Sigma_w )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% input: A model structure
%        Sigma_w

[K,K,p]=size(A);
Ar=reshape(A,K,K*p);
At=eye(K*(p-1),K*p);
Az=[Ar;At];
Sigma=zeros(K*p,K*p);
Sigma(1:K,1:K)=Sigma_w;
s=reshape(Sigma,K*K*p*p,1);
R=inv(eye(K*K*p*p)-kron(Az,Az));
R=reshape(R*s,K*p,K*p);

end

