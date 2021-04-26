function AL = A_to_f(A, nf)
%    '''Calculates A(f), in the frequency domain.
%
%     Input:
%         A(nChannels, nChannels, p) - recurrence matrix (nChannels - number of signals,
%                                      p - model order)
%         nf - frequency resolution
%
%     Output:
%         AL(nf, nChannels, nChannels)
%     '''

[nChannels,~,p] = size(A); % Dummy variable ~ is equivalent to nChannels
Jimag=sqrt(-1);
%   exponents are an array of FFT exponents, on all frequency range for
%   each lag.
%   exponents = reshape((-Jimag*pi*kron(0:(nf-1),(1:p))/nf),nf,p);
exponents = reshape((-Jimag*pi*kron(0:(nf-1),(1:p))/nf),p,nf).';
%    Af faz a multiplicacoes da exp(ar) pela matrix A, para todas frequencias
%    as funcoes repeat e transpose sao truques para possibilitar o calculo vetorial
Areshaped=reshape(A, nChannels,nChannels,1,p);

%Af = (Areshaped.repeat(nf, axis=3))
Af=zeros(nChannels,nChannels,nf,p);
for kk=1:nf
   Af(:,:,kk,:)=Areshaped;
end;
for i=1:nChannels,
   for k=1:nChannels,
      Af(i,k,:,:)=reshape(Af(i,k,:,:),nf,p).*exp(exponents);
   end;
end;

Af=permute(Af, [3,1,2,4]);
%    o fft soma o valor para todos os lags
AL=zeros(nf,nChannels,nChannels);
for kk=1:nf,
   temp=zeros(nChannels,nChannels);
   for k=1:p
      temp = temp+reshape(Af(kk,:,:,k),nChannels,nChannels);
   end;
   temp=eye(nChannels)-temp;
   AL(kk,:,:) = reshape(temp,1,nChannels,nChannels);
end;
