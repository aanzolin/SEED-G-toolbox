function Coh = coh_alg(SS)

[m,n,nFreqs] = size(SS);
Coh=zeros(size(SS));
if m == n, nChannels = m; else error('Wrong SS dimension.'); end;

for k=1:nFreqs,
    for iu=1:nChannels
        for ju=1:nChannels
            Coh(iu,ju,k)=SS(iu,ju,k)./sqrt(SS(iu,iu,k).*SS(ju,ju,k));
        end
    end
end;
