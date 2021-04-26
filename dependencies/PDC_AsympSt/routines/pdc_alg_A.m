function pdc = pdc_alg_A(A, pf, nf, metric)
%function pdc = pdc_alg_A(A, pf, nf, metric)
%     '''Generates spectral general (estatis. norm) Squared PDC matrix from AR matrix
%
%       Input:
%         A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
%         pf(n, n) - error covariance matrix
%         nf - frequency resolution
%
%       Output:
%         PDC2(n, n, nf) - PDC^2 matrix

[n,n,r] = size(A);
switch lower(metric)
    case 'euc'
        nornum = ones(n,1);
        norden = eye(n);
    case 'diag'
        nornum = 1./diag(pf);
        norden = diag(1./diag(pf));
    case 'info'
        nornum = 1./diag(pf);
        norden = inv(pf);
    otherwise
        error('Unknown metric.')
end;

Af = A_to_f(A, nf); % [ nf x n x n ]
pdc=zeros(n,n,nf);

for ff = 1:nf,
    Aff=getAff(Af,ff);
    for kj=1:n,
        Affj=Aff(:,kj).*sqrt(nornum);
        pdc(:,kj,ff) = Affj./(sqrt(abs((Aff(:,kj)')*norden*Aff(:,kj))));
    end;
end;

pdc = abs(pdc).^2;

function c = getAff(Af,ff)
%function c = getAff(C,ff)
% Input:      C [NumChannel, NumChannel, nFreqs], either PDC,Lpatnaik
%               Lv2inf, Lv2sup
%             ff - ff-th frequency component of C
% Output:     c - A[ff,:,:] element

[Nfreq,Nch,Nch] = size(Af);
c=reshape(Af(ff,:,:), Nch,Nch);
