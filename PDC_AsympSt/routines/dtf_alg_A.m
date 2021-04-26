function dtf = dtf_alg_A(A, pf, nf, metric)
%function dtf = dtf_alg_A(A, pf, nf = 64, metric = 'gen'),
%     '''Generates spectral general (estatis. norm) DTF matrix from AR matrix
%
%       Input:
%         A(nChannels, nChannels, r) - recurrence matrix (nChannels - number of signals, r - model order)
%         pf(nChannels, nChannels) - error covariance matrix
%         nf - frequency resolution
%
%       Output:
%         DTF(nChannels, nChannels, nf) - DTF matrix
%     '''

% DTF calculation extracted from asymp_dtf2.m routine.

[nChannels nChannels r] = size(A);% [nChannels,n0,p] = size(A);
Af = A_to_f(A, nf); % [ nf x nChannels x nChannels ]

dtf=zeros(nChannels,nChannels,nf);
evar_d = mdiag(pf);
evar_d_big = kron(eye(2),kron(evar_d,eye(nChannels)));
evar_big = kron(eye(2),kron(pf,eye(nChannels)));

for ff = 1:nf,
    f = (ff-1)/(2*nf); % Corrected 7/25/2011, f starting at 0.
    Af_ff = reshape(Af(ff,:,:),[nChannels, nChannels]);
    Hf = pinv(Af_ff); h = Hf(:);  % Equivalent to h = vec(Af[ff, :, :].I)
    
    h = [real(h); imag(h)];    % h = cat(h.real, h.imag, 0)
    for i = 1:nChannels,
        for j = 1:nChannels,
            Iij = fIij(i,j,nChannels);
            Ii = fIi(i,nChannels);
            switch lower(metric)
                case {'euc'}               % for DTF
                    Iije = Iij;
                    Iie  = Ii;
                case {'diag'}              % for DC
                    Iije = Iij*evar_d_big*Iij;
                    Iie  = Ii*evar_d_big*Ii;
                case {'info'}              % for iDTF
                    Iije = Iij*evar_d_big;
                    %Iie  = Ii*evar_big*Ii;
                    Iie  = Ii*evar_big;
                otherwise
                    error('Unknown metric.')
            end;
            num = h.'*Iije*h;
            den = h.'*Iie*h;
            dtf(i,j,ff) = num/den;
        end;
    end;
end;

function c = getAff(Af,ff)
%function c = getAff(C,ff)
%
% Input:      C [NumChannel, NumChannel, nf], either PDC,Lpatnaik
%               Lv2inf, Lv2sup
%             ff - ff-th frequency component of C
% Output:     c - A[ff,:,:] element

[Nfreq,Nch,Nch] = size(Af);
c = reshape(Af(ff,:,:),Nch,Nch);

%==========================================================================
function c = fIij(i,j,n)
%'''Returns Iij of the formula'''
Iij = zeros(1,n*n);
Iij(n*(j-1)+i) = 1;
Iij = diag(Iij);
c = kron(eye(2), Iij);

%==========================================================================
function c = fIi(i,n)
%    '''Returns Ii of the formula'''
Ii = zeros(1,n);
Ii(i) = 1;
Ii = diag(Ii);
Ii = kron(eye(n), Ii);
c = kron(eye(2), Ii);
%==========================================================================
function c = mdiag(a)
%  diagonal matrix
c = diag(diag(a));
