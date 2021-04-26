function c=asymp_dtf_theo(x,A,pf,nFreqs,metric,alpha)
%Compute DTF connectivity measure given by "option" from series j-->i.
%
%function c=asymp_dtf_theo(x,A,pf,nFreqs,metric,alpha)
%
% input: x - data
%        A - AR estimate matrix by MVAR
%        pf - covariance matrix provided by MVAR
%        nFreqs - number of point in [0,fs/2] frequency scale
%        metric   euc  - Euclidean     ==> original DTF
%                 diag - diagonal      ==>      DC
%                 info - information   ==>     iDTF
%        alpha = .05 default for distribution
%                if alpha = zero, do not calculate statistics
% output:  c.dtf       - |DTF|^2 estimates
%          c.th        - Threshold value with (1-avalue) significance level.
%          c.ic1,c.ic2 -  confidence interval
%          c.metric    - metric for DTF calculation
%          c.alpha     - significance level
%          c.p         - VAR model order
%          c.patdenr   -
%          c.patdfr    -
%function c=asymp_dtf_theo(x,A,pf,nFreqs,metric,alpha)

%The asymp_pdc routine, from which asymp_dtf_theo is derived, was corrected
% on 7/25/2011 to match the frequency range with plotting routine, f = 0
% was included in the frequency for-loop:
%                                for ff = 1:nFreqs,
%                                   f = (ff-1)/(2*nFreqs); %
%                                        ^?^^

if nargin<6,
   error('ASYMP_DTF_THEO requires six input arguments.')
end
[m,n]=size(x);
if m > n,
   x=x.';
end;
np = length(x);
[nChannels,n0,p] = size(A);
Af = A_to_f(A, nFreqs);

% Variables initialization
dtf = zeros(nChannels,nChannels,nFreqs);
disp('----------------------------------------------------------------------');
if alpha ~= 0,
   th   = zeros(nChannels,nChannels,nFreqs);
   ic1 = zeros(nChannels,nChannels,nFreqs);
   ic2 = zeros(nChannels,nChannels,nFreqs);
   varass1 = zeros(nChannels,nChannels,nFreqs);
   varass2 = zeros(nChannels,nChannels,nFreqs);
   patdfr = zeros(nChannels,nChannels,nFreqs);
   patdenr = zeros(nChannels,nChannels,nFreqs);
   switch lower(metric)
      case {'euc'}
         disp('               Original DTF and asymptotic statistics')
      case {'diag'}
         disp('              DC and asymptotic statistics')
      case {'info'}
         disp('             Information DTF and asymptotic statistics')
      otherwise
         error('Unknown metric.')
   end;
else
   switch lower(metric)
      case {'euc'}
         disp('                       Original DTF estimation')
      case {'diag'}
         disp('                      Generalized DTF or DC estimation')
      case {'info'}
         disp('                     Information DTF estimation')
      otherwise
         error('Unknown metric.')
   end;
end;
disp('======================================================================');

%gamma = bigautocorr(x, p);
gamma=auto_theo(A,pf);

omega = kron(pinv(gamma), pf);
%omega_evar = 2*pinv(Dup(nChannels))*kron(pf, pf)*pinv(Dup(nChannels)).'
omega_evar = 2*Dup(nChannels)*pinv(Dup(nChannels))*kron(pf, pf) ...
                            *(pinv(Dup(nChannels)).')*Dup(nChannels).';
for ff = 1:nFreqs,
   f = (ff-1)/(2*nFreqs); % Corrected 7/25/2011, f starts at 0.
   Ca = fCa(f, p, nChannels);

   Af_ff = reshape(Af(ff,:,:),[nChannels, nChannels]);
   Hf = pinv(Af_ff); h = Hf(:);  % Equivalent to h = vec(Af[ff, :, :].I)

   h = [real(h); imag(h)];    % h = cat(h.real, h.imag, 0)
   H = fdh_da(Af_ff);         % = ha; H = fdh_da(mat(Af[ff, :, :]), n)

   Omega_h = H*Ca*omega*Ca.'*H.';  % \Omega_h
   L = fChol(Omega_h); % real part only

   for i = 1:nChannels,
      for j = 1:nChannels,
         Iij = fIij(i,j,nChannels);
         Ii = fIi(i,nChannels);
         switch lower(metric)
            case {'euc'}               % for DTF
               Iije = Iij;
               Iie  = Ii;
               %evar_d_big = 1;

            case {'diag'}              % for DC
               evar_d = mdiag(pf);
               evar_d_big = kron(eye(2),kron(evar_d,eye(nChannels)));
               Iije = Iij*evar_d_big*Iij;
               Iie  = Ii*evar_d_big*Ii;

            case {'info'}              % for iDTF
               evar_d = mdiag(pf);
               evar_d_big = kron(eye(2),kron(evar_d,eye(nChannels)));
               evar_big = kron(eye(2),kron(pf,eye(nChannels)));
               %Iie  = Ii*evar_big*Ii;
               Iije = Iij*evar_d_big;
               Iie  = Ii*evar_big;

            otherwise
               error('Unknown metric.')
         end;

         num = h.'*Iije*h;
         den = h.'*Iie*h;
         dtf(i,j,ff) = num/den; % |DTF|^2
         
         % If alpha == 0, do not calculate statistics for faster DTF
         % computation.
         if alpha ~= 0,
            %'Acrescenta derivada em relacao a evar'
            switch lower(metric)
               case {'euc'}               % for DTF
                  ddtf_dev = zeros(1,nChannels^2);

               case {'diag'}              % for DC
                  %#todo: tirar partes que nao dependem de f do loop.
                  if (i == 1) && (j == 1) && (ff == 1),
                     %'derivada de vec(Ed-1) por vecE'
                     % de_deh = Dup(nChannels);
                     debig_de = fdebig_de_dtf(nChannels); %New Theta_K
                     dedinv_deh = debig_de*diag(vec(eye(nChannels)));
                  end;
                  dnum_dev = kron((Iij*h).', h.'*Iij)*dedinv_deh;
                  %'derivada do den por vecE'
                  dden_dev = kron((Ii*h).', h.'*Ii)*dedinv_deh;
                  ddtf_dev = (den*dnum_dev - num*dden_dev)/(den^2);

               case {'info'}              % for iDTF
                  if i == 1 && j == 1 && ff == 1,
                     %'derivada de vec(Ed-1) por vecE'
                     debig_de = fdebig_de_dtf(nChannels); %New Theta_K
                     dedinv_deh = debig_de*diag(vec(eye(nChannels)));
                     
                  end;
                  %'derivada do num por vecE'
                  dnum_dev = kron((Iij*h).', h.'*Iij) * dedinv_deh;
                  %'derivada do den por vecE'
                  dden_dev = kron((Ii*h).', h.'*Ii) * debig_de;
                  ddtf_dev = (den*dnum_dev - num*dden_dev)/(den^2);

               otherwise
                  error('Unknown metric.')
            end;

            G1h = 2*h.'*Iije/den - 2*num*h.'*Iie/(den^2); % Eq. (15)
            %   G1 = -G1h*H*Ca;                               % (Cont.)
            varalpha = G1h*Omega_h*G1h.';
            varevar = ddtf_dev*omega_evar*ddtf_dev.';
            varass1(i,j,ff) = (varalpha + varevar)/np; % Eq. (14)

            ic1(i,j,ff) = dtf(i,j,ff) ...
                          - sqrt(varass1(i,j,ff))*norminv(1-alpha/2.0,0,1);
            ic2(i,j,ff) = dtf(i,j,ff) ...
                          + sqrt(varass1(i,j,ff))*norminv(1-alpha/2.0,0,1);

            G2h = Iije/den; % 

            d = fEig(real(L), real(G2h)); % real() 28May2013

            patdf = (sum(d).^2)./sum(d.^2);
            patden = sum(d)./sum(d.^2);
            if (is_octave)
               th(i,j,ff) = chi2inv((1-alpha), patdf)/(patden*np);
               pvalues(i,j,ff) = 1 - chi2cdf(dtf(i,j,ff)*patden*np, patdf);
            else
               th(i,j,ff) = icdf('chi2',(1-alpha), patdf)/(patden*np);
               pvalues(i,j,ff) = 1 - cdf('chi2', dtf(i,j,ff)*patden*np, patdf);
            end;

            varass2(i,j,ff) = patdf/(patden*np).^2;
            patdfr(i,j,ff) = patdf;
            patdenr(i,j,ff) = patden;
         else % as alpha == 0, do not compute asymptotics
            %nop
         end;
      end;
   end;
end;

if alpha ~= 0,
   c.dtf=dtf;
   c.th=th;
   c.ic1=ic1;
   c.ic2=ic2;
   c.metric=metric;
   c.alpha=alpha;
   c.p=p;
   c.pvalues = pvalues; % p-values associated to DTF/DC/iDTF
   c.patden = patdenr;
   c.patdf = patdfr;
   c.varass1 = varass1;
   c.varass2 = varass2;
   
   % Statistically significant DTF on frequency scale
   dtf_temp = ((abs(dtf)-abs(th)) > 0).*dtf + ((abs(dtf)-abs(th)) <= 0)*(-1);
   dtf_temp(dtf_temp < 0) = NaN; % Octave
   c.dtf_th = dtf_temp;   
else
   c.dtf=dtf;
   c.metric=metric;
   c.alpha=0;
   c.p=p;
   c.th=[];
   c.ic1=[];
   c.ic2=[];
   c.pvalues = [];
   c.patden = [];
   c.patdf = [];
   c.varass1 = [];
   c.varass2 = [];
end;

% Power spectra and coherence calculation
c.SS = ss_alg(A, pf, nFreqs);
c.coh = coh_alg(c.SS);

%==========================================================================
function gamma = bigautocorr(x, p)
%Autocorrelation. Data in rows. From order 0 to p-1.
%Output: nxn blocks of autocorr of lags i. (Nuttall Strand matrix)'''
[n, nd] = size(x);

gamma = zeros(n*p, n*p);
for i = 1:p
   for j = 1:p
      gamma(((i-1)*n+1):i*n, ((j-1)*n+1):j*n) = xlag(x, i-1)*(xlag(x,j-1).')/nd;
   end;
end;

%==========================================================================
function c= xlag(x,tlag)
if tlag == 0,
   c=x;
else
   c = zeros(size(x));
   c(:,(tlag+1):end) = x(:,1:(end-tlag));
end;

%==========================================================================
function d = fEig(L, G2)
%'''Returns the eigenvalues'''

%L = mat(cholesky(omega, lower=1))
D = L.'*G2*L;
%    d = eigh(D, eigvals_only=True)
%disp('fEig: eig or svd?')
d = svd(D);
d1=sort(d);
%
% the two biggest eigenvalues no matter which values (non negative by
% construction
%
d=d1(length(d)-1:length(d));

if (size(d) > 2),
   disp('more than two Chi-squares in the sum:')
end;

%==========================================================================
function c = fIij(i,j,n)
%'''Returns Iij of the formula'''
Iij = zeros(1,n^2);
Iij(n*(j-1)+i) = 1;
Iij = diag(Iij);
c = kron(eye(2), Iij);

%==========================================================================
function c = fIj(j,n)
%'''Returns Ij of the formula'''
Ij = zeros(1,n);
Ij(j) = 1;
Ij = diag(Ij);
Ij = kron(Ij,eye(n));
c = kron(eye(2), Ij);

%==========================================================================
function c = fIi(i,n)
%    '''Returns Ii of the formula'''
Ii = zeros(1,n);
Ii(i) = 1;
Ii = diag(Ii);
Ii = kron(eye(n), Ii);
c = kron(eye(2), Ii);

%==========================================================================
function d = fCa(f, p, n)
%'''Returns C* of the formula'''
C1 = cos(-2*pi*f*(1:p));
S1 = sin(-2*pi*f*(1:p));
C2 = [C1; S1];
d = kron(C2, eye(n^2));

%==========================================================================
function c = fdebig_de(n)
%'''Derivative of kron(I(2n), A) by A'''
%c = kron(TT(2*n, n), eye(n*2*n)) * kron(eye(n), kron(vec(eye(2*n)), eye(n)));
A=sparse(kron(TT(2*n, n), eye(n*2*n)));
B=sparse(kron(vec(eye(2*n)), eye(n)));
c = A * kron(eye(n), B);
c=sparse(c);

%==========================================================================
function c = fdebig_de_dtf(n)
%''' New \Theta_K for DTF asymptotics'''
%c = kron(kron(eye(2*n),TT(n, 2*n)), eye(n)));
%A=sparse(kron(eye(2*n), TT(n, 2*n)));
%c=sparse(kron(A, eye(n)));
A=sparse(kron(TT(n^2,2),eye(n^2)));
A1=sparse(kron(eye(2),A));
A=sparse(kron(TT(n^2,1),eye(n)));
%A=sparse(kron(TT(1,n^2),eye(n)));
A2=sparse(kron(eye(n),A));
A3=sparse(kron(eye(n^2),vec(eye(n))));
A4=sparse(A2*A3);
A5=sparse(kron(vec(eye(2)),A4));
c=sparse(A1*A5);
% linha extra para DCsparse
%c=sparse(c*diag(vec(eye(n))));

%==========================================================================
function c = vec(x)
%vec = lambda x: mat(x.ravel('F')).T
c=x(:);

%==========================================================================
function t = TT(a,b)
%''' TT(a,b)*vec(B) = vec(B.T), where B is (a x b).'''
t = zeros(a*b);
for i = 1:a,
   for j =1:b,
      t((i-1)*b+j,(j-1)*a+i) = 1;
   end;
end;
t = sparse(t);

%==========================================================================
function L = fChol(omega)
% Try Cholesky factorization
try,
   L = chol(omega)';
   % If there's a small negative eigenvalue, diagonalize
catch,
   %   disp('linalgerror, probably IP = 1.')
   [v,d] = eig(omega);
   L = zeros(size(v));
   for i =1:length(d),
      if d(i,i)<0,
         d(i,i)=eps;
      end;
      L(:,i) = v(:,i)*sqrt(d(i,i));
   end;
end;

%==========================================================================
function c = diagtom(a)
a=sparse(a');
c=sparse(diag(a(:)));

%==========================================================================
function c = mdiag(a)
%  diagonal matrix
c=diag(diag(a));

%==========================================================================
function d=Dup(n)
%     '''D*vech(A) = vec(A), with symmetric A'''
d = zeros(n*n, (n*(n+1))/2);
count = 1;
for j=1:n,
   for i =1:n,
      if i >= j,
         d((j-1)*n+i,count)=1;
         count = count+1;
      else
         d((j-1)*n+i,:)=d((i-1)*n+j,:);
      end;
   end;
end;

%==========================================================================
function hh=fdh_da(Af)
%    '''Derivative of vec(H) by vec(A), with H = A^-1 and complex A.'''
ha = pinv(Af);
h = -kron(ha.', ha);

h1 = [real(h) -imag(h)];
h2 = [imag(h) real(h)];  %h2 = cat(h.imag, h.real, 1)
hh = -[h1; h2];           %hh = cat(h1, h2, 0)

%==========================================================================
