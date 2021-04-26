function c=asymp_pdc(x,A,pf,nFreqs,metric,alpha)
%Compute |PDC|^2 connectivity measure given by metric "option" from series j-->i.
%
%function c=asymp_pdc(x,A,pf,nFreqs,metric,alpha)
%
% input: x      - data
%        A      - AR estimate matrix by MVAR
%        pf     - covariance matrix provided by MVAR
%        nFreqs - number of point in [0,fs/2] frequency scale
%        metric - 'euc'  - Euclidean   ==> original PDC
%                 'diag' - diagonal    ==> gPDC (generalized )
%                 'info' - information ==> iPDC
%        alpha  - significance level
%                if alpha = zero, do not calculate statistics
% output:  c.pdc       - |PDC|^2 estimates
%          c.pvalues   - p-values associated to pdc estimates. 
%          c.th        - Threshold value with (1-avalue) significance level.
%          c.ic1,c.ic2 -  confidence interval
%          c.metric    - metric used for PDC calculation
%          c.alpha     - significance level
%          c.p         - VAR model order
%          c.patdenr   - 
%          c.patdfr    - degree of freedom 
%          c.SS        - power spectra
%          c.coh       - spectral coherence

%Corrected 7/25/2011 to match the frequency range with plotting
%routine, f = 0 was include in the frequency for loop:
%                                for ff = 1:nFreqs,
%                                   f = (ff-1)/(2*nFreqs); % 
%                                        ^?^^
% Optimization 7/1/2015 


if nargin<6,
   error('ASYMP_PDC requires six input arguments.')
end
[m,n]=size(x);
if m > n,
   x=x.';
end;
np = length(x);
[nChannels,dummy,p] = size(A);
Af = A_to_f(A, nFreqs);

% Variables initialization
pdc = zeros(nChannels,nChannels,nFreqs);

if alpha ~= 0,
   th   = zeros(nChannels,nChannels,nFreqs);
   ic1 = zeros(nChannels,nChannels,nFreqs);
   ic2 = zeros(nChannels,nChannels,nFreqs);
   varass1 = zeros(nChannels,nChannels,nFreqs);
   varass2 = zeros(nChannels,nChannels,nFreqs);
   patdfr = zeros(nChannels,nChannels,nFreqs);
   patdenr = zeros(nChannels,nChannels,nFreqs);
   pvalues = zeros(nChannels,nChannels,nFreqs);
   switch lower(metric)
      case {'euc'}
         disp('* Original PDC and asymptotic statistics')
      case {'diag'}
         disp('* Generalized PDC and asymptotic statistics')
      case {'info'}
         disp('* Information PDC and asymptotic statistics')
      otherwise
         error('Unknown metric.')
   end;
else
   switch lower(metric)
      case {'euc'}
         disp('* Original PDC estimation')
      case {'diag'}
         disp('* Generalized PDC estimation')
      case {'info'}
         disp('* Information PDC estimation')
      otherwise
         error('Unknown metric.')
   end;
end;
%disp('======================================================================');

switch lower(metric)
   case {'euc'}           % for original PDC
      % nop()
      
      %Iije = Iij;
      %Ije = Ij;
      
      %pu =  eye(nChannels);
      dpdc_dev = zeros(1,(nChannels*(nChannels+1))/2);
      
   case {'diag'}          % for gPDC
      evar_d = mdiag(pf);
      evar_d_big = kron(eye(2*nChannels), evar_d);
      pinv_evar_d_big = pinv(evar_d_big);
      
      % if (i == 1) && (j == 1) && (ff == 1),
%      evar_d = mdiag(pf);
%      evar_d_big = kron(eye(2*nChannels), evar_d);
%      inv_ed = pinv(evar_d_big);
      
      %'derivative of vec(Ed-1) by vecE'
      de_deh = Dup(nChannels);
      debig_de = fdebig_de(nChannels);
      
      dedinv_dev = diagtom(vec(-pinv_evar_d_big*pinv_evar_d_big));
      dedinv_deh = dedinv_dev*debig_de*de_deh;
      
   case {'info'}          % for iPDC
      evar_d = mdiag(pf);
      evar_d_big = kron(eye(2*nChannels), evar_d);
      pinv_evar_d_big = pinv(evar_d_big);
      
      evar_big = kron(eye(2*nChannels), pf);
      pinv_evar_big = pinv(evar_big);

%      inv_ed = pinv(evar_d_big);
%      evar_big = kron(eye(2*nChannels), pf);
      inv_e = sparse(pinv(evar_big));
      
      %'derivative of vec(Ed-1) by vecE'
      de_deh = Dup(nChannels);
      debig_de = fdebig_de(nChannels);
      
      dedinv_devd = sparse(diagtom(vec(-pinv_evar_d_big*pinv_evar_d_big)));
      dedinv_dehd = sparse(dedinv_devd*debig_de*de_deh);
      
      dedinv_dev = sparse(-kron(inv_e.', inv_e));
      dedinv_deh = sparse(dedinv_dev*debig_de*de_deh);

   otherwise
      error('Unknown metric.')
end;

gamma = bigautocorr(x, p);
omega = kron(inv(gamma), pf);
omega_evar = 2*pinv(Dup(nChannels))*kron(pf, pf)*pinv(Dup(nChannels)).';

%if (is_octave)
   icdf_norm_alpha = norminv(1-alpha/2.0,0,1);
%else
%   icdf_norm_alpha = icdf('norm',1-alpha/2.0,0,1);
%end

for ff = 1:nFreqs,
   f = (ff-1)/(2*nFreqs); %Corrected 7/25/2011, f starting at 0 rad/s.
   Ca = fCa(f, p, nChannels);
   
   a = Af(ff,:,:); a = a(:);    %Equivalent to a = vec(Af[ff, :, :])
   a = [real(a); imag(a)];      %a = cat(a.real, a.imag, 0)
   
   omega2 = Ca*omega*Ca.';
   L = fChol(omega2);

   for j = 1:nChannels,
      Ij = fIj(j,nChannels);
      switch lower(metric)
         case {'euc'}           % for PDC
            Ije = Ij;
         case {'diag'}          % for gPDC
            Ije = Ij*pinv_evar_d_big;
         case {'info'}          % for iPDC
            Ije = Ij*pinv_evar_big*Ij;
         otherwise
            error('Unknown metric.')
      end;
      
      for i = 1:nChannels,
         Iij = fIij(i,j,nChannels);
         %For metric = diag or info case, include evar in the expression'
         switch lower(metric)
            case {'euc'}           % for PDC
               Iije = Iij;
            case {'diag'}          % for gPDC
               Iije = Iij*pinv_evar_d_big;
            case {'info'}          % for iPDC
               Iije = Iij*pinv_evar_d_big;
            otherwise
               error('Unknown metric.')
         end;

         num = a.'*Iije*a;
         den = a.'*Ije*a;
         pdc(i,j,ff) = num/den;
         
         % If alpha == 0, no statistical calculation for faster computation.
         if alpha ~= 0,
            %'Add evar derivation'
            switch lower(metric)
               case {'euc'}
%                  dpdc_dev = zeros(1,(nChannels*(nChannels+1))/2);
               case {'diag'}
                  %'derivada do num por vecE'
                  dnum_dev = kron((Iij*a).', a.')*dedinv_deh;
                  %'derivada do den por vecE'
                  dden_dev = kron((Ij*a).', a.')*dedinv_deh;
                  dpdc_dev = (den*dnum_dev - num*dden_dev)/(den^2);
               case {'info'}
                  %'derivada do num por vecE'
                  dnum_dev = kron((Iij*a).', a.')*dedinv_dehd;
                  %'derivada do den por vecE'
                  dden_dev = kron((Ij*a).', a.'*Ij)*dedinv_deh;
                  dpdc_dev = (den*dnum_dev - num*dden_dev)/(den^2);
               otherwise
                  error('Unknown metric.')
            end;

            G1a = 2 * a.' * Iije/den - 2*num*a.'*Ije/(den^2);
            G1 = -G1a * Ca;
            varalpha = G1 * omega * G1.';
            varevar = dpdc_dev * omega_evar * dpdc_dev.';
            varass1(i,j,ff) = (varalpha + varevar)/np;
            
            ic1(i,j,ff) = pdc(i,j,ff) ...
                    - sqrt(varass1(i,j,ff))*icdf_norm_alpha;
            ic2(i,j,ff) = pdc(i,j,ff) ...
                    + sqrt(varass1(i,j,ff))*icdf_norm_alpha;

            G2a = Iije/den;
                        
            d = fEig(real(L), real(G2a)); % real() 28May2013            

            patdf = (sum(d).^2)./sum(d.^2);
            patden = sum(d)./sum(d.^2);

%             if (is_octave)
                th(i,j,ff) = chi2inv((1-alpha), patdf)/(patden*np);
                pvalues(i,j,ff) = 1 - chi2cdf(pdc(i,j,ff)*patden*np, patdf);
%             else
%                 th(i,j,ff) = icdf('chi2',(1-alpha), patdf)/(patden*np);
%                 pvalues(i,j,ff) = 1 - cdf('chi2', pdc(i,j,ff)*patden*np, patdf);
%                 
%                          
%                achi2inv = chi2inv((1-alpha), patdf)/(patden*np);
%                achi2cdf = 1 - chi2cdf(pdc(i,j,ff)*patden*np, patdf);
%                
%                th(i,j,ff) - achi2inv
%                pvalues(i,j,ff) - achi2cdf
%       
%                 
%             end
           
            varass2(i,j,ff) = patdf/(patden*np).^2;
            patdfr(i,j,ff) = patdf;
            patdenr(i,j,ff) = patden;
         else % alpha == 0, do not compute asymptotics
            % nop;
         end;
      end;
   end;
end;

if alpha ~= 0,
   c.pdc=pdc;
   c.th=th;
   c.ic1=ic1;
   c.ic2=ic2;
   c.metric=metric;
   c.alpha=alpha;
   c.p=p;
   c.pvalues = pvalues; % p-values associated to PDC/gPDC/iPDC
   c.patden = patdenr;
   c.patdf = patdfr;
   c.varass1 = varass1;
   c.varass2 = varass2;   
%  Statistically significant PDC on frequency scale
   pdc_temp = ((abs(pdc)-abs(th)) > 0).*pdc + ((abs(pdc)-abs(th)) <= 0)*(-1);
   pdc_temp(pdc_temp < 0) = NaN;
   c.pdc_th = pdc_temp;
else
   c.pdc=pdc;
   c.metric=metric;
   c.alpha=0;
   c.p=p;
   c.pvalues = [];   
   c.th=[];
   c.ic1=[];
   c.ic2=[];
   c.patdenr = [];
   c.patdfr = [];
   c.varass1 = [];
   c.varass2 = [];
end;

%%
% Power spectra and coherence calculation
c.SS = ss_alg(A, pf, nFreqs);
c.coh = coh_alg(c.SS);
%

%==========================================================================
function gamma = bigautocorr(x, p)
%Autocorrelation. Data in row-wise orientation. From order 0 to p-1.
%Output: n x n blocks of autocorr of lags i. (Nuttall Strand matrix)'''
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

c = sparse(c);  % SPARSED

%==========================================================================
function c=  fIj(j,n)
%'''Returns Ij of the formula'''
Ij = zeros(1,n);
Ij(j) = 1;
Ij = diag(Ij);
Ij = kron(Ij,eye(n));
c = kron(eye(2), Ij);

c = sparse(c);  % SPARSED

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
try
   L = chol(omega)';
   % If there's a small negative eigenvalue, diagonalize
catch
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
b = sparse(a');
c = sparse(diag(b(:)));

%==========================================================================
function c = mdiag(a)
%  diagonal matrix
c = diag(diag(a)); % UNSPARSED

%==========================================================================
function d=Dup(n)
%     '''D*vech(A) = vec(A), with symmetric A'''
d = zeros(n*n, (n*(n+1))/2);
count = 1;
for j=1:n,
   for i =1:n,
      if i >= j,
         d((j-1)*n+i,count) = 1;
         count = count+1;
      else
         d((j-1)*n+i,:) = d((i-1)*n+j,:);
      end;
   end;
end;

%==========================================================================