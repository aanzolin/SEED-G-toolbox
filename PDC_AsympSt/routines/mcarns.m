
function [pf,A,pb,B,ef,eb,ISTAT]=mcarns(u,IP)
% MVAR Nuttall-Strand algorithm
%
%[PF,A,PB,B,EF,EB,ISTAT]=MCARNS(u,IP)
%   Calculate the coeficients of multi-channel auto-regressive matrix using
%   Nuttall-Strand algorithm (a generalization of single channel harmonic
%                             method).
%
%   Input parameters:
%     IP     - Ordem of autoregressive model (integer)
%     u      - Complex matrix with NUMCHS channels of sample data
%
%   Output parameters:
%     PF     - Covariance matrix of NUMCHS x NUMCHS of linear forward
%              prediction error
%     A      - Complex array of forward linear prediction matrix
%              coefficients
%     PB     - Complex backward linear prediction error covariance array
%     B      - Complex array of backward linear prediction matrix
%              coefficients

%%   This implementation is a FORTRAN code translation from 
%   Appendix 15.B page 424 of Marple Jr.(1987). (KS 1998)
%
%   Marple Jr., SL. Digital Spectral Analysis with Application.
%   Prentice-Hall, Englewood-Cliffs, 1987. 

if nargin ~= 2, error('MCARNS requires 2 input parameters.'); end;

[lx,cx]=size(u);

if lx > cx, error('Input matrix is probably transposed.'), end;
NUMCHS=lx;          % Number of channels.
MAXORDER=200;       % Maximum order of AR model allowed for calculation.
N=max(size(u));     % N - Number of samples per channel.

%   ** Initialization **
ISTAT=0;
if (IP > MAXORDER),
   error('IP > 200.');
end;

ef=u;        % Eq. (15.91)
eb=u;        % Eq. (15.91)
%Tef=ef;     % Temporary variable ( never used );
pf=u*u';     % Eq. (15.90)
pb=pf;       % Eq. (15.90)

M=0;
%   ** Main Loop **
while 1,
   %  Update estimated covariance errors               Eq. (15.89)
   pfhat=ef(:,M+2:N)*ef(:,M+2:N)';
   pbhat=eb(:,M+1:N-1)*eb(:,M+1:N-1)';
   pfbhat=ef(:,M+2:N)*eb(:,M+1:N-1)';

   M=M+1;

   %  Calculate estimated partial correlation matrix - Eq. (15.98)
   %             (Nuttall-Strand algorithm only)
   RHO=lyap(pfhat*inv(pf),inv(pb)*pbhat,-2*pfbhat);

   %  Update forward and backward reflection coeficients
   %  Eqs. (15.73),(15.74),(15.78) (algoritmo de Nuttall-Strand)
   AM=-RHO*inv(pb);
   BM=-RHO'*inv(pf);
   %   if exist('A')~=2, A=[]; end;
   A(:,:,M)=AM;
   %   if exist('B')~=2, B=[]; end;
   B(:,:,M)=BM;
   %
   %  Update forward and backward covariance error  - Eqs. (15.75),(15.76)
   pf=pf-AM*BM*pf;
   pb=pb-BM*AM*pb;

   %  Update forward and backward predictor coeficients - Eqs.(15.84),(15.85)
   if M ~= 1,
      for K=1:M-1
         temp1=A(:,:,K);
         A(:,:,K)=A(:,:,K)+AM*B(:,:,M-K);
         B(:,:,M-K)=B(:,:,M-K)+BM*temp1;
      end;
   end;
   %  Update residues
   Tef=ef;
   ef(1:NUMCHS,N:-1:M+1)=ef(:,N:-1:M+1)+AM*eb(:,N-1:-1:M);
   eb(1:NUMCHS,N:-1:M+1)=eb(:,N-1:-1:M)+BM*Tef(:,N:-1:M+1);

   %  Verify if model order is adequate
   if M == IP, A=-A; B=-B; break, end;
end;
