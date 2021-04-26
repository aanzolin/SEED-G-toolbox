% Estimate MVAR
%
%% Syntax
%[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion)
%
%% Arguments
%
% input:
%        u     - data rows
%        maxIP - externaly defined maximum IP
%        alg   - estimation algorithm (=1 Nutall-Strand), (=2 - mlsm) ,
%                              (=3 - Vieira Morf), (=4 - QR artfit)
%        criterion for order choice - 0: MDL (not implemented)
%                                     1: AIC; 2: Hanna-Quinn, 3 Schwarz,
%                                     4: FPE, 5 - fixed order in maxIP
% output: See Data Structure Data Base
%
%% Description
%
%
%% References
%
% [1]
%
% [2]
%
%% Copyright
%
%
%%

function [IP,pf,A,pb,B,ef,eb,vaic,Vaicv]=mvar(u,maxIP,alg,criterion)

[nSegLength,nChannels] = size(u');
if nargin<3, alg=1; criterion = 1; end % defaults
if nargin<4, criterion = 1; end        % defaults
% blancos = '                                  ';
% switch alg      %  Formula from Marple Jr.
%    case 1,
%       disp([blancos 'using Nuttall-Strand algorithm.'])
%    case 2,
%       disp([blancos 'using least squares estimator.'])
%    case 3,
%       disp([blancos 'using Vieira-Morf algorithm.'])
%    case 4,
%       disp([blancos 'using ARfit estimator.'])
%    otherwise,
%       error('Unknown alg option.')
% end


if criterion==5 % Fixed order in maxIP
   if maxIP==0
      pf = u*u';     % Eq. (15.90) Equation numbering refers to Marple Jr.
      pb = pf;       % Eq. (15.90)
      ef = u;
      eb = u;
      npf = size(pf);
      A = zeros(npf,npf,0);
      B = A;
      vaic = det(pf);
      Vaicv = det(pf);
      IP = 0;
      return
   end
   IP = maxIP;
   switch alg      %  Formula from Marple Jr.
      case 1,
         [pf,A,pb,B,ef,eb] = mcarns(u,IP);
         pf = pf(:,:)/nSegLength; %Nuttall-Strand needs scaling.
      case 2,
         [pf,A,ef] = cmlsm(u,IP);
         B  = [ ]; eb = [ ]; pb = [ ];
         pf = pf(:,:)/nSegLength;
      case 3,
         [pf,A,pb,B,ef,eb] = mcarvm(u,IP);
         pf=pf(:,:)/nSegLength; %Vieira-Morf needs scaling.
      case 4,
         [pf,A,ef] = arfitcaps(u,IP);
         B  = [ ]; eb = [ ]; pb = [ ];
         % Arfit does not require scaling. pf=pf;
   end

   vaic=length(u)*log(det(pf))+2*(nChannels*nChannels)*IP;
   Vaicv=vaic;
   return
end
%
vaicv=0;
if nargin<2
   maxOrder=30;
   disp(['maxOrder limited to ' int2str(maxOrder)])
   UpperboundOrder = round(3*sqrt(nSegLength)/nChannels); % Marple Jr. p. 409
   % Suggested by Nuttall, 1976.
   UpperboundOrder = min([maxOrder UpperboundOrder]);
else
   maxOrder=maxIP;
   UpperboundOrder=maxIP;
   disp(['maxOrder limited to ' int2str(maxOrder)])
end

IP=1;
Vaicv=zeros(maxOrder+1,1);

T = length(u);
K = nChannels;


while IP <= UpperboundOrder,
   m = IP;
   switch alg
      case 1, 
         [npf,na,npb,nb,nef,neb] = mcarns(u,IP);
      case 2,
         [npf,na,nef] = cmlsm(u,IP);
      case 3,
         [npf,na,npb,nb,nef,neb] = mcarvm(u,IP);
      case 4,
         [npf,na,nef] = arfitcaps(u,IP);
   end

   switch criterion
      case 1,  % Akaike's Informaion Criterin (AIC)
         vaic = length(u)*log(det(npf)) ...
                + 2*(nChannels*nChannels)*IP;                     %(4.3.2)(Lutkepohl'05)
      case 2,  % Hannan-Quinn (HQ)
         vaic = length(u)*log(det(npf)) ...
                + 2*log(log(length(u)))*(nChannels*nChannels)*IP; %(4.3.8)(Lutkepohl'05) 
      case 3,  % Schwarz (SC) (Schwarz, 1978)
         vaic = length(u)*log(det(npf)) ...
                + log(length(u))*(nChannels*nChannels)*IP;        %(4.3.9)(Lutkepohl'05)
      case 4,  % FPE - Final prediction error (Akaike, 1970)
          vaic = log(det(npf)*((length(u)+nChannels*IP+1) ...
                 /(length(u)-nChannels*IP-1))^nChannels);         %(4.3.1)(Lutkepohl'05)
      otherwise
         %nop
   end

   Vaicv(IP+1) = vaic;
   fprintf('IP=%0.f  vaic=%f\n',IP,vaic);
   if (vaic>vaicv) && (IP~=1),
      vaic = vaicv;
      break;
   end; % Akaike condition
   vaicv = vaic;
   pf = npf;
   A  = na;
   ef = nef;
   if alg==1 || alg == 3,
      B  = nb;
      eb = neb;
      pb = npb;
   else % review status for backward prediction in clmsm
      B  = [];
      eb = [];
      pb = [];
   end
   IP = IP+1;
end;

disp(' ')
IP = IP-1;
vaic = vaicv;
Vaicv = Vaicv(2:IP+1);

switch alg
   case {1,3}, %Nuttall-Strand and Vieira-Morf need scaling.
      pf = pf(:,:)/nSegLength;
   case 4,
      % nop 
      % pf = pf; % arfit does not need scaling.
   otherwise  %
      pf = pf(:,:)/nSegLength;
end;
