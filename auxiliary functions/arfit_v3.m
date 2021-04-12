function [MAR, C, popt]=arfit_v3(Y, pmin, pmax, selector, no_const)

%INPUT:
%     Y = Multivariate data [samp channels]
%     pmin = min optimal order to test
%     pmax = max optimal order to test
%     selector = optimal order selection criteria: 'sbc'= Schwartz's Bayesian Criterion, 'fpe'= Akaike Final Prediction Error, 'zero'
%     no_const = 'zero' (assign to zero the intercept).

%OUTPUT:
%     MAR = MVAR parameters
%     C = Residual Power 
%     popt = optimal order

% ARFIT estimates multivariate autoregressive parameters
%   using MVAR with the Nuttall-Strand method [1,2].
% ARFIT is included for combatibility reasons to ARFIT [3]
%  [w, A, C, sbc, fpe] = arfit2(v, pmin, pmax, selector, no_const)
% see also: ARFIT, MVAR
% REFERENCES:
%  [1] M.S. Kay "Modern Spectral Estimation" Prentice Hall, 1988.
%  [2] S.L. Marple "Digital Spectral Analysis with Applications" Prentice Hall, 1987.
%  [3] T. Schneider and A. Neumaier, A. 2001.
%	Algorithm 808: ARFIT-a Matlab package for the estimation of parameters and eigenmodes
%	of multivariate autoregressive models. ACM-Transactions on Mathematical Software. 27, (Mar.), 58-65.

%       Revision: 1.9 $
%       Id: arfit2.m,v 1.9 2004/03/26 17:23:06 schloegl Exp $
%       Copyright (C) 1996-2004 by Alois Schloegl  <a.schloegl@ieee.org>

if (pmin ~= round(pmin) || pmax ~= round(pmax))
    error('Order must be integer.');
end
if (pmax < pmin)
    error('PMAX must be greater than or equal to PMIN.')
end

%  Set defaults and check for optional arguments

if (nargin == 3)              	% no optional arguments => set default values
    mcor = 1;                   % fit intercept vector
    selector  = 'sbc';	        % use SBC as order selection criterion
elseif (nargin == 4)          	% one optional argument
    if strcmp(selector, 'zero')
        mcor = 0;               % no intercept vector to be fitted
        selector = 'sbc';	    % default order selection
    else
        mcor = 1;	            % fit intercept vector
    end
elseif (nargin == 5)		    % two optional arguments
    if strcmp(no_const, 'zero')
        mcor = 0;               % no intercept vector to be fitted
    else
        error(['Bad argument. Usage: ', '[w,A,C,SBC,FPE,th]=AR(v,pmin,pmax,SELECTOR,''zero'')'])
    end
end

% Implementation of the MVAR estimation

[N,M]=size(Y{1});

if mcor
    [m,N] = sumskipnan(Y,2);               % calculate mean
    m = m./N;
    Y = Y - repmat(m,size(Y)./size(m));    % remove mean
end

[MAR,~,PE] = mvar_F_v3(Y, pmax, 2);         % estimate MVAR(pmax) model

ne = N-mcor-(pmin:pmax);
for p=pmin:pmax
    % Get downdated logarithm of determinant
    logdp(p-pmin+1) = log(det(PE(:,p+(1:M))*(N-p-mcor)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(selector,'sbc')

    % Schwarz's Bayesian Criterion
    sbc = logdp/M - log(ne) .* (1-(M*(pmin:pmax)+mcor)./ne);


    [val, iopt]  = min(sbc);

    % select order of model
    popt = pmin + iopt-1; % estimated optimum order

    if pmin~=pmax
        ax=[pmin:pmax];
        plot(ax,sbc);hold on;plot(popt,sbc(find(ax==popt)),'or');grid on;
    end

    if popt<pmax
        pmax=popt;
        [MAR, ~, PE] = mvar_F_v3(Y, pmax, 2);
    end

    C = PE(:,size(PE,2)+(1-M:0));

    if mcor
        I = eye(M);
        for k = 1:pmax
            I = I - MAR(:,k*M+(1-M:0));
        end
        w = I*m;
    else
        w = zeros(M,1);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(selector,'fpe')

    % logarithm of Akaike's Final Prediction Error
    fpe = logdp/M - log(ne.*(ne-M*(pmin:pmax)-mcor)./(ne+M*(pmin:pmax)+mcor));

    [val, iopt]  = min(fpe);

    % select order of model
    popt = pmin + iopt-1; % estimated optimum order

    if pmin~=pmax
        ax=[pmin:pmax];
        plot(ax,fpe);hold on;plot(popt,fpe(find(ax==popt)),'or');grid on;
    end

    if popt<pmax
        pmax=popt;
        [MAR, ~, PE] = mvar_F_v3(Y, pmax, 2);
    end

    C = PE(:,size(PE,2)+(1-M:0));

    if mcor
        I = eye(M);
        for k = 1:pmax
            I = I - MAR(:,k*M+(1-M:0));
        end
        w = I*m;
    else
        w = zeros(M,1);
    end

end

if strcmp(selector,'AIK')
    [MAR, popt, C] = AkaikeIC_v3(Y, pmin, pmax);
end
