function [MAR, Popt, Copt] = akaike_v3(Y, Pmin, Pmax)
%AkaikeIC3 modified
%graphic errors fixed for the case Pmin=Pmax

%% Input:
%   Y: Multivariate data [samples x channels]
%   Pmin: min optimal model order
%   Pmax: max optimal model order

%% Output:
%   MAR: MVAR coefficients for Popt
%   Popt: estimated optimal order
%   Copt: spectral matrix of residuals

%% Plot:
%   AIC as a function of p

%%% N rows = samples, M columns = time series
[N,M]=size(Y{1});

if (Pmin ~= round(Pmin) || Pmax ~= round(Pmax))
    error('Order must be integer.');
end

if nargin<2
    Pmin=1;
    Pmax=max(N,M)-1;
end

for p = Pmax:-1:Pmin
    [MARf,~,PE] = mvar_F_v3(Y, p, 2);
    C = PE(:,size(PE,2)+(1-M:0));
    Sigma =  eye(M).*(repmat(diag(C),1,M));
    d = det(Sigma);
    Aic(p) = ((N*reallog(d)) + (2*(M*M)*p));
    
    if p == Pmax
        Popt = p;
        MAR = MARf;
        Copt = C;
    elseif Aic(p) < Aic(Popt)
        Popt = p;
        MAR = MARf;
        Copt = C;
    end
end

%%% Graph Aic(p) in log scale and highlight the minimum of the function
if Pmin~=Pmax
    P=1:Pmax;
    B=[zeros(1,(Popt-1)),Aic(Popt),zeros(1,(Pmax-(Popt)))];
    K=Popt;
    semilogy(P,Aic,'-r',K,B,'go');
    grid on;
    xlabel('order (p)');
    ylabel('AIC(p)');
    title('Akaike Information Criterion Function');
    
    
    k=0;
    if Popt == Pmax
        for i = Pmax:-1:Pmin
            if (Aic(i)) < (Aic(i-1))
                k=k+1;
            else
                break
            end
        end
    end
    
    
    if k>0
        fprintf(1,'Attention: Estimated value for Popt could be not correct.\n');
        fprintf(1,'Aic function decreases from p = %d. \n', Pmax-k);
        fprintf(1,'Try to insert a bigger value for Pmax.\n');
    end
end