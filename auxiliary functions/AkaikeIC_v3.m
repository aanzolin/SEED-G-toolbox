function [MAR, Popt, Copt] = AkaikeIC_v3(Y, Pmin, Pmax)
%AkaikeIC3 modified
%graphic errors fixed for the case Pmin=Pmax

%INPUT:
% Y= Dati multivariati [campioni canali]
% Pmin=ordine minimo da testare
% Pmax=ordine massimo da testare

%OUTPUT:
% MAR = Parametri MVAR relativi a P=Popt
% Popt = Ordine ottimo del modello
% Copt = Matrice spettrale dei residui relativa all'ordine ottimo...
% ...CON APPROSSIMAZIONE: solo elementi della diagonale.

%PLOT:
% AIC as a function of p
% N rows = samples, M columns = time series

[N,M]=size(Y{1});

if (Pmin ~= round(Pmin) || Pmax ~= round(Pmax))
    error('Order must be integer.');
end

if nargin<2
    Pmin=1; % meno LAG di così non si può avere
    Pmax=max(N,M)-1; % più LAG di così non si puo' avere
end


% Il ciclo for procede a ritroso da Pmax a Pmin: in questo modo
% le matrici MARf ottenute come risultato da arfit avranno dimensioni
% DECRESCENTI
for p = Pmax:-1:Pmin
    % Si calcolano i parametri multivariati e le varianze dei residui
    % relativi ad ogni p considerato
    [MARf,~,PE] = mvar_F_v3(Y, p, 2);
    C = PE(:,size(PE,2)+(1-M:0));
    
    % Si ricavano le varie sigma
    Sigma =  eye(M).*(repmat(diag(C),1,M));
    
    % d contiene i determinanti delle sigma
    d = det(Sigma);
    
    % Si calcola AIC(p)
    Aic(p) = ((N*reallog(d)) + (2*(M*M)*p));
    
    % Volta per volta si memorizza il valore minimo di Aic ed
    % i parametri MAR ad esso relativi
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


% Si traccia il grafico di Aic(p) in scala semilogaritmica,
% evidenziando il punto di minimo e la relativa ascissa p = Popt

if Pmin~=Pmax
    P=1:Pmax;
    B=[zeros(1,(Popt-1)),Aic(Popt),zeros(1,(Pmax-(Popt)))];
    K=Popt;
    semilogy(P,Aic,'-r',K,B,'go');
    grid on;
    %xlim ([Pmin Pmax]);
    xlabel('order (p)');
    ylabel('AIC(p)');
    title('Akaike Information Criterion Function');
    
    
    % Capire se la funzione è decrescente
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
    
    
    % Se la funzione è decrescente si avvisa l'utente
    if k>0
        fprintf(1,'Attention: Estimated value for Popt could be not correct.\n');
        fprintf(1,'Aic function decreases from p = %d. \n', Pmax-k);
        fprintf(1,'Try to insert a bigger value for Pmax.\n');
    end
end