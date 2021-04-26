function c = dtf_analysis(u,A,pf,nFreqs,metric,alpha)

%==========================================================================
%            DTF, threshold and confidence interval calculation.
%==========================================================================

c = asymp_dtf(u,A,pf,nFreqs,metric,alpha);

% Power spectra and coherence calculation
c.SS = ss_alg(A, pf, nFreqs);
c.coh = coh_alg(c.SS);

% % Statistically significant DTF on frequency scale
% dtf_temp = ((abs(c.dtf)-c.th) > 0).*c.dtf + ((abs(c.dtf)-c.th) <= 0)*(-1);
% dtf_temp(ind2sub(size(dtf_temp),find(dtf_temp == -1))) = NaN;
% c.dtf_th = dtf_temp;

% disp('Find DTF results in "c" structured array.')
% disp(' ');