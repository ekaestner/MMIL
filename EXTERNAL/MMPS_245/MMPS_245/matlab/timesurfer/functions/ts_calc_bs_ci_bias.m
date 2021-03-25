function ci = ts_calc_bs_ci_bias(bstat,jstat,stat,alpha)
%function ci = ts_calc_bs_ci_bias(bstat,jstat,stat,alpha)
%
% Created:  11/05/13 by Don Hagler
% Last Mod: 11/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ci = [];
if ~mmil_check_nargs(nargin,4), return; end;

% compute bias-correction constant z0
z0 = ts_calc_bs_bias(bstat,stat);

% compute acceleration (see DiCiccio and Efron (1996))
acc = ts_calc_bs_acceleration(jstat);

% calculate bias-corrected percentile values
[pct1,pct2] = ts_calc_bs_pct_bias(alpha,z0,acc);

% calculate confidence intervals
n = numel(stat);
ci = zeros(n,2);
for i=1:n
  tmp = bstat(i,:);
  ci(i,1) = prctile(tmp,pct1(i),2);
  ci(i,2) = prctile(tmp,pct2(i),2);
end;
if n>1, ci = reshape(ci,[n,1,2]); end;

return

