function pval = ts_calc_bs_null_pval(vals,pval_range)
%function pval = ts_calc_bs_null_pval(vals,pval_range)
%
% Created:  11/05/13 by Don Hagler
% Last Mod: 11/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pval = [];
if ~mmil_check_nargs(nargin,2), return; end;

% reshape to row vector
vals = mmil_rowvec(vals);

% calculate percentile values
[pct1,pct2] = ts_calc_bs_pct(pval_range);

% calculate critical values
crit_vals_lo = prctile(vals,pct1,2);
crit_vals_hi = prctile(vals,pct2,2);

% calculate p-value for 0
pval = ts_calc_bs_pval(0,crit_vals_lo,crit_vals_hi,pval_range);

return;

