function pval = ts_calc_bs_val_pval(val,vals,pval_range)
%function pval = ts_calc_bs_val_pval(val,vals,pval_range)
%
% Created:  05/24/16 by Don Hagler
% Last Mod: 05/24/16 by Don Hagler
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
pval = ts_calc_bs_pval(val,crit_vals_lo,crit_vals_hi,pval_range);

return;

