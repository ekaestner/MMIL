function pval = ts_calc_bs_pval_bias(val,bstat,jstat,stat,pval_range)
%function pval = ts_calc_bs_pval_bias(val,bstat,jstat,stat,pval_range)
%
% Created:  05/24/16 by Don Hagler
% Last Mod: 05/24/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pval = [];
if ~mmil_check_nargs(nargin,4), return; end;

z0 = ts_calc_bs_bias(bstat,stat);
acc = ts_calc_bs_acceleration(jstat);
[pct1,pct2] = ts_calc_bs_pct_bias(pval_range,z0,acc);
crit_vals_lo = prctile(bstat,pct1);
crit_vals_hi = prctile(bstat,pct2);
pval = ts_calc_bs_pval(val,crit_vals_lo,crit_vals_hi,pval_range);

return;

