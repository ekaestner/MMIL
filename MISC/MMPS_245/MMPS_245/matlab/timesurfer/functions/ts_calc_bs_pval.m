function pval = ts_calc_bs_pval(val,crit_vals_lo,crit_vals_hi,pval_range)
%function pval = ts_calc_bs_pval(val,crit_vals_lo,crit_vals_hi,pval_range)
%
% Created:  11/05/13 by Don Hagler
% Last Mod: 11/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pval = [];
if ~mmil_check_nargs(nargin,4), return; end;

ind_lo = find(val < crit_vals_lo);
ind_hi = find(val > crit_vals_hi);
if ~isempty(ind_lo)
  % get critical value of smallest pval
  crit_val = crit_vals_lo(min(ind_lo));
  if crit_val == crit_vals_lo(1)
    ind_lo = 1;
  else
    % get index of largest pval with same critical value
    ind_lo = max(find(crit_val == crit_vals_lo));
  end;
  % get pval
  pval = pval_range(ind_lo);
elseif ~isempty(ind_hi)
  % get critical value of smallest pval
  crit_val = crit_vals_hi(min(ind_hi));
  if crit_val == crit_vals_hi(1)
    ind_hi = 1;
  else
    % get index of largest pval with same critical value
    ind_hi = max(find(crit_val == crit_vals_hi));
  end;
  % get pval
  pval = pval_range(ind_hi);
else
  pval = 1;
end;

return;

