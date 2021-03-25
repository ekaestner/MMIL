function ci = ts_calc_bs_ci(bstat,alpha)
%function ci = ts_calc_bs_ci(bstat,alpha)
%
% Created:  11/05/13 by Don Hagler
% Last Mod: 11/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ci = [];
if ~mmil_check_nargs(nargin,2), return; end;

% calculate percentile values
[pct1,pct2] = ts_calc_bs_pct(alpha);

% calculate confidence intervals
n = size(bstat,1);
ci = zeros(n,2);
ci(:,1) = prctile(bstat,pct1,2);
ci(:,2) = prctile(bstat,pct2,2);
if n>1, ci = reshape(ci,[n,1,2]); end;

return;

