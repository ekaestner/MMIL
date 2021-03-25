function z0 = ts_calc_bs_bias(bstat,stat)
%function z0 = ts_calc_bs_bias(bstat,stat)
%
% Created:  11/05/13 by Don Hagler
% Last Mod: 11/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z0 = [];
if ~mmil_check_nargs(nargin,2), return; end;

% compute bias-correction constant z0
z0 = norminv(mean(bsxfun(@lt,bstat,stat),2) + mean(bsxfun(@eq,bstat,stat),2)/2);

return;

