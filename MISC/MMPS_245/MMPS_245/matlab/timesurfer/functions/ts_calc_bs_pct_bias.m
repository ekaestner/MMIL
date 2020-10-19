function [pct1,pct2] = ts_calc_bs_pct_bias(alpha,z0,acc)
%function [pct1,pct2] = ts_calc_bs_pct_bias(alpha,z0,acc)
%
% Created:  11/05/13 by Don Hagler
% Last Mod: 11/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pct1 = []; pct2 = [];
if ~mmil_check_nargs(nargin,1), return; end;

z_alpha1 = norminv(alpha/2);
z_alpha2 = -z_alpha1;
pct1 = 100*normcdf(z0+(z0+z_alpha1)./(1-acc.*(z0+z_alpha1)));
pct2 = 100*normcdf(z0+(z0+z_alpha2)./(1-acc.*(z0+z_alpha2)));

return;

