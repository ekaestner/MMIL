function [pct1,pct2] = ts_calc_bs_pct(alpha)
%function [pct1,pct2] = ts_calc_bs_pct(alpha)
%
% Created:  11/05/13 by Don Hagler
% Last Mod: 11/05/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

pct1 = 100*alpha/2;
pct2 = 100-pct1;

return;
