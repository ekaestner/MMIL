function [sf1,sf2] = compute_scalefact_hist2d(valvec1,valvec2,hcnt_target,binvals1,binvals2)
%function [sf1,sf2] = compute_scalefact_hist2d(valvec1,valvec2,hcnt_target,binvals1,binvals2)
%
% Created:           by Anders Dale
% Last Mod: 08/11/12 by Don Hagler
%

hc1 = sum(hcnt_target,2);
hc2 = sum(hcnt_target,1);
[sf1_1] = compute_scalefact_hist1d(valvec1,hc1,binvals1);
[sf2_1] = compute_scalefact_hist1d(valvec2,hc2,binvals2);

sf1 = sf1_1;
sf2 = sf2_1;

return

