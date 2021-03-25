function mr_parms = mmil_serinfo_mrparms(serinfo)
%function mr_parms = mmil_serinfo_mrparms(serinfo)
%
% Required Input:
%   serinfo: series info struct created by mmil_classify_dicoms
%     must be a single struct, not an array
%
% Output:
%   mr_parms: vector containing values for [FA,TR,TE,TI,FOV]
%
% Created:  09/11/12 by Don Hagler
% Last Mod: 09/11/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
mr_parms = [];

FA = serinfo.FlipAngle;
TR = serinfo.RepetitionTime;
te_list = sort(unique(serinfo.EchoTimes));
TE = te_list(1);
TI = serinfo.InversionTime;
FOV = serinfo.readoutFOV;
mr_parms = [TR FA TE TI FOV];

