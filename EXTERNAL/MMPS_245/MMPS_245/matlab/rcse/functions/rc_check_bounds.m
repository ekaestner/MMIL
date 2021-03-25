function vals = rc_check_bounds(vals,valrange)
%function vals = rc_check_bounds(vals,valrange)
%
% Purpose: clip vals to valrange
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/04/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

vals = max(min(vals,valrange(end)),valrange(1));

