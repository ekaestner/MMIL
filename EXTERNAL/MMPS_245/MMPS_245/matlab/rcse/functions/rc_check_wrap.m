function vals = rc_check_wrap(vals,valrange)
%function vals = rc_check_wrap(vals,valrange)
%
% Purpose: wrap vals around valrange
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/04/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

ind = find(vals>valrange(2));
vals(ind) = valrange(1) + vals(ind) - valrange(2);
ind = find(vals<valrange(1));
vals(ind) = valrange(2) + vals(ind) - valrange(1);

