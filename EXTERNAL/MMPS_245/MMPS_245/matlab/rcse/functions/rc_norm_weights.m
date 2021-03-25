function nw = rc_norm_weights(w,maxflag)
%function nw = rc_norm_weights(w,maxflag)
%
% Purpose: normalize vector of weights
%
% Required Input:
%   w: vector of weights
%
% Optional Input:
%   maxflag: [0|1] whether to normalize by max value
%     otherwise normalize by sum
%     {default = 0}
%
% Early Mod: 12/22/06 by Don Hagler
% Last Mod:  02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('maxflag','var') || isempty(maxflag), maxflag = 0; end;

nw = w;

if isempty(w), return; end;

if maxflag
  f = max(w);
else
  f = sum(abs(w));
end;

if(f~=0)
  nw = w/f;
end

return;

