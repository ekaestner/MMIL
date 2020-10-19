function y = mmil_logit(x,invflag)
%function y = mmil_logit(x,invflag)
%
% Created:  01/08/15 by Don Hagler
% Last Mod: 01/08/15 by Don Hagler
%

% based on code by Anders Dale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = [];
if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('invflag','var') || isempty(invflag), invflag = 1; end;
minval = eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if invflag
  y = logit(x,invflag);
  y = (y-logit(logit(minval,0),1))/(1-2*minval);
else
  x = minval+x*(1-2*minval);
  y = logit(x,invflag);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = logit(x,invflag)
  if invflag
    y = exp(x)./(1+exp(x));
  else
    y = log(x./(1-x));
  end
return;

