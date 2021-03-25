function resp=rc_respfunc(amp,delay,rise,flat,fall,min_total)
%function resp=rc_respfunc(amp,delay,rise,flat,fall,min_total)
%
% Early Mod: 12/22/06 by Don Hagler
% Last Mod:  02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,6), return; end;

total_dur = delay+rise+flat+fall;
if (min_total>total_dur)
  total_dur = min_total;
end

resp = zeros(total_dur,1);
resp(delay:delay+rise) = linear(delay,0,delay+rise,amp);
resp(delay+rise+1:delay+rise+flat) = amp;
resp(delay+rise+flat+1:delay+rise+flat+fall) = ...
  linear(delay+rise+flat+1,amp,delay+rise+flat+fall,0);

resp = smooth(resp);

return;

function y=linear(x1,y1,x2,y2)
  if(x1>x2)
    y=[];
    return;
  end
  y=zeros(1+x2-x1,1);
  for x=1:1+x2-x1
    y(x)=x*(y2-y1)/(1+x2-x1) + y1;
  end
return;

