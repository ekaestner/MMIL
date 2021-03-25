function [data_unwarp,time_unwarp] = rc_reverse_time_warp(data,regstruct,varargin)
%function [data_unwarp,time_unwarp] = rc_reverse_time_warp(data,regstruct)
%
% Required Input:
%   data: time series data; with size [ntpoints,nconds]
%   regstruct: registration struct, output of rc_landmark_reg
%
% Optional Parameters:
%    none
%
% Output:
%   data_unwarp: matrix the size of data, after applying reverse time warp
%   time_unwarp: matrix the size of data, with relative time shift
%
% Created:  03/30/16 by Don Hagler
% Last mod: 05/26/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_unwarp = [];
time_unwarp = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'none',[],[],...
  'smf',1e-5,[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if regstruct.bspline_padflag
  ntpoints = regstruct.ntpoints + regstruct.bspline_padflag*regstruct.bspline_npad;
  data_pad = zeros(ntpoints,regstruct.nconds);
  data_pad(regstruct.bspline_npad+1:regstruct.bspline_npad+regstruct.ntpoints,:) = data;
else
  ntpoints = regstruct.ntpoints;
  data_pad = data;
end;
data_unwarp = zeros(ntpoints,regstruct.nconds);
warpfun = eval_fd(regstruct.bspline_time,regstruct.warpfun_fd); % warping function
for c=1:regstruct.nconds
  wt = warpfun(:,c);
  v = data_pad(:,c);
  % make sure that wt is monotonic
  wt0 = wt;
  if ~isempty(find(diff(wt)<=0))
    dwt = diff(wt);
    max_time = max(cumsum(dwt));
    dwt(dwt<=0) = parms.smf;
    cs_dwt = cumsum(dwt);
    cs_dwt = max_time * cs_dwt / max(cs_dwt);
    wt = [wt(1);wt(1)+cs_dwt];
  end;
  if wt(end)>regstruct.bspline_time(end)
    wt(end) = regstruct.bspline_time(end);
  end;
  if wt(1)<regstruct.bspline_time(1)
    wt(1) = regstruct.bspline_time(1);
  end;
  try
    v_fd = data2fd(wt,v,regstruct.bspline_obj);
    uv = eval_fd(regstruct.bspline_time,v_fd);
  catch me
    fprintf('%s: WARNING: error with data2fd: %s\n',mfilename,me.message);
    keyboard
    v_fd = data2fd(wt,repmat(v,[1,2]),regstruct.bspline_obj);
    uv = eval_fd(regstruct.bspline_time,v_fd);
    uv = uv(:,1);
  end;
  data_unwarp(:,c) = uv;
end;
if regstruct.bspline_padflag
  data_unwarp = data_unwarp(regstruct.bspline_npad+1:regstruct.bspline_npad+regstruct.ntpoints,:);
end;

time_unwarp = repmat(regstruct.time,[1,regstruct.nconds])-regstruct.warpfun;

return;

