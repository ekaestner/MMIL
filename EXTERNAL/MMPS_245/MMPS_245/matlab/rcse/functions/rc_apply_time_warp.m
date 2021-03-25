function [data_warp,time_warp] = rc_apply_time_warp(data,regstruct,varargin)
%function [data_warp,time_warp] = rc_apply_time_warp(data,regstruct)
%
% Required Input:
%   data: time series data; with size [ntpoints,nconds]
%   regstruct: registration struct, output of rc_landmark_reg
%
% Optional Parameters:
%    none
%
% Output:
%   data_warp: matrix the size of data, after applying time warp
%   time_warp: matrix the size of data, with relative time shift
%
% Created:  03/30/16 by Don Hagler
% Last mod: 03/30/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_warp = [];
time_warp = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'none',[],[],...
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
data_warp = zeros(ntpoints,regstruct.nconds);
warpfun = eval_fd(regstruct.bspline_time,regstruct.warpfun_fd); % warping function
for c=1:regstruct.nconds
  wt = warpfun(:,c);
  wt = max(min(wt,max(regstruct.bspline_time)),min(regstruct.bspline_time));
  v = data_pad(:,c);
  v_fd = data2fd(regstruct.bspline_time,v,regstruct.bspline_obj);
  data_warp(:,c) = eval_fd(wt,v_fd);
end;
if regstruct.bspline_padflag
  data_warp = data_warp(regstruct.bspline_npad+1:regstruct.bspline_npad+regstruct.ntpoints,:);
end;

time_warp = regstruct.warpfun-repmat(regstruct.time,[1,regstruct.nconds]);


return;

