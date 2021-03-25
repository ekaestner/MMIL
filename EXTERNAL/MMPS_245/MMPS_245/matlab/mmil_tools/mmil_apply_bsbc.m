function vol = mmil_apply_bsbc(vol,xvals,yvals,varargin)
%function vol = mmil_apply_bsbc(vol,xvals,yvals,[options])
%
% Purpose: apply intensity transformation to match
%     images between scanners / scan parameters
%   using output from mmil_estimate_bsbc
%
% Required Input:
%   vol: input volume (3D matrix)
%   xvals: vector of original intensity values
%   yvals: vector of transformed intensity values
%
% Optional Parameters:
%   'interp_method': interpolation method (interp1)
%     e.g. 'nearest','linear','spline','cubic'
%     {default = 'linear'}
%   'zclamp_flag': [0|1] set negative values to zero
%     {default = 1}
%
% Output:
%   vol: output volume with transformed intensities
%
% Created:  01/08/15 by Don Hagler
% Last Mod: 01/08/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

parms = mmil_args2parms(varargin,{...
  'interp_method','linear',{'nearest','linear','spline','cubic'},...
  'zclamp_flag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol = interp1(xvals,yvals,vol,parms.interp_method,'extrap'); 

if parms.zclamp_flag
  vol(vol<0) = 0;
end;

