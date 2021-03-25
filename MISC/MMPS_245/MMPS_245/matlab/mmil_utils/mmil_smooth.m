function v = mmil_smooth(v,sig)
%function v = mmil_smooth(v,sig)
%
% Purpose: applies Gaussian blur to 1D vector or 2D matrix
%
% Required Input:
%   v: 1D vector or 2D matrix
%   sig: Gaussian blurring kernel sigma
%
% Created:  07/22/11 by Don Hagler
% Last Mod: 07/22/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

if sig
  h = fspecial('gaussian',sig,sig);
  v = imfilter(v,h,'replicate');
end;
