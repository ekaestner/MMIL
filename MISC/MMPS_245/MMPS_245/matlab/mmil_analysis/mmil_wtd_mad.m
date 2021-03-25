function [wtd_mad,wtd_med] = mmil_wtd_mad(x,w,dim)
%function [wtd_mad,wtd_med] = mmil_wtd_mad(x,w,dim)
%
% Purpose: calculate weighted median absolute deviation (mad)
%
% Required Input:
%   x: vector of values
%
% Optional Input:
%   w: vector of weights
%     if empty, will weight all values equally
%     {default = []}
%   dim: dimension along which to calculate mad
%     if empty, will use first non-singleton dimension
%     {default = []}
%
% Output:
%   wtd_mad: weighted median absolute deviation
%   wtd_med: weighted median
%
% Created:  10/23/12 by Don Hagler
% Last Mod: 06/19/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wtd_mad = []; wtd_med = [];
if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('w','var'), w = []; end;
if ~exist('dim','var') || isempty(dim)
  dim = mmil_nonsingleton_dim(x);
end;

if isempty(w)
  w = ones(size(x));
end;

x_sz = size(x);  
w_sz = size(w);

if length(x_sz)~=length(w_sz) || any(x_sz~=w_sz)
  error('size of x does not match w');
end;

if all(w==0)
  x_sz(dim) = 1;
  wtd_mad = nan(x_sz);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wtd_med = mmil_wtd_median(x,w,dim);

ad = abs(bsxfun(@minus,x,wtd_med));

wtd_mad = mmil_wtd_median(ad,w,dim)/0.6745;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

