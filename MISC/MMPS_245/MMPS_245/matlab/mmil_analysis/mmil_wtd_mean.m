function wtd_mean = mmil_wtd_mean(x,w,dim)
%function wtd_mean = mmil_wtd_mean(x,w,dim)
%
% Purpose: calculate weighted mean
%
% Required Input:
%   x: vector of values
%   w: vector of weights
%
% Optional Input:
%   dim: dimension along which to calculate mean
%     if empty, will use first non-singleton dimension
%     {default = []}
%
% Created:  10/22/12 by Don Hagler
% Last Mod: 10/23/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wtd_mean = [];
if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('dim','var') || isempty(dim)
  dim = mmil_nonsingleton_dim(x);
end;

x_sz = size(x);
w_sz = size(w);

if length(x_sz)~=length(w_sz) || any(x_sz~=w_sz)
  error('size of x does not match w');
end;

if all(w==0)
  x_sz(dim) = 1;
  wtd_mean = nan(x_sz);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wtd_mean = sum(x.*w,dim) ./ sum(w,dim);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
