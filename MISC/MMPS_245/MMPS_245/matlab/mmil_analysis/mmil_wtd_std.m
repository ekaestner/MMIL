function wtd_std = mmil_wtd_std(x,w,dim)
%function wtd_std = mmil_wtd_std(x,w,dim)
%
% Purpose: calculate weighted standard deviation
%
% Required Input:
%   x: vector of values
%   w: vector of weights
%
% Optional Input:
%   dim: dimension along which to calculate standard deviation
%     if empty, will use first non-singleton dimension
%     {default = []}
%
% Created:  10/22/12 by Don Hagler
% Last Mod: 10/23/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wtd_std = [];
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
  wtd_std = nan(x_sz);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wtd_mean = mmil_wtd_mean(x,w,dim);

wtd_std = sqrt((sum(w,dim)./(sum(w,dim).^2 - sum(w.^2,dim))) .* ...
          sum(w.*(bsxfun(@minus,x,wtd_mean).^2),dim));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% from: http://en.wikipedia.org/wiki/Weighted_mean
% and:  http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html

% sumw = sum(w)
% sumw2 = (sum(w)).^2
% sumsqw = sum(w.^2)
%
% wtd_std = sqrt((sumw/(sumw2 - sumsqw))*(sum(w*((x-wtd_mean(x,w)).^2)))

