function wtd_med = mmil_wtd_median(x,w,dim)
%function wtd_med = mmil_wtd_median(x,w,dim)
%
% Purpose: calculate weighted median
%
% Required Input:
%   x: vector of values
%   w: vector of weights
%
% Optional Input:
%   dim: dimension along which to calculate median
%     if empty, will use first non-singleton dimension
%     {default = []}
%
% Created:  10/22/12 by Don Hagler
% Last Mod: 10/23/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wtd_med = [];
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
  wtd_med = nan(x_sz);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put dimension of interest first
permvec = [1:length(size(x))];
permvec(1) = dim;
permvec(dim) = 1;

% reshape matrices to be 2-dimensional
tmp_x = permute(x,permvec);
tmp_w = permute(w,permvec);
tmp_sz = size(tmp_x);
nvals = tmp_sz(1);
nextra = prod(tmp_sz(2:end));
tmp_x = reshape(tmp_x,[nvals,nextra]);
tmp_w = reshape(tmp_w,[nvals,nextra]);

% sort values of x
[sortx,order] = sort(tmp_x,1);
ind_sort = order(:,1);
sortw = tmp_w(ind_sort,:);

midpoint = sum(sortw,1)/2;
csumw = cumsum(sortw,1);

% find median
wtd_med = zeros(1,nextra);
for i=1:nextra
  j = find(csumw(:,i)<=midpoint(i),1,'last');
  if isempty(j)
    wtd_med(i) = sortx(1,i);
  elseif j==nvals
    wtd_med(i) = sortx(j,i);
  elseif csumw(j,i) == midpoint(i)
    wtd_med(i) = mean(sortx([j j+1],i));
  else
    wtd_med(i) = sortx(j+1,i);
  end
end;

% reshape back, ipermute
wtd_med = reshape(wtd_med,[1,tmp_sz(2:end)]);
wtd_med = ipermute(wtd_med,permvec);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% created based on a post to matlab central
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/97571

