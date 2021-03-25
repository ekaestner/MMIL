function [rho,pval] = mmil_corr_mat(D1,D2)
%function [rho,pval] = mmil_corr_mat(D1,D2)
%
% Purpose: calculate correlation between two 2D matrices
%
% Required Input:
%   D1: matrix of size = [n,p] or [n,1]
%   D2: matrix of size = [n,p]
%
% Output:
%   rho: vector of correlation coefficients with size = [1,p]
%   pval: vector of p values with size = [1,p]
%
% Created:  02/24/14 by Don Hagler
% Last Mod: 02/24/14 by Don Hagler
%
% Note: based on code by Anders Dale
%

if ~mmil_check_nargs(nargin,2), return; end;
rho = [];
pval = [];

% eliminate rows with missing or invalid values
ivec = find(isfinite(sum(D1,2)+sum(D2,2)));
D1 = D1(ivec,:); D2 = D2(ivec,:);

if size(D1,2) == 1
  D1 = repmat(D1,1,size(D2,2));
end
mu1 = mean(D1,1);
mu2 = mean(D2,1);
D1 = D1 - ones(size(D1,1),1)*mu1;
D2 = D2 - ones(size(D2,1),1)*mu2;
norm1 = sqrt(sum(D1.^2,1));
norm2 = sqrt(sum(D2.^2,1));
rho = sum(D1.*D2,1)./norm1./norm2;

if nargout>1
  n = size(D1,1);
  t = rho.*sqrt((n-2)./(1-rho.^2));
  pval = 2*tcdf(-abs(t),n-2);
end;

