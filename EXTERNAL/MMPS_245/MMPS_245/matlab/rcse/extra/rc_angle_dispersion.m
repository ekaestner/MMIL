function [amean,dmean,vmean] = rc_angle_dispersion(vmat,w)
%function [amean,dmean,vmean] = rc_angle_dispersion(vmat,w)
%
% Purpose: calculate dispersion of vectors
%
% Required Input:
%   vmat: matrix of 3D vectors with size = [3,n]
%
% Optional Input:
%   w: weights vector
%
% Output:
%   amean: mean angle between each vector and mean vector
%   dmean: mean dot product between vector and mean vector
%   vmean: mean vector
%
% Created:   01/15/13 by Don Hagler
% Last Mod:  01/16/13 by Don Hagler
%

amean = 0;
dmean = 1;
vmean = vmat;

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('w','var'), w = []; end;

if size(vmat,1) ~= 3
  error('vmat matrix must be 3 x n');
end;
n = size(vmat,2);
if isempty(w)
  w = ones(1,n);
else
  w = mmil_rowvec(w);
end;
if length(w) ~= n
  error('%s: length of w must match 2nd dim of vmat\n',mfilename);
end;

if n==1, return; end;

% normalize each vector to make sure they are all unit vectors
vn = sqrt(sum(vmat.^2,1));
if any(vn~=1)
  vmat = bsxfun(@rdivide,vmat,vn);
end;

% calculate mean vector
wmat = repmat(w,[3,1]);
vmean = mmil_wtd_mean(vmat,wmat,2);
vmean = vmean / sqrt(sum(vmean.^2));

dmean = 0;
for k=1:length(w)
  dmean = dmean + dot(vmean,vmat(:,k))*w(k);
end;
dmean = dmean/sum(w);

amean = acos(dmean)*180/pi;

return;

