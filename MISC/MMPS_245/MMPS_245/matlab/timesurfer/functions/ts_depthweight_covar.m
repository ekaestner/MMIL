function R=ts_depthweight_covar(A,p,xyz_flag)
%function R=ts_depthweight_covar(A,p,[xyz_flag])
%
% create depth weighting source covariance matrix
%
% Required Parameters:
%   A: forward matrix (sensors by sources)
%
% Optional Paramters:
%   p: depth weighting parameter
%     {default: 0.5}
%   xyz_flag: [0|1] whether A contains x, y, and z components for each source
%     {default: 1}
%
% Output:
%   R: sparse source covariance matrix square in number of sources
%
% Created:  06/11/09 by Don Hagler
% Last Mod: 06/11/09 by Don Hagler
%

R = [];
if (~mmil_check_nargs(nargin,1)) return; end;

if ~exist('p','var') | isempty(p), p=0.5; end;
if ~exist('xyz_flag','var') | isempty(xyz_flag), xyz_flag=true; end;

[num_sensors,num_sources]=size(A);

R = speye(num_sources);

if xyz_flag
  num_dips = floor(num_sources/3);
  for k=1:num_dips
    a1 = A(:,3*k-2);
    a2 = A(:,3*k-1);
    a3 = A(:,3*k);
    f = (a1'*a1 + a2'*a2 + a3'*a3)^-p;
    R(3*k-2,3*k-2) = f;
    R(3*k-1,3*k-1) = f;
    R(3*k,3*k) = f;
  end;
else
  for k=1:num_sources
    a = A(:,k);
    f = (a'*a)^-p;
    R(k,k) = f;
  end;
end;


