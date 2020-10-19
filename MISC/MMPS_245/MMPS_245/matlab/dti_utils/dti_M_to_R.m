function R = dti_M_to_R(M)
%function R = dti_M_to_R(M)
%
% Required Input:
%   M: 3x3 or 4x4 matrix
%     (if 4x4 matrix is supplied, only the scale invariant
%      3x3 rotation component will be used)
%
% Output:
%   R: 3x3 scale invariant rotation matrix
%
% Created:  11/05/10 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
R = [];

% check matrix
if all(size(M)==[3 3])
  R = M;
elseif all(size(M)==[4 4])
  R = M(1:3,1:3);
else
  error('M must be either 3x3 or 4x4');
end;

% scale rotation part of transformation matrix by voxel sizes
voxres = sqrt(sum(R.^2,1));
R = R ./ (ones(3,1)*voxres);


