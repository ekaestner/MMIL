function volV = dti_rotate_vectors(volV,R)
%function volV_rot = dti_rotate_vectors(volV,R)
%
% Required Input:
%   volV: volume containing vector at each voxel
%     size = [nx,ny,nz,3]
%   R: 3x3 rotation matrix
%     (if 4x4 matrix is supplied, only the scale invariant
%      3x3 rotation component will be used)
%
% Output:
%   volV_rot: volume with tensors rotated by R
%     size = [nx,ny,nz,3]
%
% Created:  05/02/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

% check dimensions
nx = size(volV,1);
ny = size(volV,2);
nz = size(volV,3);
if length(size(volV))~=4 | any(size(volV)~=[nx,ny,nz,3])
  fprintf('%s: ERROR: volV must have size = [%d,%d,%d,3]\n',...
    mfilename,nx,ny,nz);
  return;
end;
nvox = nx*ny*nz;
volV = reshape(volV,nvox,3);

% make sure rotation matrix is 3x3 scale invariant
R = dti_M_to_R(R);
% account for RAS to LPH flips
M = dti_M_to_R(M_RAS_TO_LPH);
R = inv(M)*R*M;

% apply rotation and reshape to volume
volV = (R*volV')';
volV = reshape(volV,nx,ny,nz,3);

return;

