function volT_rot = dti_rotate_tensors(volT,R,volmask)
%function volT_rot = dti_rotate_tensors(volT,R,[volmask])
%
% Required Input:
%   volT: volume containing tensor at each voxel
%     size = [nx,ny,nz,3,3]
%   R: 3x3 rotation matrix
%     (if 4x4 matrix is supplied, only the scale invariant
%      3x3 rotation component will be used)
%
% Output:
%   volT_rot: volume with tensors rotated by R
%     size = [nx,ny,nz,3,3]
%
% Created:  04/11/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

volT_rot = [];

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('volmask','var'), volmask = []; end;
smf = 10^-10;

% check dimensions
nx = size(volT,1);
ny = size(volT,2);
nz = size(volT,3);
nvox = nx*ny*nz;
ndims = length(size(volT));

switch ndims
  case 5
    if any(size(volT)~=[nx,ny,nz,3,3])
      error('volT must have size = [%d,%d,%d,3,3]',nx,ny,nz);
    end;
    num_tensor_components = 9;
  case 4
    num_tensor_components = size(volT,4);
    switch num_tensor_components
      case {6,9}
      otherwise
        error('number of tensor components for volT must be 6 or 9 (is %d)',...
          num_tensor_components);
    end;
  otherwise
    error('wrong number of dimensions (%d) for volT',ndims);
end;

% make sure rotation matrix is 3x3 scale invariant
R = dti_M_to_R(R);
% account for RAS to LPH flips
M = dti_M_to_R(M_RAS_TO_LPH);
R = inv(M)*R*M;

volT = reshape(volT,[nvox,num_tensor_components]);

% find masked, non-zero voxels
tmp = sum(abs(volT),2);
if ~isempty(volmask)
  voxels = find(tmp>smf & volmask);
else
  voxels = find(tmp>smf);
end;
zvox = setdiff([1:nvox],voxels);
nvox_ROI = length(voxels);
fprintf('%s: %d selected voxels\n',mfilename,nvox_ROI);

volT(zvox,:,:) = 0;

% for each voxel, rotate tensor with R
for m=1:nvox_ROI
  k=voxels(m);
  T = squeeze(dti_components2tensor(volT(k,:)));
  T = R*T*R';
  volT(k,:) = dti_tensor2components(T,num_tensor_components);
end;

switch ndims
  case 5
    volT = reshape(volT,[nx,ny,nz,3,3]);
  case 4
    volT = reshape(volT,[nx,ny,nz,num_tensor_components]);
end;    
volT_rot = volT;

return;

