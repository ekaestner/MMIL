function qmat_rot = dti_rotate_qmat(qmat,R)
%function qmat_rot = dti_rotate_qmat(qmat,R)
% Required Input:
%  qmat: nx3 matrix of diffusion vectors
%
%   R: 3x3 rotation matrix to be applied to qmat
%     (if 4x4 matrix is supplied, only the scale invariant
%      3x3 rotation component will be used)
%
% Output:
%   qmat_rot: qmat rotated by R
%
% Notes:
%   Works for GE Scanner data acquired in RAS format
%   for other data format it needs testing.
%
% Created:  11/05/10 by Vijay Venkatraman
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

% make rotation matrix 3x3 if 4x4 was supplied
R = dti_M_to_R(R);

% adjust qmat for rotation component of registration
R = dti_M_to_R(M_LPH_TO_RAS)*R*dti_M_to_R(M_RAS_TO_LPH);

qmat_rot = (R*qmat')';

return;
