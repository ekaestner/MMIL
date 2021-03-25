function voxsize = fs_voxsize(M)
%function voxsize = fs_voxsize(M)
%
% Required Input:
%   M: 4x4 vox2ras matrix
%
% Created:  05/05/10 by Don Hagler
% Last Mod: 05/25/12 by Don Hagler
%

voxsize=[];

if ~mmil_check_nargs(nargin, 1), return; end;

if any(size(M)~=[4,4])
  error('M matrix must be 4x4');
end;

voxsize = sqrt(sum(M(1:3,1:3).^2,1));

