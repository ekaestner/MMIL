function RASv = fs_mni2ras(MNIv,T)
%function RASv = fs_mni2ras(MNIv,T)
%
% Purpose: convert MNI Talairach coordinates to FreeSurfer RAS
%
% Required Input:
%   MNIv: Nx3 matrix of MNI Talairach coordinates
%   T: 3x4 MNI Talairach transformation matrix (see fs_read_talxfm)
%
% Output:
%   RASv: RAS coordinates
%
% Note: to adjust true Talairach coordinates to MNI Talairach with
%   Matthew Brett's "best guess" approximation, use fs_tal2mni
%
%
% Created:   11/03/11 by Don Hagler
% Last Mod:  02/24/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

% check size of MNIv
if size(MNIv,2)~=3, error('size of MNIv must be Nx3'); end;
N = size(MNIv,1);

% transpose MNIv and add ones to the matrix
MNIv = MNIv';
MNIv(4,:) = ones(1,N);

% invert T
tmp = eye(4);
tmp(1:3,:) = T;
Tinv = inv(tmp);
Tinv = Tinv(1:3,:);

% Convert MNI Talairach to FreeSurfer RAS coordinates
RASv = Tinv * MNIv;

% return Nx3 matrix
RASv = RASv';

return


