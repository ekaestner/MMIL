function MNIv = fs_ras2mni(RASv,T)
%function MNIv = fs_ras2mni(RASv,T)
%
% Purpose: convert FreeSurfer RAS to MNI "Talairach" coordinates
%
% Required Input:
%   RASv: Nx3 matrix of vertex coordinates in RAS space
%   T: 3x4 MNI Talairach transformation matrix (see fs_read_talxfm)
%
% Output:
%   MNIv: MNI "Talairach" coordinates
%
% Note: to adjust MNI Talairach coordinates to true Talairach with
%   Matthew Brett's "best guess" approximation, use fs_mni2tal
%   This function was based on Darren Webber's freesurfer_surf2tal
%
% Created:   11/03/11 by Don Hagler
% Last Mod:  11/03/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

% check size of RASv
if size(RASv,2)~=3, error('size of RASv must be Nx3'); end;
N = size(RASv,1);

% transpose RASv and add ones to the matrix
RASv = RASv';
RASv(4,:) = ones(1,N);

% Convert FreeSurfer RAS to Talairach coordinates
MNIv = T * RASv;

% return Nx3 matrix
MNIv = MNIv';

return

