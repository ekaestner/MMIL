function TALv = fs_mni2tal(MNIv)
%function TALv = fs_mni2tal(MNIv)
%
% Purpose: transform MNI "Talairach" coordinates
%   to true Talairach space (Matthew Brett's best guess)
%
% Input:
%   MNIv: Nx3 matrix of  MNI Talairach coordinates
%
% Output:
%   TALv: approximate true Talairach coordinates
%
% Created:   10/08/99 by Matthew Brett
% Early Mod: 02/01/03 by Darren Webber
% Last Mod:  11/03/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml
% $Revision: 1.2 $ $Date: 2004/11/17 21:04:27 $
% Licence:  GNU GPL, no express or implied warranties
% Matthew Brett 10/8/99, matthew.brett@mrc-cbu.cam.ac.uk
% modified 02/2003, Darren.Weber_at_radiology.ucsf.edu
%                   - removed dependence on spm_matrix and
%                     abstracted the matrix in case it changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~mmil_check_nargs(nargin,1), return; end;

% check size of MNIv
if size(MNIv,2)~=3, error('size of MNIv must be Nx3'); end;
N = size(MNIv,1);

% transpose MNIv and add ones to the matrix
MNIv = MNIv';
MNIv(4,:) = ones(1,N);

% Transformation matrices, different zooms above/below AC
M2T = fs_mni2tal_matrix;
tmp = MNIv(3,:) < 0;  % 1 if below AC
MNIv(:,  tmp) = (M2T.rotn * M2T.downZ) * MNIv(:,  tmp);
MNIv(:, ~tmp) = (M2T.rotn * M2T.upZ  ) * MNIv(:, ~tmp);
TALv = MNIv(1:3, :);

% return Nx3 matrix
TALv = TALv';

return
