function MNIv = fs_tal2mni(TALv)
%function MNIv = fs_tal2mni(TALv)
%
% Purpose: transform Talairach coordinates to
%   MNI "Talairach" space (Matthew Brett's best guess)
%
% Input:
%   TALv: Nx3 matrix of Talairach coordinates
%
% Output:
%   MNIv: corresponding MNI coordinates
%
% Created:   02/02/01 by Matthew Brett
% Early Mod: 11/17/04 by Darren Webber
% Last Mod:  11/03/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Revision: 1.2 $ $Date: 2004/11/17 21:04:32 $
% Licence:  GNU GPL, no express or implied warranties
% Matthew Brett 2/2/01, matthew.brett@mrc-cbu.cam.ac.uk
% modified 02/2003, Darren.Weber_at_radiology.ucsf.edu
%                   - swapped inv() for slash equivalent
%                   - removed dependence on spm_matrix
% See also, MNI2TAL & the best guess discussion at
% http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% check size of TALv
if size(TALv,2)~=3, error('size of TALv must be Nx3'); end;
N = size(TALv,1);

% transpose TALv and add ones to the matrix
TALv = TALv';
TALv(4,:) = ones(1,N);

% Transformation matrices, different zooms above/below AC
M2T = fs_mni2tal_matrix;
tmp = TALv(3,:) < 0;  % 1 if below AC
TALv(:,  tmp) = (M2T.rotn * M2T.downZ) \ TALv(:,  tmp);
TALv(:, ~tmp) = (M2T.rotn * M2T.upZ  ) \ TALv(:, ~tmp);
MNIv = TALv(1:3, :);

% return Nx3 matrix
MNIv = MNIv';

return

