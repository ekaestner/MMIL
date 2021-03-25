function [MNIv,TALv] = freesurfer_surf2tal(SURFv,TalairachXFM)

% freesurfer_surf2tal - convert freesurfer RAS to Talairach coordinates
%
% [MNIv,TALv] = freesurfer_surf2tal(SURFv,TalairachXFM)
%
% The input 'SURFv' are from freesurfer_read_surf; they must be Nx3
% vertex coordinates in freesurfer RAS coordinates. The second
% input is the matrix from the talairach.xfm file (see
% freesurfer_read_talxfm).
%
% This function converts the RAS coordinates into the MNI Talairach
% coordinates (MNIv).  It can also return the adjusted Talairach
% coordinates (TALv), based on the notes by Matthew Brett. This function
% calls mni2tal to apply the conversion from MNI template space to 'true'
% Talairach space.  See Matthew Brett's online discussion of this topic,
% http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml
%
% The surface RAS coordinates are arranged into a column vector,
% SURFv = [R A S 1]', which is multiplied by the talairach.xfm
% matrix, TalairachXFM, to give the MNI Talairach coordinates:
%
% MNIv = TalairachXFM * SURFv;
% TALv = mni2tal(MNIv);
%
% The return matrices are Nx3
%

% $Revision: 1.3 $ $Date: 2005/07/05 23:40:14 $

% Licence:  GNU GPL, no express or implied warranties
% History:  11/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% transpose SURFv and add ones to the matrix

Nvertices = size(SURFv,1);

SURFv = SURFv';
SURFv(4,:) = ones(1, Nvertices);

%--------------------------------------------------------------------------
% Convert FreeSurfer RAS to Talairach coordinates

ver = '$Revision: 1.3 $';
fprintf('FREESURFER_SURF2TAL [v %s]\n',ver(11:15));

MNIv = TalairachXFM * SURFv;
MNIv = MNIv(1:3,:)'; % return Nx3 matrix

if nargout > 1,
    TALv = mni2tal(MNIv')';
end

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Email from
% Yasunari Tosa, Ph.D.               Email: tosa@nmr.mgh.harvard.edu
% NMR Ctr, Mass. General Hospital    TEL: 617-726-4050 Building 149, 13th Street
% Charlestown, MA 02129
%
%
% I think that the label file uses "surfaceRAS" values, not "RAS".
% Therefore you have to translate in two steps
%
%          surfaceRAS --------->  RAS --------------> talairachRAS
%                           translation             talairach.xfm
%
% where translation affine xfm is given by (I hope you can read
% it as 4x1 vector = (4 x 4 translation matrix) x (4 x 1 vector) )
%
%         RAS              translation     label_position
%           [ x_r ]     =   [ 1 0 0 c_r  ]  [label_r]
%           [ x_a ]          [ 0 1 0 c_a ]  [label_a]
%           [ x_s ]          [ 0 0 1 c_s  ]  [label_s]
%           [ 1    ]          [ 0 0 0   1    ]  [   1      ]
%
% i.e.  x_r = label_r  + c_r     (actual RAS position is shifted by c_(ras) value).
%       x_a = label_a + c_a
%       x_s = label_s + c_s
%
% Then multiply RAS by talairach.xfm (4 x 4) gives the talairachRAS position.