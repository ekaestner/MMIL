function [lh_dip_info,rh_dip_info]=ts_dip_info(subj,subjdir,surf,angle_flag)
%function [lh_dip_info,rh_dip_info]=ts_dip_info(subj,[subjdir],[surf],[angle_flag])
%
% Purpose: load FreeSurfer surface and return lh_dip_info and rh_dip_info
%   containing x,y,z location and nx,ny,nz normal vector for each vertex
%   of left and right cortical hemispheres
%
% Required Parameters:
%   subj: FreeSurfer subject name (name of recon directory)
%
% Optional Parameters:
%   subjdir: root directory containing FreeSurfer subjects
%     {default = getenv('SUBJECTS_DIR')}
%   surf: name of surface file (e.g. 'white', 'pial')
%     {default = 'white'}
%   angle_flag: [0|1] weight face normals by face angles
%     {default = 1}
%
% Created:    04/05/10 by Don Hagler
% Last Mod:   06/25/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lh_dip_info = [];
rh_dip_info = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('subjdir','var') | isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end;
end;
if ~exist('surf','var') | isempty(surf), surf = 'white'; end;
if ~exist('angle_flag','var') | isempty(angle_flag), angle_flag = 1; end;
hemilist = {'lh','rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for h=1:length(hemilist)
  hemi = hemilist{h};
  fname_surf = [subjdir '/' subj '/surf/' hemi '.white'];

  surf = fs_read_surf(fname_surf);
  normals = ts_calc_surf_normals(surf,angle_flag);

  dip_info = cat(2,surf.vertices,normals)';
  switch hemi
    case 'lh'
      lh_dip_info = dip_info;
    case 'rh'
      rh_dip_info = dip_info;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: calc_surf_normals_old requires "Projects" toolbox
%       no longer used because of a problem with ico-resampled surfaces

function normals = calc_surf_normals_old(surf)
  msurf = preprocessQ(surf);
  normals = computeNormals(msurf);
return;

