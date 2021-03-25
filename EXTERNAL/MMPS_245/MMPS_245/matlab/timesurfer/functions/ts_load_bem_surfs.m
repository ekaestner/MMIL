function [verts,faces]=ts_load_bem_surfs(bem_surf_files);
%function [verts,faces]=ts_load_bem_surfs(bem_surf_files);
%
% ts_load_bem_surfs: load FreeSurfer compatible tri files for BEM calculations
%
% Input:
%   bem_surf_files: cell array of 3 names of FreeSurfer compatible tri files
%     files should be in this order:
%       (1) inner skull
%       (2) outer skull
%       (3) outer scalp
%     optionally, can supply a single file name for the inner skull if
%       1-shell BEM is used
%
% Output:
%   verts: cell array of 3 matrices containing vertex coordinates
%     nverts x 3 (x,y,z coordinates)
%     If bem_surf_files is not a cell array, but instead just a file name
%       for the inner skull, verts will be just a matrix of coordinates
%       for the inner skull.
%   faces: cell array of 3 matrices containing triangle vertex numbers
%
% Created:  08/02/06 by Don Hagler
% Last Mod: 08/05/09 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;

verts = [];
faces = [];

if isempty(bem_surf_files)
  fprintf('%s: error: bem_surf_files is empty\n',mfilename);
  return;
end;
if iscell(bem_surf_files)
  num_bem_surfs=length(bem_surf_files);
  not_cell_flag=0;
else
  bem_surf_files={bem_surf_files};
  num_bem_surfs=1;
  not_cell_flag=1;
end;

for s=1:num_bem_surfs
  fname = bem_surf_files{s};
  if exist(fname,'file')
    [tmp_verts,tmp_faces]=ts_load_tri_file(fname);
    verts{s} = tmp_verts;
    faces{s} = tmp_faces;
  else
    error('file %s not found',fname);
  end;
end;

if not_cell_flag
  verts = verts{1};
  faces = faces{1};
end;

return;


