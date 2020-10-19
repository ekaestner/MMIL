function fs_convert_fsurf2tri(fname_fsurf,fname_tri,nverts)
%function fs_convert_fsurf2tri(fname_fsurf,fname_tri,[nverts])
%
% Purpose: convert freesurfer binary format surface file
%   to ascii tri file format
%
% Required:
%   fname_fsurf: input freesurfer format surface file
%   fname_tri: output tri format surface file
%
% Optional:
%   nverts: number of vertices in output surface file
%
% created:        07/16/08 Don Hagler
% last modified:  02/24/10 Don Hagler
%
% see also: fs_read_surf, fs_write_trisurf
%

if (~mmil_check_nargs(nargin, 2)), return; end;
if ~exist('nverts','var'), nverts = []; end;

surf = fs_read_surf(fname_fsurf);

if ~isempty(nverts)
  % reduce number of vertices to nverts (approximately)
  R = nverts/surf.nverts;
  [nf,nv] = reducepatch(surf.faces,surf.vertices,R);
  surf.faces = nf;
  surf.vertices = nv;
  surf.nverts = size(nv,1);
  surf.nfaces = size(nf,1);
end;

fs_write_trisurf(surf,fname_tri);


