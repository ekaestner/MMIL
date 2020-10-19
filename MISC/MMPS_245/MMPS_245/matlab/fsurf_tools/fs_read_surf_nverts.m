function [nverts,nfaces] = fs_read_surf_nverts(fname)
%function [nverts,nfaces] = fs_read_surf_nverts(fname)
%
% Purpose: read number of vertices from a freesurfer surface file
% 
% Required Input:
%   fname: name of freesurfer binary surface file
%
% Output:
%   nverts: number of vertices
%   nfaces: number of faces (triangles)
%
% code for reading surfaces taken from Darren Weber's freesurfer_read_surf
%
% see also: fs_read_surf
%
% created:        06/11/06 Don Hagler
% last modified:  10/31/11 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRIANGLE_FILE_MAGIC_NUMBER  =  16777214 ;
QUAD_FILE_MAGIC_NUMBER      =  16777215 ;

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0),
  error(sprintf('could not open surface file %s.',fname));
end

magic = fs_fread3(fid) ;

if (magic == QUAD_FILE_MAGIC_NUMBER),
  nverts = fs_fread3(fid) ;
  nfaces = fs_fread3(fid) ;
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER),
  tline = fgets(fid); % read creation date text line
  tline = fgets(fid); % read info text line
  nverts = fread(fid, 1, 'int32') ; % number of vertices
  nfaces = fread(fid, 1, 'int32') ; % number of faces
else
  error(sprintf('unknown magic number in surface file %s.',fname));
end

return
