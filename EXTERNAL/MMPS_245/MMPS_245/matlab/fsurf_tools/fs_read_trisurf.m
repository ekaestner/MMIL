function [surf] = fs_read_trisurf(fname)
% fs_read_trisurf - read a freesurfer surface tri file
% 
% [surf] = fs_read_trisurf(fname)
% 
% Reads the vertex coordinates (mm) and face lists from a surface file
%
% surf is a structure containing:
%   nverts: number of vertices
%   nfaces: number of faces (triangles)
%   faces:  vertex numbers for each face (3 corners)
%   vertices: x,y,z coordinates for each vertex
%
% code for reading tri file from Nitin Bangera's tri2mat
%
% see also: fs_read_surf, fs_find_neighbors, fs_calc_triarea
%
% created:        05/09/06 Don Hagler
% last modified:  08/17/14 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


surf = [];

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

fid=fopen(fname,'rt');
surf.nverts=str2num(fgetl(fid));
surf.vertices=zeros(surf.nverts,3);
while (~feof(fid))
  for i=1:surf.nverts,  
    temp=str2num(fgetl(fid));
    surf.vertices(i,:)=temp(1,2:end);
  end
  surf.nfaces=str2num(fgetl(fid));
  surf.faces=zeros(surf.nfaces,3);
  for j=1:surf.nfaces,  
    temp1=str2num(fgetl(fid));
    surf.faces(j,:)=temp1(1,2:end);
  end
end

fclose(fid);

return
