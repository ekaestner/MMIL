function fs_write_trisurf(surf,fname)
% fs_write_trisurf - write a freesurfer surface tri file
% 
% fs_write_trisurf(surf,fname)
% 
% Writes the vertex coordinates (mm) and face lists to a surface file
%
% surf is a structure containing:
%   nverts: number of vertices
%   nfaces: number of faces (triangles)
%   faces:  vertex numbers for each face (3 corners)
%   vertices: x,y,z coordinates for each vertex
%
% created:        07/13/06 Don Hagler
% last modified:  01/08/10 Don Hagler
%
% see also: fs_read_trisurf, fs_write_surf
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(fname,'wt');
if fid==-1
  fprintf('%s: ERROR: unable to create file %s\n',mfilename,fname);
  return;
end;

fprintf(fid,'%d\n',surf.nverts);
for i=1:surf.nverts
  fprintf(fid,'%d %f %f %f\n',...
    i,surf.vertices(i,1),surf.vertices(i,2),surf.vertices(i,3));
end;

fprintf(fid,'%d\n',surf.nfaces);
for i=1:surf.nfaces
  fprintf(fid,'%d %d %d %d\n',...
    i,surf.faces(i,1),surf.faces(i,2),surf.faces(i,3));
end;

fclose(fid);

return
