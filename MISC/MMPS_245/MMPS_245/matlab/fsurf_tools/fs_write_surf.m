function fs_write_surf(surf,fname)
% fs_write_surf - write a freesurfer surface file
% 
% fs_write_surf(surf,fname)
% 
% Writes the vertex coordinates (mm) and face lists to a surface file
%
% surf is a structure containing:
%   nverts: number of vertices
%   nfaces: number of faces (triangles)
%   faces:  vertex numbers for each face (3 corners)
%   vertices: x,y,z coordinates for each vertex
%
% created:        06/20/07 Don Hagler
% last modified:  01/08/10 Don Hagler
%
% see also: fs_read_surf, fs_write_trisurf
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRIANGLE_FILE_MAGIC_NUMBER  =  16777214;

coords = reshape(surf.vertices',3*surf.nverts,1);
faces = reshape(surf.faces',3*surf.nfaces,1) - 1; % revert to 0-based indices

user = getenv('USER');
if isempty(user), user = getenv('LOGNAME'); end;
if isempty(user), user = 'UNKNOWN'; end;

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;
if (fid < 0),
    str = sprintf('%s: could not open surface file %s.',mfilename,fname) ;
    error(str) ;
end

fprintf('%s: writing triangle file...',mfilename); tic;
fs_fwrite3(fid,TRIANGLE_FILE_MAGIC_NUMBER) ;
fprintf(fid,'created by %s on %s\n\n',user,datestr(now)); % creation date text line

fwrite(fid,surf.nverts,'int32'); % number of vertices
fwrite(fid,surf.nfaces,'int32'); % number of faces
fwrite(fid,coords,'float32');
fwrite(fid,faces,'int32');
t=toc; fprintf('done (%0.2f sec)\n',t);

fclose(fid) ;

return
