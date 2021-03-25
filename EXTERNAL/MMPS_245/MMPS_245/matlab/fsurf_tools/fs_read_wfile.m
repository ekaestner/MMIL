function [weights,vertices] = fs_read_wfile(fname)
%
% [weights,vertices] = read_wfile(fname)
% reads from a binary 'w' file
%	fname - name of file to read from
%
%	weights - vector of values for each non-zero vertex
%	vertices - vector of vertex indices for each non-zero vertex
%
% modified from freesurfer's fast_read_wfile.m
%
% created:     ?      by Don Hagler
% last mod: 08/20/07  by Don Hagler
%

weights = [];
vertices = [];

if(nargin ~= 1)
  help(mfilename);
  return;
end

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
  error('failed to open w file %s',fname);
end

fread(fid, 1, 'int16') ;  % Skip ilat
vnum = fs_fread3(fid) ; % Number of non-zero values
vertices = zeros(vnum,1) ;
weights = zeros(vnum,1) ;
for i=1:vnum
  vertices(i) = fs_fread3(fid) ;
  weights(i) = fread(fid, 1, 'float') ;
end

fclose(fid) ;

vertices=vertices+1; % vertex numbers start at 0, but for matlab, index from 1

return;

