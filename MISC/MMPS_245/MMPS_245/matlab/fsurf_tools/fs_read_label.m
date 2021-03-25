function vertices = fs_read_label(fname)
%function vertices = fs_read_label(fname)
%
% purpose: reads from an ascii 'label' file
%
% Input:
%   fname: file name of label file
%
% Output:
%	  vertices: vector of vertex indices (1-based) for each vertex in label
%
% modified from freesurfer's read_label_old.m
%
% created:  09/07/06 Don Hagler
% last mod: 07/14/08 Don Hagler
%

vertices = [];
if (~mmil_check_nargs(nargin, 1)), return; end;

fid = fopen(fname, 'rt') ;
if (fid < 0)
  error('unable to open label file %s',fname);
end

s = fgetl(fid);
s = fgetl(fid);
[num_vert_label] = sscanf(s,'%d');
vertices = zeros(num_vert_label,1);
vert_data = zeros(5);
for vert = 1:1:num_vert_label,
  s = fgetl(fid);
  vert_data = sscanf(s,'%d %f %f %f %f');
  vertices(vert) = vert_data(1:1);
end;
fclose(fid);

vertices=vertices+1; % vertex numbers start at 0, but for matlab, index from 1

return;

