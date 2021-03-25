function fs_write_label(vertices,fname,subj)
%function fs_write_label(vertices,fname,[subj])
%
% Required Input:
%   vertices: vector 1-based vertex numbers
%   fname: output file name
%
% Optional Input:
%   subj: FreeSurfer subject name
%     {default = 'fsaverage'}
%
% Created:   10/20/10 by Don Hagler
% Last Mod:  03/16/11 by Don Hagler
%

% based on write_label distributed with freesurfer

if ~mmil_check_nargs(nargin, 2), return; end;
if ~exist('subj','var') || isempty(subj), subj = 'fsaverage'; end;

npoints = length(vertices); 
lxyz = zeros(npoints,3); 
lvals  = zeros(npoints,1); 


% open as an ascii file
fid = fopen(fname, 'w') ;
if(fid == -1)
  error('could not open %s',fname);
end

fprintf(fid,'#!ascii label, from subject %s \n',subj);
fprintf(fid,'%d\n',npoints);

% Make sure they are npoints by 1 %
vertices = reshape(vertices,[npoints 1]) - 1; % change to zero base
lxyz   = reshape(lxyz,[npoints 3]);
lvals  = reshape(lvals,[npoints 1]);

l = [vertices lxyz lvals];
fprintf(fid,'%d %f %f %f %f\n',l') ;

fclose(fid) ;

