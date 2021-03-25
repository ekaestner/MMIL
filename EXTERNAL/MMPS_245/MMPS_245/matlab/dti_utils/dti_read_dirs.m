function qmat = dti_read_dirs(dirfile,nb0)
%function qmat = dti_read_dirs(dirfile,[nb0])
%
% Input:
%  dirfile: ascii file containing 3 columns
%   with a row for each diffusion direction
%  nb0 (optional): number of b0 images
%   default: 1
%
% qmat will be returned with 0 0 0 prepended as
%   first nb0 rows, so dirfile should not include
%   this
%
% Created:  10/21/06 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

maxndirs = 200;

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('nb0','var')
  nb0=1;
end;

qmat = [];

fid=fopen(dirfile,'rt');
qmat = zeros(maxndirs,3);
d=0;
while (~feof(fid))
  temp=fgetl(fid);
  A = sscanf(temp,'%f');
  if length(A)>3, break; end;
  d = d+1;
  qmat(d,:)=A';
end
fclose(fid);
ndirs = d;
qmat = qmat(1:ndirs,:);
qmat = cat(1,zeros(nb0,3),qmat);

