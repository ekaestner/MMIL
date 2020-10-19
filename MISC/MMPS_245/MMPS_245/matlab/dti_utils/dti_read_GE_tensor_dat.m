function qmat = dti_read_GE_tensor_dat(ndirs,nb0,fname)
%function qmat = dti_read_GE_tensor_dat(ndirs,[nb0],[fname])
%
% Purpose: read diffusion vectors from GE tensor.dat file
%
% Input:
%  ndirs: number of diffusion directions
%  nb0 (optional): number of b0 images
%   default: 1
%  fname (optional): file name containing diffusion vectors
%   default: 'GE_tensor_14x.dat'
%
% Output:
%  qmat: matrix of diffusion direction vectors
%   size will be ndirs+nb0 x 3
%   (0 0 0 will be prepended as first nb0 rows)
%
% Created:  08/04/07 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

qmat = [];

if ~mmil_check_nargs(nargin,1), return; end;

min_ndirs = 6;
max_ndirs = 129;

if ~exist('nb0','var') nb0=[]; end;
if isempty(nb0) nb0=1; end;
if ~exist('fname','var'), fname=[]; end;
if isempty(fname) fname='GE_tensor_14x.dat'; end;

if ndirs>max_ndirs
  error('ndirs must be less than %d',max_ndirs);
end;
if ndirs<min_ndirs
  error('ndirs must be greater than %d',min_ndirs);
end;

fid=fopen(fname,'rt');
if fid==-1
  error('failed to open %s for reading',fname);
end;

qmat = zeros(ndirs+nb0,3);
qmat(1:nb0,:) = zeros(nb0,3);
N=min_ndirs;
while (~feof(fid))
  temp=fgetl(fid);
  if strcmp(temp(1),'#'), continue; end;
  N = sscanf(temp,'%d');
  if length(N)==1 && N==ndirs, break; end;
end

for i=1:ndirs
  temp=fgetl(fid);
  A = sscanf(temp,'%f');
  if length(A)~=3
    error('failed to read direction vector from %s',fname);
  end;
  qmat(nb0+i,:)=A';
end;

fclose(fid);

