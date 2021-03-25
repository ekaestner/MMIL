function vol = mmil_norm_vol(vol,normval,thresh)
%function vol = mmil_norm_vol(vol,[normval],[thresh])
%
% Purpose: normalize each voxel of input time-series volume to mean
%
% Required Input:
%   vol: timeseries (4D) volume
%     size should be [nx,ny,nz,nt]
%
% Optional Parameters:
%   normval: value to be assigned to mean
%     {defalut = 100}
%   thresh: when normalize to mean, set voxels with original values
%    less than this to 0
%     {default = 10}
%
% Created:  07/13/12 Don Hagler
% Last Mod: 11/01/12 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('normval','var') || isempty(normval), normval = 100; end;
if ~exist('thresh','var') || isempty(thresh), thresh = 10; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx,ny,nz,nt] = size(vol);
vec = reshape(vol,[nx*ny*nz,nt]);
vec_mean = mean(vec,2);
vec = normval*vec./(vec_mean*ones(1,nt)+eps);
vec(vec_mean<thresh,:) = 0;
vol = reshape(vec,size(vol));

