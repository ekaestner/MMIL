function vol2r = vol_resample_pad(vol2,vol1,M_vol1_to_vol2,interpm,bclamp)
%function vol2r = vol_resample_pad(vol2,vol1,M_vol1_to_vol2,[interpm],[bclamp])
%
% purpose: calls vol_resample with slice padding on ends of images
%
%   interpm : 0: Nearest 1:Linear 2:cubic (default)
%             3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%   bclamp : set negative value to zero default = true
%
% Created:  05/09/07 by Don Hagler
% Last Mod: 07/27/10 by Don Hagler
%

vol2r = [];

if nargin<3
  help(mfilename);
  return;
end;
if ~exist('interpm','var') | isempty(interpm), interpm = 2; end;
if ~exist('bclamp','var') | isempty(bclamp), bclamp = 1; end;

if length(size(vol2.imgs))~=3 || length(size(vol1.imgs))~=3
  error('input volumes must be 3D');
end;

npad = 5;

vol2r = [];

tmp_vol2 = vol2;
tmp_vol2_mask = vol2;
nx2 = size(vol2.imgs,1);
ny2 = size(vol2.imgs,2);
nz2 = size(vol2.imgs,3);
tmp_vol2.imgs = zeros(nx2,ny2,nz2+2*npad);
tmp_vol2_mask.imgs = zeros(nx2,ny2,nz2+2*npad);

nx1 = size(vol1.imgs,1);
ny1 = size(vol1.imgs,2);
nz1 = size(vol1.imgs,3);
vol1.imgs = zeros(nx1,ny1,nz1+2*npad);

tmp_vol2.imgs(:,:,npad+1:nz2+npad) = vol2.imgs;
tmp_vol2_mask.imgs(:,:,npad+1:nz2+npad) = ones(size(vol2.imgs));
for i=1:npad
  tmp_vol2.imgs(:,:,i) = vol2.imgs(:,:,1);
end;
for i=nz2+npad+1:nz2+2*npad
  tmp_vol2.imgs(:,:,i) = vol2.imgs(:,:,nz2);
end;

vol1.Mvxl2lph_orig = vol1.Mvxl2lph;
vol1.Mvxl2lph = shift_M(vol1.Mvxl2lph,[nx1,ny1,nz1],[nx1,ny1,nz1+2*npad]);
tmp_vol2.Mvxl2lph = shift_M(tmp_vol2.Mvxl2lph,[nx2,ny2,nz2],[nx2,ny2,nz2+2*npad]);
tmp_vol2_mask.Mvxl2lph = tmp_vol2.Mvxl2lph;

tmp_vol2r = vol_resample(tmp_vol2,vol1,M_vol1_to_vol2,interpm,3,bclamp);
tmp_vol2r_mask = vol_resample(tmp_vol2_mask,vol1,M_vol1_to_vol2,0);

vol2r = tmp_vol2r;
vol2r.Mvxl2lph = vol1.Mvxl2lph_orig;

vol2r.imgs(tmp_vol2r_mask.imgs==0) = 0;
vol2r.imgs = squeeze(vol2r.imgs(:,:,npad+1:nz1+npad));

vol2r.imgs(isnan(vol2r.imgs)) = 0;

return;


function M_shift = shift_M(M,nvox,nvox_shift);
  M_shift = M;
  M_shift(1:3,4) = M_shift(1:3,4) + M(1:3,:)*[nvox/2+1 1]' - M_shift(1:3,:)*[nvox_shift/2+1 1]';
return;


