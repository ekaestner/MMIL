function [vol_res,M_res] = mmil_resample_vol(vol,M,varargin)
%function [vol_res,M_res] = mmil_resample_vol(vol,M,[options])
%
% Purpose: resamples a single-frame volume with specified
%   registration matrix and output dimensions and vox2ras matrix
%
% Required Input:
%   vol: volume matrix
%   M: 4x4 vox2ras matrix
%
% Optional Parameters:
%   'M_ref': 4x4 vox2ras matrix for reference volume
%     if not specified, use M
%     {default = []}
%   'nvox_ref': number of voxels [nx,ny,nz] for reference volume
%     if not specified, use size(vol)
%     {default = []}
%   'M_reg': registration matrix in LPH (scanner) coordinates applied to volume
%     {default = eye(4)}
%   'interpm': [0|1|2|3|4|5] interpolation method
%      0:nearest  1:linear  2:cubic  3:Key's spline  4:cubic spline  5: hamming sinc
%     {default = 2}
%   'bclamp' : [0|1] set negative values to zero
%     {default = 1}
%
% Created:  10/10/12 by Don Hagler
% Last Mod: 10/10/12 by Don Hagler
%

% based on vol_resample_pad, created 05/09/07 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
vol_res = []; M_res = [];
parms = mmil_args2parms(varargin, { ...
  'M',M,[],...
  'M_ref',[],[],...
  'nvox_ref',[],[],...
  'M_reg',eye(4),[],...
  'interpm',2,[0:5],...
  'bclamp',true,[false,true],...
...
  'npad',5,[],...
});

if isempty(parms.M_ref), parms.M_ref = parms.M; end;
if isempty(parms.nvox_ref), parms.nvox_ref = size(vol); end;

if length(size(vol))~=3 || length(parms.nvox_ref)<3
  error('input volumes must be 3D');
end;

if length(parms.nvox_ref)>3
  parms.nvox_ref = parms.nvox_ref(1:3);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol1 = ctx_mgh2ctx(zeros(parms.nvox_ref),parms.M_ref);
vol2 = ctx_mgh2ctx(vol,parms.M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol2r = [];

% pad vol2 with extra slices at top and bottom
tmp_vol2 = vol2;
tmp_vol2_mask = vol2;
nx2 = size(vol2.imgs,1);
ny2 = size(vol2.imgs,2);
nz2 = size(vol2.imgs,3);
tmp_vol2.imgs = zeros(nx2,ny2,nz2+2*parms.npad);
tmp_vol2_mask.imgs = zeros(nx2,ny2,nz2+2*parms.npad);

% pad vol1 with extra slices at top and bottom
nx1 = size(vol1.imgs,1);
ny1 = size(vol1.imgs,2);
nz1 = size(vol1.imgs,3);
vol1.imgs = zeros(nx1,ny1,nz1+2*parms.npad);

% duplicate edge slices
tmp_vol2.imgs(:,:,parms.npad+1:nz2+parms.npad) = vol2.imgs;
tmp_vol2_mask.imgs(:,:,parms.npad+1:nz2+parms.npad) = ones(size(vol2.imgs));
for i=1:parms.npad
  tmp_vol2.imgs(:,:,i) = vol2.imgs(:,:,1);
end;
for i=nz2+parms.npad+1:nz2+2*parms.npad
  tmp_vol2.imgs(:,:,i) = vol2.imgs(:,:,nz2);
end;

% account for difference in center of slices
vol1.Mvxl2lph_orig = vol1.Mvxl2lph;
vol1.Mvxl2lph = shift_M(vol1.Mvxl2lph,[nx1,ny1,nz1],[nx1,ny1,nz1+2*parms.npad]);
tmp_vol2.Mvxl2lph = shift_M(tmp_vol2.Mvxl2lph,...
                            [nx2,ny2,nz2],[nx2,ny2,nz2+2*parms.npad]);
tmp_vol2_mask.Mvxl2lph = tmp_vol2.Mvxl2lph;

% resample volumes
tmp_vol2r = vol_resample(tmp_vol2,vol1,parms.M_reg,parms.interpm,3,parms.bclamp);
tmp_vol2r_mask = vol_resample(tmp_vol2_mask,vol1,parms.M_reg,0);
vol2r = tmp_vol2r;

% remove padding slices
vol2r.imgs(tmp_vol2r_mask.imgs==0) = 0;
vol2r.imgs = squeeze(vol2r.imgs(:,:,parms.npad+1:nz1+parms.npad));
vol2r.Mvxl2lph = vol1.Mvxl2lph_orig;

% set NaNs to 0
vol2r.imgs(isnan(vol2r.imgs)) = 0;

[vol_res,M_res] = ctx_ctx2mgh(vol2r);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_shift = shift_M(M,nvox,nvox_shift);
  M_shift = M;
  M_shift(1:3,4) = M_shift(1:3,4) + M(1:3,:)*[nvox/2+1 1]' -...
                                    M_shift(1:3,:)*[nvox_shift/2+1 1]';
return;


