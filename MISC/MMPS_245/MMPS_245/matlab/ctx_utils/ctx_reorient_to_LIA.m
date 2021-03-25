function ctx_vol_out = ctx_reorient_to_LIA(ctx_vol,trans_vec)
%function ctx_vol_res = ctx_reorient_to_LIA(ctx_vol,[trans_vec])
%
% Purpose: reorient volume to LIA to match freesurfer's mri_convert
%   without chopping off the back of the brain
%
% Required:
%   ctx_vol: volume in ctx structure
%
% Optional:
%   trans_vec: vector of [x,y,z] translation in vox2ras matrix
%     if empty, calculate from input vox2ras
%     {default = []}
%
% Last Mod: 07/02/08 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;

if ~exist('trans_vec','var'), trans_vec=[]; end;

volsz = size(ctx_vol.imgs);

M = ctx_vol.Mvxl2lph(1:3,1:3);
%ctr_lph = ctx_vol.Mvxl2lph(1:3,:)*[(volsz+1)/2 1]';
%ctr_lph = ctx_vol.Mvxl2lph(1:3,:)*[(volsz-1)/2 1]'; % compensate for mri_convert wrongness
ctr_lph = ctx_vol.Mvxl2lph(1:3,:)*[(volsz)/2 1]';

Lcorrvec = [1 0 0]*M./sqrt(sum(abs(M).^2,1));
Icorrvec = [0 0 -1]*M./sqrt(sum(abs(M).^2,1));
Acorrvec = [0 -1 0]*M./sqrt(sum(abs(M).^2,1));

[Lcorr,Ldim] = max(abs(Lcorrvec));
[Icorr,Idim] = max(abs(Icorrvec));
[Acorr,Adim] = max(abs(Acorrvec));

ctx_vol_out.imgs = ctx_vol.imgs;
ctx_vol_out.Mvxl2lph = ctx_vol.Mvxl2lph;
ctx_vol_out.imgs = permute(ctx_vol.imgs,[Ldim,Idim,Adim]);
ctx_vol_out.Mvxl2lph = ctx_vol.Mvxl2lph(:,[Ldim,Idim,Adim,4]);

if Lcorrvec(Ldim)<0
  ctx_vol_out.imgs = flipdim(ctx_vol_out.imgs,1);
  ctx_vol_out.Mvxl2lph(:,1) = -1*ctx_vol_out.Mvxl2lph(:,1);
end
if Icorrvec(Idim)<0
  ctx_vol_out.imgs = flipdim(ctx_vol_out.imgs,2);
  ctx_vol_out.Mvxl2lph(:,2) = -1*ctx_vol_out.Mvxl2lph(:,2);
end
if Acorrvec(Adim)<0
  ctx_vol_out.imgs = flipdim(ctx_vol_out.imgs,3);  
  ctx_vol_out.Mvxl2lph(:,3) = -1*ctx_vol_out.Mvxl2lph(:,3);
end

% compensate for 1 voxel shift by mri_convert
ctx_vol_out.imgs(:,:,2:end) = ctx_vol_out.imgs(:,:,1:end-1);
ctx_vol_out.imgs(:,:,1) = 0;

if isempty(trans_vec)
  %ctr_lph2 = ctx_vol_out.Mvxl2lph(1:3,:)*[(volsz+1)/2 1]';
  %ctr_lph2 = ctx_vol_out.Mvxl2lph(1:3,:)*[(volsz-1)/2 1]';
  ctr_lph2 = ctx_vol_out.Mvxl2lph(1:3,:)*[(volsz)/2 1]';
  ctx_vol_out.Mvxl2lph(1:3,4) = ctx_vol_out.Mvxl2lph(1:3,4) - (ctr_lph2-ctr_lph);
else
  ctx_vol_out.Mvxl2lph(1:3,4) = trans_vec;
end;
