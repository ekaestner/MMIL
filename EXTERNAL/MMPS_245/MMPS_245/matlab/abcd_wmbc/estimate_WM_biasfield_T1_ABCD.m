function [vol_corr,vol_ratio,vol_pred] = estimate_WM_biasfield_T1_ABCD(vol,vol_wm,vol_bm,vol_gm,lambda2)

if ~exist('lambda2','var') | isempty(lambda2), lambda2 = 10; end;

%vol_bm.imgs(40:210,40:end,20:end) = 1; % Should use anatomical bounding box
vol_bm.imgs(:) = 1; % Include all voxels (may be slow)

vol_corr = vol; vol_ratio = vol; vol_pred = vol;
[vol_corr.imgs,vol_ratio.imgs,vol_pred.imgs] = mmil_corr_wm_bias_spsm_T1_ABCD(vol.imgs,vol_wm.imgs,vol_bm.imgs,vol_gm.imgs,'vol_target',110,'lambda2',lambda2);


