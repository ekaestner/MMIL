function [vol_corr,vol_ratio,vol_pred] = ...
  mmil_corr_wm_bias_spsm_amd(vol,vol_wm,vol_bm,vol_gm,varargin)
%function [vol_corr,vol_ratio,vol_pred] = ...
%  mmil_corr_wm_bias_spsm(vol,vol_wm,vol_bm,[options])
%
% Purpose: estimate B1 bias field using white matter segmentation mask
%   and sparse smoothing in volume
%
% Required Input:
%   vol: brain image volume
%   vol_wm: white matter mask volume
%   vol_bm: brain mask volume
%
% Optional Input ('key',value pairs):
%   'vol_target': target volume for white matter intensity
%     may be a scalar value
%     if empty, will use median white matter intensity value
%     {default = []}
%   'ratio_range': range of values to constrain ratio
%     relative to medianaratio within brain mask
%     if empty, no clipping
%     {default = []}
%   'image_range': range of values to clip output image
%     if empty, no clipping
%     {default = []}
%   'lambda0': weighting factor for difference volume
%     {default = 1}
%   'lambda1': weighting factor for smooth volume
%     {default = 0}
%   'lambda2': weighting factor for Laplacian of smooth volume
%     {default = 100}
%
% Output:
%   vol_corr: corrected volume
%   vol_ratio: bias field ratio volume
%   vol_pre: predicted intensity volume
%
% Created:  02/23/15 by Don Hagler
% Prev Mod: 11/28/15 by Don Hagler
% Prev Mod: 11/06/17 by Anders Dale
% Last Mod: 11/06/17 by Don Hagler
%
% Based on code by Anders Dale
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_corr = [];
vol_ratio = [];
vol_pred = [];
if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms(varargin,{...
  'vol_target',[],[],...
  'ratio_range',[],[],...
  'image_range',[],[],...
  'lambda0',1,[0,1e10],...
  'lambda1',0,[0,1e10],...
  'lambda2',1e2,[0,1e10],...
});

vol_wm = 1.0*(vol_wm>0);
vol_bm = 1.0*(vol_bm>0);

if isempty(parms.vol_target)
  parms.vol_target = median(vol(vol_wm>0));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[i1max i1min] = maxmin(find(sum(sum(vol_bm,3),2)));
[i2max i2min] = maxmin(find(sum(sum(vol_bm,3),1)));
[i3max i3min] = maxmin(find(sum(sum(vol_bm,2),1)));

vol_val_wm = vol.*vol_wm; 

%showVol(ctx_mgh2ctx(vol,eye(4)),ctx_mgh2ctx(vol_wm,eye(4)),ctx_mgh2ctx(vol_val_wm,eye(4)));

% Remove edge voxels from mask if the intensity is much less than local maximnum within mask -- should look for a way to use imdilate 
vol_val_locmax = vol_val_wm;
niter = 5;
for iter = 1:niter
% vol_val_locmax = vol_wm.*max(cat(4,circshift(vol_val_locmax,-1,1),circshift(vol_val_locmax,1,1),circshift(vol_val_locmax,-1,2),circshift(vol_val_locmax,1,2),circshift(vol_val_locmax,-1,3),circshift(vol_val_locmax,1,3)),[],4);
 vol_val_locmax = vol_wm.*max(cat(4,circshift(vol_val_locmax,[-1,0,0]),circshift(vol_val_locmax,[1,0,0]),circshift(vol_val_locmax,[0,-1,0]),circshift(vol_val_locmax,[0,1,0]),circshift(vol_val_locmax,[0,0,-1]),circshift(vol_val_locmax,[0,0,1])),[],4);
end
vol_wm_masked = vol_wm.*(vol_val_wm./vol_val_locmax>0.8);

niter = 1; % Not clear that higher number of iterations improves things
for iter = 1:niter
  vol_tmp = FastSmoothVolume(vol(i1min:i1max,i2min:i2max,i3min:i3max),parms.lambda0,parms.lambda1,parms.lambda2,vol_wm_masked(i1min:i1max,i2min:i2max,i3min:i3max));
  vol_pred = nan(size(vol));
  vol_pred(i1min:i1max,i2min:i2max,i3min:i3max) = vol_tmp;
  vol_wm_masked = vol_wm_masked.*(vol./vol_pred>0.9);
end

% udpate vol_wm_masked based on vol_val_locmax << vol_pred, or vol << vol_pred, then run FastSmoothVolume(vol

% calculate bias field ratio
vol_ratio = parms.vol_target./vol_pred;

vol_ratio_sm = vol_ratio; vol_ratio_sm(~vol_wm_masked) = NaN;
vol_ratio_sm = smoothn(vol_ratio_sm,0.1);

% apply correction
vol_corr = vol.*vol_ratio_sm;
vol_corr(isnan(vol_corr)) = 0;
vol_pred(isnan(vol_pred)) = 0;

%keyboard

return;

showVol(ctx_mgh2ctx(vol,eye(4)),ctx_mgh2ctx(vol_ratio,eye(4)),ctx_mgh2ctx(vol_ratio_sm,eye(4)),ctx_mgh2ctx(vol_corr,eye(4)),ctx_mgh2ctx(vol_corr.*vol_wm_masked,eye(4)),ctx_mgh2ctx(vol_val_locmax.*vol_wm_masked,eye(4)),ctx_mgh2ctx(vol_pred.*vol_wm_masked,eye(4)));

% ToDo
%   Look into using smoothn instead of FastSmoothVolume (should be much faster)
%   Should include only voxels that have a high a posteriori probability of being WM (look inside corSeg)
%   Could use previous vol_tmp as initial estimate, to save time
%   Smoothness constraint should be on vol_ratio, not vol!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function volsm = FastSmoothVolume(vol,lambda0,lambda1,lambda2,wvol)
  if ~exist('wvol','var')
    wvol = ones(size(vol));
  end
  dims = size(vol);
  c1 = round(dims(1)/2);
  c2 = round(dims(2)/2);
  c3 = round(dims(3)/2);
  niter = 2;
  volsm = vol;
  diffvol = volsm-vol;
  for iter = 0:niter
    tic
    niters = nan;
    if iter > 0
      g = lambda0*wvol.*diffvol + lambda1*volsm + lambda2*volLaplacian(volsm);
      H = @ (x) SparseSmoothHfun(x,lambda0,lambda1,lambda2,wvol);
      [dx,flag,relres,niters,resvec] = bicgstab(H,-mmil_colvec(g),1e-2,1000);
      dx = reshape(dx,dims);
      volsm = volsm+dx;
    end
    diffvol = volsm-vol;
    costvol0 = wvol.*diffvol;
    costvol1 = sqrt(volsm.^2);
    costvol2 = sqrt(cat(1,diff(volsm,1,1),zeros(1,dims(2),dims(3))).^2 + cat(2,diff(volsm,1,2),zeros(dims(1),1,dims(3))).^2 + cat(3,diff(volsm,1,3),zeros(dims(1),dims(2),1)).^2);
    cost = 1/2*lambda0*sum(costvol0(:).^2) + 1/2*lambda1*sum(costvol1(:).^2) + 1/2*lambda2*sum(costvol2(:).^2);
    fprintf('%s: iter=%d: cost=%f, niters=%.1f\n',mfilename,iter,cost,niters)
    toc
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lapvol = volLaplacian(vol)
  dims = size(vol);
  lapvol = 0;
  D = sparse(diff(eye(dims(1))));
  DtD = D'*D;
  lapvol = lapvol + reshape(DtD*reshape(double(vol),[dims(1) dims(2)*dims(3)]),dims);
  D = sparse(diff(eye(dims(2))));
  DtD = D'*D;
  lapvol = lapvol + permute(reshape(DtD*reshape(permute(double(vol),[2 1 3]),[dims(2) dims(1)*dims(3)]),[dims(2) dims(1) dims(3)]),[2 1 3]);
  D = sparse(diff(eye(dims(3))));
  DtD = D'*D;
  lapvol = lapvol + permute(reshape(DtD*reshape(permute(double(vol),[3 2 1]),[dims(3) dims(2)*dims(1)]),[dims(3) dims(2) dims(1)]),[3 2 1]);
  % should have opposite sign??
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = SparseSmoothHfun(x,lambda0,lambda1,lambda2,w)
  dims = size(w);
  x = reshape(x,dims);
  y = lambda0*w.*x + lambda1*x + lambda2*volLaplacian(x);
  y = mmil_colvec(y);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


