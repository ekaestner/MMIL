function [vol_corr,vol_ratio,vol_pred] = ...
  mmil_corr_wm_bias_spsm(vol,vol_wm,vol_bm,varargin)
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
%     relative to median ratio within brain mask
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
%     {default = 5}
%   'niters_spsm': number of iterations of sparse smoothing
%     {default = 2}
%
% Output:
%   vol_corr: corrected volume
%   vol_ratio: bias field ratio volume
%   vol_pre: predicted intensity volume
%
% Created:  02/23/15 by Don Hagler
% Prev Mod: 11/28/15 by Don Hagler
% Last Mod: 08/02/17 by Don Hagler
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
  'lambda2',5,[0,1e10],...
  'niters_spsm',3,[1,10],...
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

vol_wm_edge = abs(del2(vol_wm)).*vol_wm>0;
vol_wm_masked = vol_wm.*(~vol_wm_edge);

vol_tmp = FastSmoothVolume(vol(i1min:i1max,i2min:i2max,i3min:i3max),...
                           parms.lambda0,parms.lambda1,parms.lambda2,....
                           vol_wm_masked(i1min:i1max,i2min:i2max,i3min:i3max),...
                           parms.niters_spsm);
vol_pred = nan(size(vol));
vol_pred(i1min:i1max,i2min:i2max,i3min:i3max) = vol_tmp;

% calculate bias field ratio
vol_ratio = parms.vol_target./vol_pred;
med_ratio = median(vol_ratio(vol_bm>0));
vol_ratio(~vol_bm) = med_ratio;

% restrict range relative to median
if ~isempty(parms.ratio_range)
  vol_ratio(vol_ratio<med_ratio*parms.ratio_range(1)) = ...
    med_ratio*parms.ratio_range(1);
  vol_ratio(vol_ratio>med_ratio*parms.ratio_range(2)) = ...
    med_ratio*parms.ratio_range(2);
end;

% apply correction
vol_corr = vol.*vol_ratio;
vol_corr(isnan(vol_corr)) = 0;
vol_pred(isnan(vol_pred)) = 0;

if ~isempty(parms.image_range)
  vol_corr(vol_corr<parms.image_range(1)) = parms.image_range(1);
  vol_corr(vol_corr>parms.image_range(2)) = parms.image_range(2);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function volsm = FastSmoothVolume(vol,lambda0,lambda1,lambda2,wvol,niters)
  if ~exist('wvol','var')
    wvol = ones(size(vol));
  end
  dims = size(vol);
  c1 = round(dims(1)/2);
  c2 = round(dims(2)/2);
  c3 = round(dims(3)/2);
  volsm = vol;
  diffvol = volsm-vol;
  for iter = 0:niters
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


