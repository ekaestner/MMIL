function [vol_corr,vol_ratio,vol_pred] = ...
  mmil_corr_wm_bias_dct(vol,vol_wm,vol_bm,varargin)
%function [vol_corr,vol_ratio,vol_pred] = ...
%  mmil_corr_wm_bias_dct(vol,vol_wm,vol_bm,[options])
%
% Purpose: estimate B1 bias field using white matter segmentation mask
%   and Discrete Cosine Transform basis functions
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
%   'ndct': number of DCT basis functions
%     {default = 4}
%   'ratio_range': range of values to constrain ratio
%     if empty, no clipping
%     {default = []}
%   'image_range': range of values to clip output image
%     if empty, no clipping
%     {default = []}
%
% Output:
%   vol_corr: corrected volume
%   vol_ratio: bias field ratio volume
%   vol_pre: predicted intensity volume
%
% Created:  02/23/15 by Don Hagler
% Last Mod: 11/28/15 by Don Hagler
%
% Based on code by Anders Dale and Nate White
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_corr = [];
vol_ratio = [];
vol_pred = [];
if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms(varargin,{...
  'vol_target',[],[],...
  'ndct',4,[],...
  'ratio_range',[],[],...
  'image_range',[],[],...
});

vol_wm = 1.0*(vol_wm>0);
vol_bm = 1.0*(vol_bm>0);

if isempty(parms.vol_target)
  parms.vol_target = median(vol(vol_wm>0));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dims = size(vol);
indvec1 = [1:dims(1)]; indvec2 = [1:dims(2)]; indvec3 = [1:dims(3)];
[C1 C2 C3] = ...
  ndgrid(indvec1-mean(indvec1),indvec2-mean(indvec2),indvec3-mean(indvec3));

v1 = mmil_rowvec(sum(sum(vol_bm,3),2));
ivec1 = find(v1>0);
v2 = mmil_rowvec(sum(sum(vol_bm,3),1));
ivec2 = find(v2>0);
v3 = mmil_rowvec(sum(sum(vol_bm,2),1));
ivec3 = find(v3>0);

T = dctmtx(parms.ndct)';
T1 = nan(prod(dims),parms.ndct);
T2 = nan(prod(dims),parms.ndct);
T3 = nan(prod(dims),parms.ndct); 
for i = 1:parms.ndct
  v = interp1(linspace(ivec1(1),ivec1(end),parms.ndct),T(:,i),indvec1,'pchip','extrap');
  T1(:,i) = mmil_colvec(repmat(reshape(v,[dims(1) 1 1]),[1 dims(2) dims(3)]));
  v = interp1(linspace(ivec2(1),ivec2(end),parms.ndct),T(:,i),indvec2,'pchip','extrap');
  T2(:,i) = mmil_colvec(repmat(reshape(v,[1 dims(2) 1]),[dims(1) 1 dims(3)]));
  v = interp1(linspace(ivec3(1),ivec3(end),parms.ndct),T(:,i),indvec3,'pchip','extrap');
  T3(:,i) = mmil_colvec(repmat(reshape(v,[1 1 dims(3)]),[dims(1) dims(2) 1]));
end

A = nan(prod(dims),parms.ndct^3);
cnt = 0;
for i1 = 1:parms.ndct
  for i2 = 1:parms.ndct
    for i3 = 1:parms.ndct
      cnt = cnt+1;
      A(:,cnt) = parms.vol_target(:).*T1(:,i1).*T2(:,i2).*T3(:,i3);
    end;
  end;
end;

defvec = (vol_wm>0);
yvec = vol(:);
betahat = A(defvec,:)\yvec(defvec);
vol_pred = reshape(A*betahat,size(vol));

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

