function ctx_volrat = mmil_estimate_B1RxRatio_ADNI(indir,BCmstem,BCpstem,...
                                                         HCmstem,HCpstem,ext)
%function ctx_volrat = mmil_estimate_B1RxRatio_ADNI(indir,BCmstem,BCpstem,...
%                                                         HCmstem,HCpstem,[ext])
%
% Created:  03/09/12 by Don Hagler
% Last Mod: 02/28/14 by Don Hagler
%

if ~mmil_check_nargs(nargin,5), return; end;

if ~exist('ext','var') || isempty(ext), ext = '.mgz'; end;

sm = 32;
thresh = 2;

fname1 = sprintf('%s/%s%s',indir,HCmstem,ext);
fname2 = sprintf('%s/%s%s',indir,HCpstem,ext);
if exist(fname1,'file')
  ctx_vol_abs = ctx_load_mgh(fname1);
  if exist(fname2,'file')
    ctx_vol_phi = ctx_load_mgh(fname2);
    vol = complex(ctx_vol_abs.imgs.*cos(ctx_vol_phi.imgs),ctx_vol_abs.imgs.*sin(ctx_vol_phi.imgs));
  else
    vol = ctx_vol_abs.imgs;
  end
else
  ctx_volrat = [];
  return;
end
vol = abs(mmil_smooth3d(vol,sm,sm,sm));
ctx_vol_hc = ctx_vol_abs;
ctx_vol_hc.imgs = vol;
clear('ctx_vol_abs','ctx_vol_phi','vol');

fname1 = sprintf('%s/%s%s',indir,BCmstem,ext);
fname2 = sprintf('%s/%s%s',indir,BCpstem,ext);
if exist(fname1,'file')
  ctx_vol_abs = ctx_load_mgh(fname1);
  if exist(fname2,'file')
    ctx_vol_phi = ctx_load_mgh(fname2);
    vol = complex(ctx_vol_abs.imgs.*cos(ctx_vol_phi.imgs),ctx_vol_abs.imgs.*sin(ctx_vol_phi.imgs));
  else
    vol = ctx_vol_abs.imgs;
  end
else
  ctx_volrat = [];
  return;
end
vol = abs(mmil_smooth3d(vol,sm,sm,sm));
ctx_vol_bc = ctx_vol_abs;
ctx_vol_bc.imgs = vol;
clear('ctx_vol_abs','ctx_vol_phi','vol');

[noiselevel_hc,signallevel_hc] = ...
  ctx_ComputeSignalAndNoiseLevels(ctx_vol_hc.imgs);
[noiselevel_bc,signallevel_bc] = ...
  ctx_ComputeSignalAndNoiseLevels(ctx_vol_bc.imgs);

fprintf('%s: registering B1BC and B1HC images...\n',mfilename);
[M_v1_to_v2, min_cost] = rbreg_vol2vol_jpdf(ctx_vol_bc,ctx_vol_hc);
ctx_vol_hc = vol_resample(ctx_vol_hc, ctx_vol_bc, M_v1_to_v2, 2);

maskindx = find(ctx_vol_hc.imgs>thresh*noiselevel_hc);
mu_hc = mean(ctx_vol_hc.imgs(maskindx));
mu_bc = mean(ctx_vol_bc.imgs(maskindx));
murat = mu_bc/mu_hc;

ctx_volrat = ctx_vol_hc;
ctx_volrat.imgs = sqrt(abs((ctx_vol_bc.imgs.^2-noiselevel_bc^2) ./ ...
                      (ctx_vol_hc.imgs.^2-noiselevel_hc^2)));
voltmp = max(0.5*murat,min(2*murat,ctx_volrat.imgs));
ctx_volrat.imgs = zeros(size(ctx_volrat.imgs));
ctx_volrat.imgs(maskindx) = voltmp(maskindx);

