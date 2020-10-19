function [ctx_volrat,mr_parms] = mmil_estimate_B1RxRatio_GECal(indir,fstem_cal,varargin)
%function [ctx_volrat,mr_parms] = mmil_estimate_B1RxRatio_GECal(indir,fstem_cal,[options])...
%
% Required Input:
%  indir: directory containing mgh/mgz files
%  fstem_cal: file stem of calibration scan
%
% Optional Input:
%  'fstem_hires': file stem of high-res structural scan
%     if specified, ctx_volrat will be resampled
%     to orientation and resolution of structural scan
%    {default = []}
%  'stype': type of high-res structural scan
%        either 'MPR', 'FLASHhi', or 'FLASHlo'
%    {default = 'MPR'}
%  'regflag': [0|1] rigid body registration between
%     cal and hires scans using jpdf method
%    {default = 1}
%  'thresh': threshold for noise, multiplied by estimated noise level
%    {default = 4}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%    {default = '.mgz'}
%
% Output:
%  ctx_volrat: volume of ratio between head coil and body coil scans
%          (ctx format)
%
% Created:  03/09/12 by Don Hagler
% Rcnt Mod: 10/01/13 by Chun C. Fan
% Last Mod: 02/28/14 by Don Hagler
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ctx_volrat = []; mr_parms = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'indir',indir,[],...
  'fstem_cal',fstem_cal,[],...
  ...
  'fstem_hires',[],[],...
  'stype','MPR',{'MPR','FLASHhi','FLASHlo'},...
  'regflag',true,[false true],...
  'thresh',4,[],...
  'ext','.mgz',{'.mgh','.mgz'},...  
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ctx_volrat = [];

fname = sprintf('%s/%s%s',parms.indir,parms.fstem_cal,parms.ext);
if ~exist(fname,'file')
  error('file %s not found',fname);
end
[vol_cal,M_cal,mr_parms] = fs_load_mgh(fname);
nf = size(vol_cal,4);
nx = size(vol_cal,1);
smoothscale = round(nx/8);
if nf~=2
  fprintf('%s: WARNING: calibration file %s does not contain two frames\n',...
    mfilename,fname);
  return;
end;
vol_hc = vol_cal(:,:,:,1);
vol_bc = vol_cal(:,:,:,2);

if ~isempty(parms.fstem_hires)
  fname = sprintf('%s/%s.mgz',parms.indir,parms.fstem_hires);
  if ~exist(fname,'file')
    fname = sprintf('%s/%s.mgh',parms.indir,parms.fstem_hires);
    if ~exist(fname,'file')
      fprintf('%s: WARNING: high-res file %s not found\n',mfilename,fname);
      return;
    end;
  end
  [vol_hres,M_hres] = fs_load_mgh(fname);
  nx = size(vol_hres,1);
  smoothscale = round(nx/8);
  ctx_vol_hc = ctx_mgh2ctx(vol_hc,M_cal);
  ctx_vol_bc = ctx_mgh2ctx(vol_bc,M_cal);
  ctx_vol_hres = ctx_mgh2ctx(vol_hres,M_hres);
  if ~parms.regflag
    M_hres2cal = eye(4);  
  else
    fprintf('%s: registering calibration scan to %s (type %s)...\n',...
      mfilename,fname,parms.stype);
    M_hres2cal = mmil_atlas_jpdfreg_func2struct(ctx_vol_hc,ctx_vol_hres,...
        'ftype','GECAL','stype',parms.stype);
  end;
  ctx_vol_hc_res = vol_resample(ctx_vol_hc,ctx_vol_hres,M_hres2cal,1);
  ctx_vol_bc_res = vol_resample(ctx_vol_bc,ctx_vol_hres,M_hres2cal,1);
  [vol_hc,M_cal] = ctx_ctx2mgh(ctx_vol_hc_res);
  [vol_bc,M_cal] = ctx_ctx2mgh(ctx_vol_bc_res);
end;

if 0
  tmp1 = squeeze(vol_hres(:,:,90));
  tmp2 = squeeze(vol_hc(:,:,90));
  for i=1:100
    imagesc(tmp1); axis image;
    pause(0.1);  
    imagesc(tmp2); axis image;
    pause(0.1);  
  end;
end;

fprintf('%s: smoothing...\n',mfilename);
vol_hc = mmil_smooth3d(vol_hc,smoothscale,smoothscale,smoothscale);
vol_bc = mmil_smooth3d(vol_bc,smoothscale,smoothscale,smoothscale);
ctx_vol_hc = ctx_mgh2ctx(vol_hc,M_cal);
ctx_vol_bc = ctx_mgh2ctx(vol_bc,M_cal);

fprintf('%s: calculating ratio...\n',mfilename,fname);
[noiselevel_hc,signallevel_hc] = ...
  ctx_ComputeSignalAndNoiseLevels(ctx_vol_hc.imgs);
[noiselevel_bc,signallevel_bc] = ...
  ctx_ComputeSignalAndNoiseLevels(ctx_vol_bc.imgs);
maskindx = find(ctx_vol_hc.imgs>parms.thresh*noiselevel_hc);
mu_hc = mean(ctx_vol_hc.imgs(maskindx));
mu_bc = mean(ctx_vol_bc.imgs(maskindx));
murat = mu_bc/mu_hc;
ctx_volrat = ctx_vol_hc;
ctx_volrat.imgs = ...
  sqrt(abs((ctx_vol_bc.imgs.^2-noiselevel_bc^2)./...
  (ctx_vol_hc.imgs.^2-noiselevel_hc^2)));
voltmp = max(0.5*murat,min(2*murat,ctx_volrat.imgs));
ctx_volrat.imgs = zeros(size(ctx_volrat.imgs));
ctx_volrat.imgs(maskindx) = voltmp(maskindx);

