function MMIL_Calculate_T1Ratio(ContainerPath,forceflag);
%function MMIL_Calculate_T1Ratio(ContainerPath,forceflag);
%
%
% Created:   07/01/08 by Don Hagler
% Rcnt Mod:  03/09/12 by Don Hagler
% Last Mod:  09/14/12 by Don Hagler
%

if ~exist('forceflag','var'), forceflag = 0; end;
if ~exist('scaleflag','var'), scaleflag = 0; end;
ext = '.mgz';

%% todo: Should be based on average of larger #s of subjects
fname_hist2d = [getenv('MMPS_DIR') '/atlases/T1_Atlas/hist2d_target.mat'];

fname_hi_in = sprintf('%s/hiFA_res%s',ContainerPath,ext);
fname_lo_in = sprintf('%s/loFA_res%s',ContainerPath,ext);

if ~exist(fname_hi_in,'file') || ~exist(fname_lo_in,'file'), return; end;

if scaleflag
  fname_hi_out = sprintf('%s/hiFA_res_scaled%s',ContainerPath,ext);
  fname_lo_out = sprintf('%s/loFA_res_scaled%s',ContainerPath,ext);

  % calculate scale factors for hiFA and loFA images
  if ~exist(fname_hi_out,'file') & ~exist(fname_lo_out,'file') || forceflag
    fprintf('%s: computing scale factors for hiFA and loFA...\n',mfilename);
    load(fname_hist2d);
    [vol_hiFA,mr_parms_hiFA] = ctx_load_mgh(fname_hi_in);
    [vol_loFA,mr_parms_loFA] = ctx_load_mgh(fname_lo_in);
    [sf_hiFA,sf_loFA] = compute_scalefact_hist2d(...
      vol_hiFA.imgs(find(volmask.imgs)),...
      vol_loFA.imgs(find(volmask.imgs)),...
      hcnt_target,binvals1,binvals2);
    vol_hiFA.imgs = sf_hiFA*vol_hiFA.imgs;
    vol_loFA.imgs = sf_loFA*vol_loFA.imgs;
    fprintf('%s: saving %s...\n',mfilename,fname_hi_out);
    ctx_save_mgh(vol_hiFA,fname_hi_out,mr_parms_hiFA);
    fprintf('%s: saving %s...\n',mfilename,fname_lo_out);
    ctx_save_mgh(vol_loFA,fname_lo_out,mr_parms_loFA);
  end;

  fname_hi_in = fname_hi_out;
  fname_lo_in = fname_lo_out;

  fname_ratio = sprintf('%s/T1ratio_res_scaled%s',ContainerPath,ext);
else
  fname_ratio = sprintf('%s/T1ratio_res%s',ContainerPath,ext);
end;


if ~exist(fname_ratio,'file') || forceflag
  fprintf('%s: computing T1 ratio...\n',mfilename);
  [vol_hiFA,mr_parms_hiFA] = ctx_load_mgh(fname_hi_in);
  [vol_loFA,mr_parms_loFA] = ctx_load_mgh(fname_lo_in);
  vol_ratio = vol_hiFA;
  vol_ratio.imgs = vol_hiFA.imgs./(vol_hiFA.imgs+vol_loFA.imgs+eps);
  fprintf('%s: saving %s...\n',mfilename,fname_ratio);
  ctx_save_mgh(vol_ratio,fname_ratio,mr_parms_hiFA);
end;
 
