function MMIL_Register_MRI_to_MRI(ContainerRootDir,ContainerDir1,ContainerDir2,prevflag,sparseflag,forceflag);
%function MMIL_Register_MRI_to_MRI(ContainerRootDir,ContainerDir1,ContainerDir2,[prevflag],[sparseflag],[forceflag]);
%
% prevflag: look for MPR_res_reg.mgz (register to previous scan that was registered to baseline)
%   default=0
% sparseflag: renormalize intensities ignoring voxels where normalization changes things too much
%   default=0
% forceflag: repeat calculations even if output files exist
%   default=0
%
%
% Early Mod: 01/24/07 by Don Hagler
% Rcnt Mod:  05/21/12 by Don Hagler
% Last Mod:  09/14/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

if ~exist('prevflag','var'), prevflag = []; end
if isempty(prevflag), prevflag=0; end;
if ~exist('sparseflag','var'), sparseflag = []; end
if isempty(sparseflag), sparseflag=0; end;
if ~exist('forceflag','var'), forceflag = 0; end

ext = '.mgz';
diff_thresh = 0.75;
spsm_niter = 50;

fprintf('%s(''%s'',''%s'',''%s'',%d)\n',...
  mfilename,ContainerRootDir,ContainerDir1,ContainerDir2,forceflag);

fullCDir1 = sprintf('%s/%s',ContainerRootDir,ContainerDir1);
fullCDir2 = sprintf('%s/%s',ContainerRootDir,ContainerDir2);

if ~prevflag
  fname1 = sprintf('%s/MPR_res%s',fullCDir1,ext);
else
  fname1 = sprintf('%s/MPR_res_reg%s',fullCDir1,ext);
end;
fname2 = sprintf('%s/MPR_res%s',fullCDir2,ext);
fname_out = sprintf('%s/MPR_res_reg%s',fullCDir2,ext);

if ~exist(fname_out) | forceflag
  if ~exist(fname1,'file')
    fprintf('%s: file %s not found\n',mfilename,fname1);
  end
  if ~exist(fname2,'file')
    fprintf('%s: file %s not found\n',mfilename,fname2);
  end
  try
    vol1 = ctx_load_mgh(fname1);
    vol2 = ctx_load_mgh(fname2);
  catch
    fprintf(1,'%s: error loading files\n',mfilename);
    return;
  end
  
  % find brain mask for target volume
  volmask = mmil_dct_brainmask(vol1);

  % register volumes
  [M_v1_to_v2, min_cost] = rbreg_vol2vol_icorr_mask(vol1, vol2, volmask, 0, 4);

  % normalize intensities
  vol2_res = vol_resample(vol2, vol1, M_v1_to_v2, 2);
  vol1_sm = vol1;
  vol2_res_sm = vol2_res;
  vol1_sm.imgs = mmil_smooth3d(vol1.imgs,32,32,32);
  vol2_res_sm.imgs = mmil_smooth3d(vol2_res.imgs,32,32,32);
  mean1 = mean(abs(vol1_sm.imgs(find(volmask.imgs(:)))));
  mean2 = mean(abs(vol2_res_sm.imgs(find(volmask.imgs(:)))));
  valreg1 = 1e-6*mean1; % regularization
  valreg2 = 1e-6*mean2;
  volrat = (abs(vol1_sm.imgs)+valreg1)./(abs(vol2_res_sm.imgs)+valreg2);
  vol2_res_norm = vol2_res;
  vol2_res_norm.imgs = vol2_res.imgs.*volrat;

  if sparseflag
    % calculate normalized difference in intensities
    tmp_vol2_res = vol2_res;
    tmp_vol2_res.imgs = vol2_res.imgs*(mean1/mean2+eps);
    vol_diff = sqrt((vol2_res_norm.imgs - tmp_vol2_res.imgs).^2)./(abs(tmp_vol2_res.imgs)+eps);
    vol_diff_mask = vol_diff;
    vol_diff_mask(find(vol_diff>diff_thresh))=0;
    vol_diff_mask(find(vol_diff<diff_thresh))=1;

    % renormalize intensities - ignore voxels where normalization changed things too much
    vol1_spsm = vol1;
    vol2_res_spsm = vol2_res;
    % sparse smoothing to fill in gaps (lesions)
    vol1_spsm.imgs = double(SparseSmoothScalarVol(single(vol1.imgs),vol_diff_mask,spsm_niter));
    vol2_res_spsm.imgs = double(SparseSmoothScalarVol(single(vol2_res.imgs),vol_diff_mask,spsm_niter));
    % isotropic smoothing for massive blurring
    vol1_sm.imgs = mmil_smooth3d(vol1_spsm.imgs,32,32,32);
    vol2_res_sm.imgs = mmil_smooth3d(vol2_res_spsm.imgs,32,32,32);
    % recalculate ratio
    volrat = (abs(vol1_sm.imgs)+valreg1)./(abs(vol2_res_sm.imgs)+valreg2);
    vol2_res_norm.imgs = vol2_res.imgs.*volrat;
  end;

  [vol,M] = ctx_ctx2mgh(vol2_res_norm);
  fs_save_mgh(vol,fname_out,M);
end;


