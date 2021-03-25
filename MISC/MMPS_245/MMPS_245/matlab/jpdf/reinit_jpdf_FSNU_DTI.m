function reinit_jpdf_FSNU_DTI(ProjID)
%
% Purpose: creating multi-subject jpdf for T1 to DTI registration
%
% Early Mod:03/18/10 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

dir_atlas = [getenv('MMPS_DIR') 'atlases'];
fname_T1_atlas = sprintf('%s/DTI_reg/DTI_reg_targ_T1.mgh',dir_atlas);
fname_DTI_atlas = sprintf('%s/DTI_reg/DTI_reg_targ_DTI.mgh',dir_atlas);
stype = 'FSNU';
ftype = 'DTI';
fname_jpdf = sprintf('%s/jpdfs/%s_%s_jpdf.mat',dir_atlas,stype,ftype);

min_ndirs = 6;
min_bval = 1000;
flex_flag = 0;
min_nb0 = 1;

plotflag = 1;

numbins1 = 128;
numbins2 = 128;

revflag = 0;
partial_infix = 'corr_resDTI';

subject_list = {...
'ma01'...
'cm01'...
'ld01'...
'dk01'...
'kg01'...
'mg01'...
'tb01'...
'kc01'...
'jj01'...
'rt01'...
'dw01'...
'js01'...
'ls01'...
'lm01'...
'mm01'...
'me01'...
'vl01'...
};

volsz_atlas = [256,256,256];
M_atlas = [-1 0 0 129; 0 0 1 -129; 0 -1 0 129; 0 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsubs = length(subject_list);
if ~nsubs, return; end;

[StudyInfo,RootDirs,ProjInfo] = MMIL_Get_StudyInfo(ProjID);

if isempty(StudyInfo)
  fprintf('%s: ERROR reading StudyInfo\n',mfilename);
  return;
end;
SubjIDs = {StudyInfo.SubjID};

if exist(fname_jpdf,'file')
  fname_jpdf_bak = sprintf('%s/jpdfs/%s_%s_jpdf_bak_%s.mat',...
    dir_atlas,stype,ftype,datestr(date,'yymmdd'));
  cmd = sprintf('mv %s %s',fname_jpdf,fname_jpdf_bak);
  [status,result] = unix(cmd);
  if status
    fprintf('%s: ERROR: failed to backup jpdf file %s to %s:\n',...
      mfilename,fname_jpdf,fname_jpdf_bak);
    disp(result);
    return;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input files
DTI_fnames = [];
T1_fnames = [];
mask_fnames = [];
for i=1:nsubs
  SubjID = subject_list{i};
  j= find(strcmp(SubjID,SubjIDs));
  if isempty(j)
    fprintf('%s: ERROR: bad SubjID: %s\n',mfilename,SubjID);
    return;
  end;

  snums = StudyInfo(j).DTIScanNums;
  PROCContainerPath = sprintf('%s/%s',RootDirs.proc,StudyInfo(j).proc);
  FSContainerPath = sprintf('%s/%s',RootDirs.fsurf,StudyInfo(j).fsurf);

  % select infix
  infix = DTI_MMIL_Get_Infix(PROCContainerPath,snums,partial_infix);
  if isempty(infix)
    fprintf('%s: ERROR: failed to find infix for %s with %s\n',...
      mfilename,PROCContainerPath,partial_infix);
    return;
  end;

  % b=0 DTI volume registered to freesurfer
  DT_fstem = DTI_MMIL_Set_DT_fstem(PROCContainerPath,
    'snums',snums,'infix',infix,'revflag',revflag,...
    'min_ndirs',min_ndirs,'min_bval',min_bval,...
    'flex_flag',flex_flag,'min_nb0',min_nb0);
  fname_dti = sprintf('%s_b0_resT1.mgz',DT_fstem);
  if ~exist(fname_dti,'file')
    fprintf('%s: ERROR: DTI file %s not found\n',mfilename,fname_DTI);
    return;
  end;
  
  % T1-weighted freesurfer volume
  fname_T1 = sprintf('%s/mri/nu.mgz',FSContainerPath);
  if ~exist(fname_T1,'file')
    fprintf('%s: ERROR: T1 file %s not found\n',mfilename,fname_T1);
    return;
  end;

  % brain mask volume
  fname_mask = sprintf('%s/mri/nu-brainmask.mgz',FSContainerPath);
  if ~exist(fname_mask,'file')
    fprintf('%s: ERROR: mask file %s not found\n',mfilename,fname_mask);
    return;
  end;

  DTI_fnames{i} = fname_dti;
  T1_fnames{i} = fname_T1;
  mask_fnames{i} = fname_mask;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_avg_T1 = ctx_mgh2ctx(zeros(volsz_atlas),M_atlas);
vol_avg_DTI = vol_avg_T1;
jpdf_mask_sum = 0;
for i=1:nsubs
  SubjID = subject_list{i};
  fprintf('%s: loading data for %s...\n',mfilename,SubjID);
  vol_DTI = ctx_load_mgh(DTI_fnames{i});
  vol_T1 = ctx_load_mgh(T1_fnames{i});
  vol_mask = ctx_load_mgh(mask_fnames{i});

  % compute scalefactor for vol_T1
  [hc,bv] = hist(vol_T1.imgs(:),1000);
  if i==1
    sf_T1 = 1;
    v1_hc = hc;
    v1_bv = bv;
    v1_bn0 = 50;
  else
    sf_T1 = ComputeHistScalefactor(hc,bv,v1_hc,v1_bv,v1_bn0);
    vol_T1.imgs = sf_T1*vol_T1.imgs;
    fprintf('%s: sf_T1 = %f\n',mfilename,sf_T1);
  end

  % compute scalefactor for vol_DTI
  [hc,bv] = hist(vol_DTI.imgs(:),1000);
  if i==1
    sf_DTI = 1;
    v2_hc = hc;
    v2_bv = bv;
    v2_bn0 = 50;
  else
    sf_DTI = ComputeHistScalefactor(hc,bv,v2_hc,v2_bv,v2_bn0);
    vol_DTI.imgs = sf_DTI*vol_DTI.imgs;
    fprintf('%s: sf_DTI = %f\n',mfilename,sf_DTI);
  end

  % calculate jpdf
  [sum_log10_val, jpdf_mask, bins1_mask, bins2_mask, jentropy]=...
    vols_jhist_mask_amd(vol_T1,vol_DTI,vol_mask,1,eye(4),eye(4),numbins1,numbins2);
  jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function
  if plotflag
    imagesc(jpdf_mask);
    drawnow;
    pause(0.1);
  end;

  % average jpdf's across subjects
  jpdf_mask_sum = jpdf_mask_sum + jpdf_mask;

  % NOTE: T1 images are pre-registered to atlas by
  %   MMIL_Register_and_Resample_MRI_Volumes
  vol_avg_DTI.imgs = vol_avg_DTI.imgs + vol_DTI.imgs;
  vol_avg_T1.imgs = vol_avg_T1.imgs + vol_T1.imgs;
end;
jpdf_mask = jpdf_mask_sum / nsubs;

if plotflag
  fprintf('%s: plotting average jpdf...\n',mfilename);
  imagesc(jpdf_mask);
  drawnow;
  pause(0.1);
end;

fprintf('%s: saving results...\n',mfilename);
vol_avg_T1.imgs = vol_avg_T1.imgs / nsubs;
vol_avg_DTI.imgs = vol_avg_DTI.imgs / nsubs;
ctx_save_mgh(vol_avg_T1,fname_T1_atlas);
ctx_save_mgh(vol_avg_DTI,fname_DTI_atlas);
save(fname_jpdf,'jpdf_mask','bins1_mask','bins2_mask',...
  'v2_hc','v2_bv','v2_bn0','v1_hc','v1_bv','v1_bn0');

