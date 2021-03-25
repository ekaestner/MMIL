% early mod: 03/09/12 by Don Hagler
% last mod:  04/12/12 by Don Hagler

dir_atlas = [getenv('MMPS_DIR') '/atlases'];
fname_T1_atlas = sprintf('%s/BOLD_reg/BOLD_reg_targ_T1.mgh',dir_atlas);
fname_BOLD_atlas = sprintf('%s/BOLD_reg/BOLD_reg_targ_BOLD.mgh',dir_atlas);
fname_BOLD_mask_atlas = sprintf('%s/BOLD_reg/BOLD_reg_targ_BOLD_mask.mgh',dir_atlas);
fname_nvals_BOLD = sprintf('%s/BOLD_reg/BOLD_nvals.mgh',dir_atlas);
fname_mask_atlas = sprintf('%s/BOLD_reg/BOLD_reg_targ_mask.mgh',dir_atlas);
stype = 'FSNU';
ftype = 'BOLD';
fname_jpdf = sprintf('%s/jpdfs/%s_%s_jpdf.mat',dir_atlas,stype,ftype);
dir_tmp = sprintf('%s/tmp_init_jpdf_FSNU_BOLD',dir_atlas);

forceflag = 0;
ndrop = 255;
numbins1 = 128;
numbins2 = 128;
hatl_offset_RAS = [0 -15 -10]; % displace resampled image from atlas (preserve nose)
M_hatl_offset = eye(4);
M_hatl_offset(1:3,4) = hatl_offset_RAS;
smf = 10^-5;

plotflag = 0;

FSRootDir = '/space/md3/2/halgdev/projects/megret/subjects';
PROCRootDir = '/space/md3/2/halgdev/projects/megret/fspace';

subject_list = {...
'amyj'...
'cletea'...
'ellas'...
'jdm'...
'larap'...
'taos'...
};

%'donh'...


func_list = {...
'060312AJ'...
'051126CA'...
'060528ES'...
'060219JM'...
'060312LP'...
'060129TS'...
};

%'060101DH'...

func_subdir_list = {...
'1-pol-dfix-ccw'...
'1-bw-check-pol-ccw'...
'1-pol-dfix-ccw'...
'1-pol-dfix-ccw'...
'1-pol-dfix-ccw'...
'1-bw-check-pol-ccw'...
};

%'2-bw-check-pol-ccw'...

func_fname_list = {...
'pol-dfix-ccw-vreg-unwarp_fmri+orig.BRIK'...
'bw-check-pol-ccw-strip-vreg-unwarp+orig.BRIK'...
'pol-dfix-ccw-vreg-unwarp_fmri+orig.BRIK'...
'pol-dfix-ccw-vreg-unwarp+orig.BRIK'...
'pol-dfix-ccw-vreg-unwarp_fmri+orig.BRIK'...
'bw-check-pol-ccw-vreg-unwarp+orig.BRIK'...
};

%'bw-check-pol-ccw-vreg-unwarp+orig.BRIK'...

volsz_atlas = [256,256,256];
M_atlas =  [-1     0     0   129;...
             0     0     1  -129;...
             0    -1     0   129;...
             0     0     0     1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsubs = length(subject_list);
if ~nsubs, return; end;

if exist(fname_jpdf,'file')
  fname_jpdf_bak = sprintf('%s/%s_%s_jpdf_bak_%s.mat',...
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

[status,msg] = mkdir(dir_tmp);
if ~status
  fprintf('%s: ERROR: failed to create tmp dir %s:\n',mfilename,dir_tmp);
  disp(msg);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input files
BOLD_fnames = [];
T1_fnames = [];
reg_fnames = [];
for i=1:nsubs
  % T1-weighted freesurfer volume
  fname_T1 = sprintf('%s/%s/mri/nu.mgz',...
    FSRootDir,subject_list{i});
  if ~exist(fname_T1,'file')
    fprintf('%s: ERROR: T1 file %s not found\n',mfilename,fname_T1);
    return;
  end;

  % BOLD volume registered to freesurfer
  fname_BOLD = sprintf('%s/%s/image/%s/%s',...
    PROCRootDir,func_list{i},func_subdir_list{i},func_fname_list{i});
  if ~exist(fname_BOLD,'file')
    fprintf('%s: ERROR: BOLD file %s not found\n',mfilename,fname_BOLD);
    return;
  end;
  
  % registration matrix file
  fname_reg = sprintf('%s/%s/image/%s/register.dat',...
    PROCRootDir,func_list{i},func_subdir_list{i});
  if ~exist(fname_reg,'file')
    fprintf('%s: ERROR: registration file %s not found\n',mfilename,fname_reg);
    return;
  end;

  BOLD_fnames{i} = fname_BOLD;
  T1_fnames{i} = fname_T1;
  reg_fnames{i} = fname_reg;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_avg_T1 = ctx_mgh2ctx(zeros(volsz_atlas),M_atlas);
vol_avg_BOLD = vol_avg_T1;
vol_avg_mask = vol_avg_T1;
nvals_T1 = 0;
nvals_BOLD = 0;
jpdf_mask_sum = 0;
for i=1:nsubs
  SubjID = subject_list{i};
  fprintf('%s: processing subject %s...\n',mfilename,SubjID);

  fname_BOLD = BOLD_fnames{i};
  fname_T1 = T1_fnames{i};
  fname_reg = reg_fnames{i};

  % convert BOLD from BRIK to mgh
  fname_BOLD_tmp = sprintf('%s/%s_BOLD.mgh',dir_tmp,SubjID);
  if ~exist(fname_BOLD_tmp,'file') || forceflag
    fprintf('%s:   converting input BOLD BRIK to mgh...\n',mfilename);
    cmd = sprintf('mri_convert --ndrop %d %s %s',...
      ndrop,fname_BOLD,fname_BOLD_tmp);
    [status,result] = unix(cmd);
    if status
      fprintf('%s: ERROR: failed to convert file %s:\n',mfilename,fname_BOLD);
      disp(result);
      return;
    end;
  end;
  fname_BOLD = fname_BOLD_tmp;

  % convert T1 from mgz to mgh
  fname_T1_tmp = sprintf('%s/%s_T1.mgh',dir_tmp,SubjID);
  if ~exist(fname_T1_tmp,'file') || forceflag
    fprintf('%s:   converting input T1 to mgh...\n',mfilename);
    cmd = sprintf('mri_convert %s %s',...
      fname_T1,fname_T1_tmp);
    [status,result] = unix(cmd);
    if status
      fprintf('%s: ERROR: failed to convert file %s:\n',mfilename,fname_T1);
      disp(result);
      return;
    end;
  end;
  fname_T1 = fname_T1_tmp;
  
  % make local copy of register.dat file
  fname_reg_tmp = sprintf('%s/%s_register.dat',dir_tmp,SubjID);
  if ~exist(fname_reg_tmp,'file') || forceflag
    cmd = sprintf('cp %s %s',fname_reg,fname_reg_tmp);
    [status,result] = unix(cmd);
    if status
      fprintf('%s: ERROR: failed to copy file %s:\n',mfilename,fname_reg);
      disp(result);
      return;
    end;
  end;
  fname_reg = fname_reg_tmp;

  % convert register.dat to fsl mat format
  fname_fslreg = sprintf('%s/%s_regfsl.mat',dir_tmp,SubjID);
  if ~exist(fname_fslreg,'file') || forceflag
    fprintf('%s:   converting register.dat file to fsl reg file...\n',mfilename);
    cmd = sprintf('tkregister2 --targ %s --mov %s',...
      fname_T1,fname_BOLD);
    cmd = sprintf('%s --reg %s --fslregout %s',...
      cmd,fname_reg,fname_fslreg);
    cmd = sprintf('%s --float2int floor --nofix --noedit',cmd);
    [status,result] = unix(cmd);
    if status
      fprintf('%s: cmd = %s\n',mfilename,cmd);
      fprintf('%s: ERROR: failed to generate fsl reg file:\n',mfilename);
      disp(result);
      return;
    end;
  end;

  % resample BOLD to T1 with fsl_rigid_register
  fname_BOLD_tmp = sprintf('%s/%s_BOLD_res_T1.mgh',dir_tmp,SubjID);
  if ~exist(fname_BOLD_tmp,'file')
    fprintf('%s:   resampling BOLD to T1...\n',mfilename);
    cmd = sprintf('fsl_rigid_register -r %s -i %s -o %s -applyxfm %s',...
      fname_T1,fname_BOLD,fname_BOLD_tmp,fname_fslreg);
    [status,result] = unix(cmd);
    if status
      fprintf('%s: ERROR: failed to generate fsl reg file:\n',mfilename);
      disp(result);
      return;
    end;
  end;
  fname_BOLD = fname_BOLD_tmp;

  % register T1 to atlas, create brain mask
  fname_T1_tmp = sprintf('%s/%s_T1_atlas.mgh',dir_tmp,SubjID);
  fname_mask = sprintf('%s/%s_brainmask_atlas.mgh',dir_tmp,SubjID);
  fname_reg_tmp = sprintf('%s/%s_reg_atlas.mat',dir_tmp,SubjID);
  if ~exist(fname_T1_tmp,'file') || ~exist(fname_mask,'file') |...
     ~exist(fname_reg_tmp,'file') || forceflag
    fprintf('%s:   registering to atlas...\n',mfilename);
    [vol_T1,mr_parms] = ctx_load_mgh(fname_T1);
    vol_mask = mmil_dct_brainmask(vol_T1);
    % rigid body register to atlas and resample to 256^3, 1mm^3, coronal slices
    vol_T1 = ctx_resample_to_LIA(vol_T1,M_atlas_to_T1);
    vol_mask = ctx_resample_to_LIA(vol_mask,M_atlas_to_T1);
    ctx_save_mgh(vol_T1,fname_T1_tmp,mr_parms);
    ctx_save_mgh(vol_mask,fname_mask,mr_parms);
    save(fname_reg_tmp,'M_atlas_to_T1');
  else
    vol_T1 = ctx_load_mgh(fname_T1_tmp);
    vol_mask = ctx_load_mgh(fname_mask);
    load(fname_reg_tmp);
  end;
  fname_T1 = fname_T1_tmp;
  
  fname_BOLD_tmp = sprintf('%s/%s_BOLD_atlas.mgh',dir_tmp,SubjID);
  fname_BOLD_mask = sprintf('%s/%s_BOLD_mask_atlas.mgh',dir_tmp,SubjID);
  if ~exist(fname_BOLD_tmp,'file') || ...
     ~exist(fname_BOLD_mask,'file') || forceflag
    fprintf('%s:   loading BOLD files %s...\n',mfilename,fname_BOLD);
    [vol_BOLD,mr_parms] = ctx_load_mgh(fname_BOLD);

    fprintf('%s:   creating BOLD mask %s...\n',mfilename,fname_BOLD);
    vol_BOLD_mask = vol_BOLD;
    vol_BOLD_mask.imgs = 1.0*(vol_BOLD_mask.imgs>smf);

    fprintf('%s:   resampling BOLD mask to atlas...\n',mfilename);
    vol_BOLD_mask = ctx_resample_to_LIA(vol_BOLD_mask,M_atlas_to_T1);
    ctx_save_mgh(vol_BOLD_mask,fname_BOLD_mask,mr_parms);

    fprintf('%s:   resampling BOLD to atlas...\n',mfilename);
    vol_BOLD = ctx_resample_to_LIA(vol_BOLD,M_atlas_to_T1);
    ctx_save_mgh(vol_BOLD,fname_BOLD_tmp,mr_parms);
  else
    vol_BOLD = ctx_load_mgh(fname_BOLD_tmp);
    vol_BOLD_mask = ctx_load_mgh(fname_BOLD_mask);
  end;
  fname_BOLD = fname_BOLD_tmp;

  fname_jpdf_tmp = sprintf('%s/%s_jpdf.mat',dir_tmp,SubjID);
  if ~exist(fname_jpdf_tmp,'file') || forceflag
    fprintf('%s:   combining masks...\n',mfilename);
    vol_comb_mask = vol_mask;
    vol_comb_mask.imgs = 1.0*(vol_comb_mask.imgs>smf);
    vol_comb_mask.imgs = vol_BOLD_mask.imgs .* vol_comb_mask.imgs;

    fprintf('%s:   calculating jpdf...\n',mfilename);
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
      fprintf('%s:     sf_T1 = %f\n',mfilename,sf_T1);
    end

    % compute scalefactor for vol_BOLD
    [hc,bv] = hist(vol_BOLD.imgs(:),1000);
    if i==1
      sf_BOLD = 1;
      v2_hc = hc;
      v2_bv = bv;
      v2_bn0 = 50;
    else
      sf_BOLD = ComputeHistScalefactor(hc,bv,v2_hc,v2_bv,v2_bn0);
      vol_BOLD.imgs = sf_BOLD*vol_BOLD.imgs;
      fprintf('%s:     sf_BOLD = %f\n',mfilename,sf_BOLD);
    end

    % calculate jpdf
    [sum_log10_val, jpdf_mask, bins1_mask, bins2_mask, jentropy]=...
      vols_jhist_mask_amd(vol_T1,vol_BOLD,vol_comb_mask,1,eye(4),eye(4),numbins1,numbins2);
    jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function

    save(fname_jpdf_tmp,'jpdf_mask','bins1_mask','bins2_mask',...
      'v2_hc','v2_bv','v2_bn0','v1_hc','v1_bv','v1_bn0');
  else
    fprintf('%s:   loading jpdf file %s...\n',mfilename,fname_jpdf_tmp);
    load(fname_jpdf_tmp);
  end;

  if plotflag
    imagesc(jpdf_mask);
    title(sprintf('jpdf for %s',SubjID));
    drawnow;
    pause(0.1);
    fname_plot = sprintf('%s/%s_%s_%s_jpdf_plot.jpg',...
      dir_tmp,SubjID,stype,ftype);
    print(gcf,'-djpeg',fname_plot);
  end;

  fprintf('%s:   adding to cross-subject sum...\n',mfilename);

  % track number of values in each voxel
  nvals_BOLD = nvals_BOLD + vol_BOLD_mask.imgs;

  % sum jpdf's across subjects
  jpdf_mask_sum = jpdf_mask_sum + jpdf_mask;

  vol_avg_BOLD.imgs = vol_avg_BOLD.imgs + vol_BOLD.imgs;
  vol_avg_T1.imgs = vol_avg_T1.imgs + vol_T1.imgs;
  vol_avg_mask.imgs = vol_avg_mask.imgs + vol_mask.imgs;
end;

fprintf('%s: averaging across subjects...\n',mfilename);
% average across subjects
jpdf_mask = jpdf_mask_sum / nsubs;
vol_avg_T1.imgs = vol_avg_T1.imgs / nsubs;
vol_avg_mask.imgs = vol_avg_mask.imgs / nsubs;
vol_avg_mask.imgs = 1.0*(vol_avg_mask.imgs>0.5);
vol_avg_BOLD.imgs = vol_avg_BOLD.imgs / nsubs;
vol_avg_BOLD.imgs(nvals_BOLD<nsubs-0.5) = 0;
vol_avg_BOLD_mask = vol_avg_mask;
vol_avg_BOLD_mask.imgs = 1.0*(nvals_BOLD>nsubs-0.5).*vol_avg_mask.imgs;

vol_nvals_BOLD = vol_avg_BOLD;
vol_nvals_BOLD.imgs = nvals_BOLD;

if plotflag
  fprintf('%s: plotting average jpdf...\n',mfilename);
  imagesc(jpdf_mask);
  title('average jpdf');
  drawnow;
  pause(0.1);
  fname_plot = sprintf('%s/%s_%s_jpdf_plot.jpg',dir_atlas,stype,ftype);
  print(gcf,'-djpeg',fname_plot);
end;

fprintf('%s: saving results...\n',mfilename);
ctx_save_mgh(vol_avg_T1,fname_T1_atlas);
ctx_save_mgh(vol_avg_BOLD,fname_BOLD_atlas);
ctx_save_mgh(vol_avg_BOLD_mask,fname_BOLD_mask_atlas);
%ctx_save_mgh(vol_nvals_BOLD,fname_nvals_BOLD);
ctx_save_mgh(vol_avg_mask,fname_mask_atlas);
save(fname_jpdf,'jpdf_mask','bins1_mask','bins2_mask',...
  'v2_hc','v2_bv','v2_bn0','v1_hc','v1_bv','v1_bn0');

