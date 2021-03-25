function abcd_reg2mni_taskBOLD(cpath,fspath,varargin)
%function abcd_reg2mni_taskBOLD(cpath,fspath,[options])
%
% Required Input:
%   cpath: full path of proc_bold container
%   fspath: full path of fsurf container
%
% Optional Parameters:
%   'tasknames': task name or cell array of task names
%       allowed: {'MID','SST','nBack'}
%     {default = {'MID','SST','nBack'}}
%   'infix': string attached to processing BOLD data files
%     {default = 'corr_resBOLD'}
%   'concat_flag': [0|1|2] analyze concatenated across scans
%     0: analyze each scan individually
%     1: analyze concatenated scans
%     2: analyze individually and concatenated
%     {default = 1}
%   'outdir': output directory
%     {default = [pwd '/reg2mni_taskBOLD']}
%   'verbose': [0|1] display status updates
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  11/16/17 by Don Hagler
% Last Mod: 11/16/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'cpath',cpath,[],...
  'fspath',fspath,[],...
...
  'tasknames',{'MID','SST','nBack'},{'MID','SST','nBack'},...
  'infix','corr_resBOLD',[],...
  'concat_flag',1,[0:2],...
  'outdir',[pwd '/reg2mni_taskBOLD'],[],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
...
  'fext','.mgh',[],...
  'orient','RAS',[],...
  'brain_thresh',1,[],...
  'smooth_fwhm',4,[],...
  'reg_infix','reg2mni',[],...
});

if ~iscell(parms.tasknames)
  parms.tasknames = {parms.tasknames};
end;
parms.ntasks = length(parms.tasknames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load RegInfo
[RegInfo,fname_reg,errcode] = BOLD_MMIL_Load_RegInfo(parms.cpath,...
  'infix',parms.infix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

% resample head to fMRI-space
fname_head = sprintf('%s/head.nii',parms.outdir);
if ~exist(fname_head,'file') || parms.forceflag
  fname_orig = sprintf('%s/mri/nu.mgz',fspath);
  if parms.verbose
    fprintf('%s: resampling %s to fMRI space, saving as %s...\n',...
      mfilename,fname_orig,fname_head);
  end;
  [vol,M] = fs_load_mgh(fname_orig);
  [vol_res,M_res] = mmil_resample_vol(vol,M,...
    'M_ref',RegInfo.M_T2,'nvox_ref',RegInfo.volsz_T2(1:3),...
    'M_reg',inv(RegInfo.M_T1_to_T2));
  fs_save_nifti(vol_res,fname_head,M_res,[],parms.orient);
end;

% resample brain to fMRI-space
fname_brain_orig = sprintf('%s/brain.mgh',parms.outdir);
fname_brain = sprintf('%s/brain.nii',parms.outdir);
if ~exist(fname_brain,'file') || ~exist(fname_brain_orig,'file') || parms.forceflag
  if parms.verbose
    fprintf('%s: resampling brain from aseg to fMRI space, saving as %s...\n',...
      mfilename,fname_brain);
  end;
  fname_brain_tmp = sprintf('%s.mgh',mmil_tempfname('brain',parms.outdir));
  mmil_aseg_brainmask(fspath,'fname_mask',fname_brain_tmp,'forceflag',parms.forceflag);
  [vol,M] = fs_load_mgh(fname_brain_tmp);
  [vol_res,M_res] = mmil_resample_vol(vol,M,...
    'M_ref',RegInfo.M_T2,'nvox_ref',RegInfo.volsz_T2(1:3),...
    'M_reg',inv(RegInfo.M_T1_to_T2));
  fs_save_mgh(vol_res,fname_brain_orig,M_res);
  fs_mri_convert(fname_brain_orig,fname_brain,...
    'out_orient',parms.orient,'forceflag',parms.forceflag);
  delete(fname_brain_tmp);
end;
  
if parms.verbose
  fprintf('%s: registering to mni space...\n',mfilename);
end;

% register T1 to atlas with flirt and fnirt
tp = [];
tp.outdir = [parms.outdir '/' parms.reg_infix];
tp.verbose = parms.verbose;
tp.forceflag = parms.forceflag;
args = mmil_parms2args(tp);
mmil_reg2mni(fname_head,fname_brain,args{:});

% apply warp to T1
tp.fname_in = fname_head;
tp.fname_out = sprintf('%s/head_%s.nii',parms.outdir,parms.reg_infix);
args = mmil_parms2args(tp);
fname_head_atlas = mmil_reg2mni(fname_head,fname_brain,args{:});

% apply warp to mask
tp.fname_in = fname_brain;
tp.fname_out = sprintf('%s/brain_%s.nii',parms.outdir,parms.reg_infix);
args = mmil_parms2args(tp);
fname_brain_atlas = mmil_reg2mni(fname_head,fname_brain,args{:});

for j=1:parms.ntasks
  taskname = parms.tasknames{j};

  % get condition and contrast names
  switch taskname
    case 'MID'
      [stim_labels,conds_contrast] = abcd_set_contrasts_mid;
    case 'SST'
      [stim_labels,conds_contrast] = abcd_set_contrasts_sst;
    case 'nBack'
      [stim_labels,conds_contrast] = abcd_set_contrasts_nback;
  end;
  nstims = length(stim_labels);
  ncontrasts = length(conds_contrast);
  contrast_names = {conds_contrast.name};

  % find analysis subdirectories in proc_bold container
  switch parms.concat_flag
    case 0
      dpat = sprintf('%s/taskBOLD_%s_scan_*_%s_analysis',...
        parms.cpath,taskname,parms.infix);
    case 1
      dpat = sprintf('%s/taskBOLD_%s_scans_*_%s_analysis',...
        parms.cpath,taskname,parms.infix);
    case 2
      dpat = sprintf('%s/taskBOLD_%s_scan*_%s_analysis',...
        parms.cpath,taskname,parms.infix);
  end;
  dlist = dir(dpat);
  ndirs = length(dlist);
  for k=1:ndirs
    taskdir = dlist(k).name;
    outdir = [parms.outdir '/' taskdir];
    flist = dir(sprintf('%s/%s/taskBOLD_%s_scan*_%s_3dDeconv%s',...
      parms.cpath,taskdir,taskname,parms.infix,parms.fext));
    if isempty(flist)      
      fprintf('%s: WARNING: no 3dDeconv results files in %s\n',...
        mfilename,taskdir);
      continue;
    elseif length(flist)>1
      fprintf('%s: WARNING: multiple 3dDeconv results files in %s\n',...
        mfilename,taskdir);
      flist = flist(1);
    end;
    fname_in = sprintf('%s/%s/%s',parms.cpath,taskdir,flist.name);

    % split deconv results
    if parms.verbose
      fprintf('%s: splitting deconv results for %s...\n',mfilename,taskdir);
    end;
    vol = [];
    vol_brain = [];
    for i=1:ncontrasts
      cname = contrast_names{i};
      fstem = sprintf('%s/%s_%s',outdir,taskname,cname);
      fname_out = sprintf('%s.nii',fstem);
      if ~exist(fname_out,'file') || parms.forceflag
        if isempty(vol)
          ind_coef = 1 + 2*nstims + [1:2:2*ncontrasts];
          [vol,M,mrp,volsz] = fs_load_mgh(fname_in);
        end;
        vol_tmp = vol(:,:,:,ind_coef(i));
        % apply brain mask
        if isempty(vol_brain)
          vol_brain = fs_load_mgh(fname_brain_orig);
        end;
        vol_tmp = vol_tmp .* (vol_brain >= parms.brain_thresh);
        fs_save_nifti(vol_tmp,fname_out,M,[],parms.orient,...
          sprintf('--fwhm %d',parms.smooth_fwhm));
      end;
      % apply warp to coef map
      tp.fname_in = fname_out;
      tp.fname_out = sprintf('%s_%s.nii',fstem,parms.reg_infix);
      args = mmil_parms2args(tp);
      fname_reg = mmil_reg2mni(fname_head,fname_brain,args{:});
    end;
  end;
end;


