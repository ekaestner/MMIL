function ABCD_Analyze_taskBOLD_Exam(ContainerPath,FSContainerPath,fname_eprime,varargin)
%function ABCD_Analyze_taskBOLD_Exam(ContainerPath,FSContainerPath,fname_eprime,[options])
%
% Usage:
%  ABCD_Analyze_taskBOLD_Exam(ContainerPath,FSContainerPath,fname_eprime,'key1', value1,...)
%
% Required Parameters:
%   ContainerPath: full path of BOLDPROC directory containing processed BOLD data
%   FSContainerPath: full path of directory containing freesurfer recon
%   fname_eprime: full path of csv or txt file derived from eprime output
%
% Optional GLM Analysis Parameters:
%   'outstem': output file stem used for output directories
%     if empty, will construct from 'BOLD' and taskname
%     {default = []}
%   'analysis_outfix': string attached to analysis output subdirectories
%     {default = 'analysis'}
%   'snums': vector of scan numbers on which to run GLM analysis
%     if empty, will skip GLM analysis
%     {default = []}
%   'stim_times_flag': [0|1|2] use stimulus time files
%     0 = .1D, 1 = .txt and 2 = _block.txt 
%     {default = 1}
%   'stim_times_model': name of stimulus model
%     model names in 3dDeconvolve for -stim_times
%     allowed: 'SPMG','GAM','TENT', 'BLOCK'
%     {default = 'SPMG'}
%   'stim_times_nbasis': number of basis functions for TENT
%     {default = 10}
%   'concat_flag': [0|1|2] analyze concatenated across scans
%     0: analyze each scan individually
%     1: analyze concatenated scans only
%     2: analyze individually and concatenated
%     {default = 1}
%   'skipTRs': number of TRs at beginning of each run to ignore
%     {default = 0}
%   'minfrac': minimum fraction of a TR to register an event
%     {default = 0.1}
%   'minlag': number of TRs for minimum "lag" between stimulus and response
%     {default = 0}
%   'maxlag': number of TRs for maximum "lag" between stimulus and response
%     {default = 14}
%   'smoothsteps': smoothing steps on surface after painting
%     {default = 0}
%   'sphere_flag': [0|1] whether to sample to icosohedral sphere after painting
%     {default = 0}
%   'sphsmoothsteps': smoothing steps on spherical surface
%     applied separately from smoothsteps
%     {default = 0}
%   'pthresh': probability threshold applied to f-stats for each contrast
%     if set to 0 or 1, no thresholding applied
%     {default = 0.05}
%   'contrasts_flag': [0|1] calculate glt contrasts between each condition
%     {default = 0}
%   'iresp_flag': [0|1] output impulse response functions for each condition
%     {default = 0}
%
% Optional Parameters:
%   'snums_valid': vector of scan numbers that were processed
%     if empty, will assume all were processed
%     {default = []}
%   'infix': BOLD file name infix (e.g. '', 'corr', 'corr_resBOLD')
%     {default = []}
%   'mc_flag': [0|1] whether within-scan motion correction was done
%     Will use motion.1D files as nuisance regressors for GLM
%     {default = 1}
%   'mc_inter_flag': [0|1] whether between-scan motion correction was done
%     Allows for a single reference scan to be used for registration to FS
%     {default = 1}
%   'regFS_flag': [0|1] whether to register BOLD scans to FreeSurfer nu.mgz
%     if registration already done, will not redo
%     necessary for painting to surface
%     if 0 and not already done, paint_flag is ignored
%     {default = 0}
%   'paint_flag': [0|1|2] whether to paint results to surface
%     if 0, no painting to surface
%     if 1, paint results to surface after averaging across scans
%       (requires motion correction and register.dat file for first scan)
%     if 2, paint results to surface before averaging across scans
%       (requires register.dat file for each scan if not motion corrected)
%     {default = 1}
%   'aseg_flag': [0|1] extract ROI averages for subcortical segmentation
%     {default = 0}
%   'aseg_erode_flag': [0|1] "erode" aseg ROIs by smoothing and thresholding
%     after resampling to BOLD resolution
%     {default = 0}
%   'aseg_erode_nvoxels': smoothing kernel sigma for aseg erosion (# of voxels)
%     {default = 1}
%   'fname_aseg': name of aseg file
%      if empty, will use aseg.mgz in subjdir/subj/mri
%     {default = []}
%   'aparc_flag': [0|1] extract ROI averages for cortical surface parcellation
%     {default = 0}
%   'fnames_aparc': cell array of annotation files
%     if empty, will use ?h.aparc.annot files in subjdir/subj/label
%     otherwise, must use ?h.<name>.annot naming convention
%     {default = []}
%   'fparc_flag': extract time series for fparc cortical surface ROIs
%     ignored if fnames_fparc and fname_points are empty
%     {default = 1}
%   'fnames_fparc': cell array of annotation files in fsaverage space
%     will be resampled to to individual subject space before use
%     {default = []}
%   'render_flag': [0|1] whether to generate images of surface maps
%     {default = 0}
%   'projdist': dist (mm) to project along normal when painting
%     {default = 1}
%   'projfrac': fractional dist to project along normal when painting
%     {default = 0.5}
%   'projfrac_flag': [0|1] whether to use projfrac (1) or projdist (0)
%     {default = 0}
%   'projfrac_avg': vector of [min max del] for averaging multiple samples
%     If empty, use projfrac instead if projfrac_flag=1
%     {default = []}
%   'projdist_avg': vector of [min max del] for averaging multiple samples
%     If empty, use projdist instead
%     {default = []}
%   'resamp_flag': [0|1] whether to resample results to 1x1x1 before painting
%     {default = 0}
%   'force_repaint_flag': [0|1] whether to repaint even if output files exist
%     (do not redo volume calculations)
%     {default = 0}
%   'surfname': name of surface onto which to sample volume data
%     {default = white}
%   'forceflag': [0|1] whether to run calculations even if output files exist
%     {default = 0}
%
% Created:  12/02/16 by Don Hagler
% Last Mod: 10/27/17 by Dani Cornejo
% Last Mod: 11/03/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: check eprime time difference
%%       if only one scan, is it run 1 or run 2?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not modify this section:

% based on MMIL_Analyze_taskBOLD_Exam
% Created:  09/01/08 by Don Hagler
% Last Mod: 03/07/13 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

% check input parameters
parms = check_input(ContainerPath,FSContainerPath,fname_eprime,varargin);

% get number of scans, etc.
[parms,errcode] = check_sess(parms);
if errcode, return; end;

% check freesurfer, regFS
parms = check_fsurf(parms);

% get file names for each scan
[parms,errcode] = check_scans(parms);
if errcode, return; end;

% create stimulus time course files
parms = prep_stims(parms);

% create condition contrast structure
parms = prep_contrasts(parms);

% set stimulus file names
parms = set_stim_fnames(parms);

% create cell array with sets of scan numbers
parms = set_snum_sets(parms);

% create output directory / directories
parms = prep_outdir(parms);

% resample fparc ROIs to subject
if parms.fparc_flag
  if ~isempty(parms.fnames_fparc)
    if parms.verbose==2
%      fprintf('%s: resample_fparc_annot...\n',mfilename);
    end;
    parms = resample_fparc_annot(parms);
  else
    parms.fparc_flag = 0;
  end;
end;

% run GLM analysis
parms = run_analysis(parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_input(ContainerPath,FSContainerPath,fname_eprime,options)
  errcode = 0;
  parms = mmil_args2parms(options,{...
    'cpath',ContainerPath,[],...
    'fspath',FSContainerPath,[],...
    'fname_eprime',fname_eprime,[],...
  ...
    'outstem',[],[],...
    'analysis_outfix','analysis',[],...
    'snums',[],[],...
    'stim_times_flag',1,[0:2],...
    'stim_times_model','SPMG',{'SPMG','TENT','GAM','BLOCK'},...
    'stim_times_nbasis',10,[1,100],...
    'concat_flag',1,[0:2],...
    'contrasts_flag',false,[false true],...
    'iresp_flag',false,[false true],...
    'skipTRs',0,[0,1000],...
    'minfrac',0.1,[0,1],...
    'minlag',0,[0,10],...
    'maxlag',14,[0,30],...
    'tksmooth',10,[0,1000],...
    'smoothsteps',0,[0,1000],...
    'sphere_flag',false,[false true],...
    'sphsmoothsteps',0,[0,1000],...
    'pthresh',0.05,[0,1],...
  ...
    'stim_fnames',[],[],...
    'stim_labels',[],[],...
    'snums_valid',[],[],...
    'infix',[],[],...
    'mc_flag',true,[false true],...
    'mc_inter_flag',true,[false true],...
    'regFS_flag',false,[false true],...
    'paint_flag',1,[0,2],...
    'aseg_flag',false,[false true],...
    'aseg_erode_flag',false,[false true],...
    'aseg_erode_nvoxels',1,[1:100],...
    'fname_aseg',[],[],...
    'aparc_flag',false,[false true],...
    'fnames_aparc',[],[],...
    'fparc_flag',true,[false true],...
    'fnames_fparc',[],[],...
    'render_flag',false,[false true],...
    'projdist',1,[0,10],...
    'projfrac',0.5,[0,2],...
    'projfrac_flag',false,[false true],...
    'projdist_avg',[],[],...
    'projfrac_avg',[],[],...
    'paint_surf','white',[],...
    'resamp_flag',0,[0,1,2],...
    'force_repaint_flag',false,[false true],...
    'surfname','white',[],...
    'mask_midbrain_flag',false,[false true],...
    'forceflag',false,[false true],...    
  ... % parameters for rendering output
    'curvflag',true,[false true],...
    'curvfact',0.2,[0,1],...
    'zoom',1,[0.01,100],...
    'tif_flag',true,[false true],...
    'tif_dpi',300,[10,1000],...
    'fmin',0.1,[0,100],...
    'fmid',1,[0.001,1000],...
    'fmax',2,[0.001,1000],...
    'view_surfname','inflated',{'inflated','white','pial'},...
    'viewlist',{'lat','med'},{'left','right','lat','med','pos',...
                              'ant','sup','ven','infr','infl'},...
  ...
    'verbose',2,[0:2],...
    'motion_radius',50,[],... % for calculating distance from angle
    'min_numTRs',300,[],...
    'motion_absflag',true,[false true],...
    'fnamestem','BOLD',[],...
    'out_ext','.mgh',{'.mgh','.mgz'},...
    'revflag',false,[false true],... % only applys for pep and ipp
    'hemilist',{'lh','rh'},{'lh','rh'},...
  ... % overridden later
    'conds_contrast',[],[],...
  ... % parameter names
    'paint_tags',{'regfile' 'projfrac_flag' 'projdist' 'projfrac'...
                  'projdist_avg' 'projfrac_avg' 'mask_midbrain_flag' ...
                  'subjdir' 'surfname' 'overwrite_flag'},[],...
    'sview_tags',{'frames','subjdir','surfname','surfdir','hemi',...
                  'colorscale','cmap','fmax','fmin','fmid',...
                  'curvflag','curvfact','view','zoom',...
                  'tif_flag','tif_dpi','outstem','outdir',...
                  'visible_flag','pause_dur','forceflag'},[],...
    'deconv_tags',{'fname_motion','outdir','outstem','skipTRs','out_ext',...
                   'stim_times_flag','stim_times_model','stim_times_nbasis',...
                   'minlag','maxlag','stim_labels','contrasts_flag',...
                   'conds_contrast','iresp_flag','TR','forceflag'},[],...
    'aseg_tags',{'outdir','outstem','fname_out','csv_flag','fname_aseg',...
                 'aseg_aparc_flag','fname_vals','dispvec','disp_roicodes',...
                 'disp_scalefact','disp_suffix','dispfact','erode_flag',...
                 'erode_nvoxels','scalefact','minval','M_reg',...
                 'res_outfix','fname_colorlut','verbose','forceflag',...
                 'aseg_roilist','aparc_roilist','exclude_roilist','frames'},[],...
    'surf_roi_tags',{'fname_aparc','fname_label','fname_weights','frames',...
                     'minval','scalefact','fname_colorlut','hemi',...
                     'weights_thresh','verbose','annot_name'},[],...
  });

  parms.nhemi = length(parms.hemilist);

  if isempty(ContainerPath)
    fprintf('%s: ERROR: ContainerPath is empty\n',mfilename);
    errcode = 1;
    return;
  elseif ~exist(ContainerPath,'file')
    fprintf('%s: ERROR: ContainerPath %s not found\n',mfilename,ContainerPath);
    errcode = 1;
    return;
  end;
  if isempty(FSContainerPath)
    fprintf('%s: ERROR: FSContainerPath is empty\n',mfilename);
    errcode = 1;
    return;
  elseif ~exist(FSContainerPath,'file')
    fprintf('%s: ERROR: FSContainerPath %s not found\n',mfilename,FSContainerPath);
    errcode = 1;
    return;
  end;
    
  % set task name based on eprime file name
  [fpath,fstem] = fileparts(parms.fname_eprime);
  if ~isempty(regexpi(fstem,'MID')) 
    parms.taskname = 'MID';
  elseif ~isempty(regexpi(fstem,'SST')) 
    parms.taskname = 'SST';
  elseif ~isempty(regexpi(fstem,'WM')) 
    parms.taskname = 'nBack';
  else
    k = regexp(fstem,'\w+_(?<task>\w+)$','names');
    parms.taskname = k.task;
  end;
  if strcmp(upper(parms.taskname),'WM')
    parms.taskname = 'nBack';
  end;

  % for fs_paint
  parms.overwrite_flag = (parms.forceflag || parms.force_repaint_flag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_sess(parms)
  errcode = 0;
  % get number of scans, etc.
  [parms.ScanInfo,parms.SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(parms.cpath,...
    'snums',parms.snums_valid,'fnamestem',parms.fnamestem);
  if errcode || isempty(parms.ScanInfo)
    errcode = 1;
    return;
  end;
  fprintf('%s: %d BOLD scans in %s\n',...
    mfilename,parms.SessInfo.nscans,parms.cpath);
  if isempty(parms.snums)
    parms.snums = [1:parms.SessInfo.nscans];
  end;
  if isempty(parms.snums)
    fprintf('%s: ERROR: no valid scans\n',mfilename);
    errcode = 1;
    return;
  end;  
  % get numTRs and TR
  parms.numTRs = [parms.ScanInfo(parms.snums).nreps];   
  parms.TR = [parms.ScanInfo(parms.snums).TR]/1000;
  % exclude scans with numTRs < min_numTRs, but use the rest
  ind_keep = find(parms.numTRs>parms.min_numTRs); 
  if isempty(ind_keep)
    fprintf('%s: ERROR: numTRs below minimum (%d): %d: scans %s\n',...
            mfilename,parms.min_numTRs,parms.numTRs,strtrim(sprintf('%d ',parms.snums)));
    errcode = 1;
    return;
  end      
  parms.snums = parms.snums(ind_keep);
  parms.numTRs = parms.numTRs(ind_keep);
  parms.TR = parms.TR(ind_keep);
  % use maximum numTRs and TR if multiple scans
  if length(parms.snums)>1
    % check for mismatch in numTRs
    if length(unique(parms.numTRs))>1
      fprintf('%s: WARNING: mismatch in numTRs: scans %s\n',...
        mfilename,strtrim(sprintf('%d ',parms.snums)));
    end;
    % check for mismatch in TRs
    if length(unique(parms.TR))>1
      fprintf('%s: WARNING: mismatch in TRs: scans %s\n',...
        mfilename,strtrim(sprintf('%d ',parms.snums)));
    end;
    parms.numTRs = max(parms.numTRs);
    parms.TR = max(parms.TR);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_fsurf(parms)
  % for ROI analysis
  if parms.aseg_flag, parms.resamp_flag = 1; end;
  if parms.aparc_flag || parms.fparc_flag, parms.paint_flag = 1; end;
  % check FreeSurfer path and regT1
  if parms.paint_flag || parms.resamp_flag
    if isempty(parms.fspath)
      fprintf('%s: WARNING: setting paint_flag = 0 because fspath is empty',mfilename);
    elseif ~exist(parms.fspath,'dir')
      fprintf('%s: WARNING: setting paint_flag = 0 because fspath %s not found',...
        mfilename,parms.fspath);
    end;
    if isempty(parms.fspath) || ~exist(parms.fspath,'dir')
      parms.resamp_flag = 0; parms.aseg_flag = 0;
      parms.paint_flag = 0; parms.aparc_flag = 0; parms.fparc_flag = 0;
      return;
    end;
    % registration to T1
    if parms.regFS_flag
      fname_T1 = sprintf('%s/mri/nu.mgz',parms.fspath);
      BOLD_MMIL_Register_to_T1(...
        parms.cpath,...
        'fname_T1',fname_T1,...
        'snums',parms.snums_valid,...
        'infix',parms.infix,...
        'forceflag',parms.forceflag);
    end;
    % get subjdir from fspath
    [parms.subjdir,parms.subj,ext] = fileparts(parms.fspath);
    parms.subj = [parms.subj ext]; % in case of .'s and such
    setenv('SUBJECTS_DIR',parms.subjdir); % NOTE: not sure this is necessary
    fprintf('%s: using scan %d as reference\n',...
      mfilename,parms.SessInfo.regT1_ref);
    parms.fstem_ref = [parms.cpath '/' parms.SessInfo.fstem_regT1_ref];
    if ~isempty(parms.infix)
      parms.fstem_ref = [parms.fstem_ref '_' parms.infix];
    end;
  end;
  % check fname_aseg
  if parms.aseg_flag
    if isempty(parms.fname_aseg)
      parms.fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
    end;
    if ~exist(parms.fname_aseg,'file')
      fprintf('%s: ERROR: aseg file %s not found\n',mfilename,parms.fname_aseg);
      errcode = 1;
      return;
    end;
  else
    parms.fname_aseg = [];
  end;
  % check fnames_aparc
  parms.aparc_hemis = [];
  parms.fnames_aparc = [];
  parms.num_aparcs = 0;
  if parms.aparc_flag
    if isempty(parms.fnames_aparc)
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        parms.aparc_hemis{h} = hemi;
        parms.aparc_names{h} = 'aparc';
        parms.fnames_aparc{h} = sprintf('%s/label/%s.aparc.annot',...
          parms.fspath,hemi);
      end;
    else
      if ~iscell(parms.fnames_aparc), parms.fnames_aparc = {parms.fnames_aparc}; end;
      for f=1:length(parms.fnames_aparc)
        fname = parms.fnames_aparc{f};
        n = regexp(fname,'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
        if isempty(n)
          error('unexpected naming convention for aparc file %s\n',fname);
        end;
        parms.aparc_hemis{f} = n.hemi;
        parms.aparc_names{f} = n.name;
      end;
    end;
    parms.num_aparcs = length(parms.fnames_aparc);
    for f=1:parms.num_aparcs
      if ~exist(parms.fnames_aparc{f},'file')
        fprintf('%s: ERROR: aparc annot file %s not found\n',...
          mfilename,parms.fnames_aparc{f});
        errcode = 1;
        return;
      end;
    end;
  end;
  % check fnames_fparc
  if parms.fparc_flag && ~isempty(parms.fnames_fparc)
    if ~iscell(parms.fnames_fparc)
      parms.fnames_fparc = {parms.fnames_fparc};
    end;
    for f=1:length(parms.fnames_fparc)
      if ~exist(parms.fnames_fparc{f},'file')
        fprintf('%s: ERROR: fparc annot file %s not found\n',...
          mfilename,parms.fnames_fparc{f});
        errcode = 1;
        return;
      end;
    end;    
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_scans(parms)
  errcode = 0;
  % check input, save input file names
  parms.data_fnames = [];
  parms.data_fstems = [];
  parms.reg_fnames = [];
  parms.motion_fnames = [];
  for i=1:length(parms.snums)
    s = parms.snums(i);
    % check for bad scan num
    if ~ismember(s,parms.SessInfo.snums_valid)
      fprintf('%s: ERROR: bad BOLD Scan Num (%d)\n',mfilename,s);
      errcode = 1;
      return;
    end;
    % check motion and data files exist
    fstem = sprintf('%s/%s',parms.cpath,parms.ScanInfo(s).fstem);
    if ~isempty(parms.infix)
      fstem = [fstem '_' parms.infix];
    end;
    if parms.mc_flag
      parms.motion_fnames{i} = [fstem '_motion.1D'];
      if ~exist(parms.motion_fnames{i},'file')
        fprintf('%s: ERROR: motion 1D file %s not found\n',...
          mfilename,parms.motion_fnames{i});
        errcode = 1;
        return;
      end;
    end;
    parms.data_fstems{i} = fstem;
    parms.data_fnames{i} = sprintf('%s.mgz',fstem);
    if ~exist(parms.data_fnames{i},'file')
      fprintf('%s: ERROR: data file %s not found\n',mfilename,parms.data_fnames{i});
      errcode =1 ;
      return;
    end;
    if parms.paint_flag || parms.resamp_flag
      % check register.dat exists
      if parms.mc_inter_flag
        parms.reg_fnames{i} = sprintf('%s_register.dat',parms.fstem_ref);
      else
        parms.reg_fnames{i} = sprintf('%s_register.dat',fstem);
      end;
      if ~exist(parms.reg_fnames{i},'file')
        fprintf('%s: WARNING: reg file %s not found (cannot resample to T1 or paint to surface)\n',...
          mfilename,parms.reg_fnames{i});
        parms.resamp_flag = 0; parms.aseg_flag = 0;
        parms.paint_flag = 0; parms.aparc_flag = 0;
        parms.fparc_flag = 0;
      end;
    end;
  end;
  parms.nscans = length(parms.snums);
  if ~parms.nscans, return; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_stims(parms)
  % create 1D files from eprime file
  switch upper(parms.taskname)
    case 'MID'
      parms = prep_stims_mid(parms);
    case 'SST'
      parms = prep_stims_sst(parms);
    case 'NBACK'
      parms = prep_stims_nback(parms);
    otherwise
      error('prep_stims not implemented for %s',parms.taskname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_stims_mid(parms)
  tparms = [];
  tparms.outstem = [];
  tparms.outdir = [parms.cpath '/stim1D_MID'];
  tparms.numTRs = parms.numTRs;
  tparms.TR = parms.TR;
  tparms.minfrac = parms.minfrac;
  tparms.nskipTRs = 0; % stimuli start after skipped TRs
  tparms.forceflag = parms.forceflag;
  args = mmil_parms2args(tparms);
  abcd_extract_eprime_mid(parms.fname_eprime,args{:});
  parms.stimdir = tparms.outdir;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_stims_sst(parms)
  tparms = [];
  tparms.outstem = [];
  tparms.outdir = [parms.cpath '/stim1D_SST'];
  tparms.numTRs = parms.numTRs;
  tparms.TR = parms.TR;
  tparms.minfrac = parms.minfrac;
  tparms.nskipTRs = 0; % stimuli start after skipped TRs
  tparms.forceflag = parms.forceflag;
  args = mmil_parms2args(tparms);
  abcd_extract_eprime_sst(parms.fname_eprime,args{:});
  parms.stimdir = tparms.outdir;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_stims_nback(parms)
  tparms = [];
  tparms.outstem = [];
  tparms.outdir = [parms.cpath '/stim1D_nBack'];
  tparms.numTRs = parms.numTRs;
  tparms.TR = parms.TR;
  tparms.minfrac = parms.minfrac;
  tparms.nskipTRs = 0; % stimuli start after skipped TRs
  tparms.forceflag = parms.forceflag;
  args = mmil_parms2args(tparms);
  abcd_extract_eprime_nback(parms.fname_eprime,args{:});
  parms.stimdir = tparms.outdir;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_stim_fnames(parms)
  if any(~ismember(parms.snums,1:parms.SessInfo.nscans))
    error('bad scan numbers in snums');
  end;
  nscans = length(parms.snums);
  fprintf('%s: number of scans = %d\n',mfilename,nscans);
  if nscans==0, return; end;
  nstims = length(parms.stim_labels);
  fprintf('%s: number of stimulus conditions = %d\n',mfilename,nstims);
  parms.stim_fnames = cell(nscans,1);
  if parms.stim_times_flag == 1
    fext = '.txt';
  elseif parms.stim_times_flag == 2
    fext = '_block.txt';
  else 
    fext = '.1D';
  end;
  [tmp,fstem_eprime] = fileparts(parms.fname_eprime);
  for s=1:nscans
    stim_fnames = cell(nstims,1);
    for c=1:nstims
      stim_label = parms.stim_labels{c};
      stim_fname = sprintf('%s/%s_scan%d_%s%s',...
        parms.stimdir,fstem_eprime,s,stim_label,fext);
      if ~exist(stim_fname,'file')
        fprintf('%s: WARNING: stim file %s not found, creating null stim file...\n',... 
          mfilename,stim_fname);
        % create file with all zeros (1D) or * (txt)
        switch fext
          case '.txt'
            % write to txt file
            fid = fopen(stim_fname,'wt');
            if fid<0, error('failed to open %s for writing',stim_fname); end;
            fprintf(fid,'*\n');
            fclose(fid);
          case '.1D'
            % write to 1D file
            fid = fopen(stim_fname,'wt');
            if fid<0, error('failed to open %s for writing',stim_fname); end;
            for k=1:parms.numTRs
              fprintf(fid,'%d\n',0);
            end;
            fclose(fid);
        end;
      end;
      stim_fnames{c} = stim_fname;
    end;
    parms.stim_fnames{s} = stim_fnames;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_contrasts(parms)
  switch upper(parms.taskname)
    case 'MID'
      [parms.stim_labels,parms.conds_contrast] = abcd_set_contrasts_mid();
    case 'SST'
      [parms.stim_labels,parms.conds_contrast] = abcd_set_contrasts_sst();
    case 'NBACK'
      [parms.stim_labels,parms.conds_contrast] = abcd_set_contrasts_nback();
  end;
  parms.nstims = length(parms.stim_labels);
  parms.ncontrasts = length(parms.conds_contrast);
  parms.condnames = cat(2,parms.stim_labels,{parms.conds_contrast.name});
  parms.nconds = parms.nstims + parms.ncontrasts;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_snum_sets(parms)
  parms.snum_sets = {};
  if parms.nscans==1 || ismember(parms.concat_flag,[0,2])
    parms.snum_sets = num2cell(parms.snums);
  end;
  if parms.nscans>1 && parms.concat_flag>0
    parms.snum_sets{end+1} = parms.snums;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_outdir(parms)
  if isempty(parms.outstem)
    parms.outstem = sprintf('%s_%s',parms.fnamestem,parms.taskname);
  end;
  % create output dirs
  parms.outdirs = cell(size(parms.snum_sets));
  parms.outstems = cell(size(parms.snum_sets));
  for i=1:length(parms.snum_sets)
    s = parms.snum_sets{i};
    if length(s)>1
      parms.outstems{i} = sprintf('%s_scans',parms.outstem);
      for j=1:length(s)
        parms.outstems{i} = sprintf('%s_%d',parms.outstems{i},s(j));
      end;
    else
      parms.outstems{i} = sprintf('%s_scan_%d',parms.outstem,s);
    end;
    if ~isempty(parms.infix)
      parms.outstems{i} = sprintf('%s_%s',parms.outstems{i},parms.infix);
    end;
    parms.outdirs{i} = sprintf('%s/%s_%s',...
      parms.cpath,parms.outstems{i},parms.analysis_outfix);
    mmil_mkdir(parms.outdirs{i});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_fparc_annot(parms)
  for f=1:length(parms.fnames_fparc)
    fname_in = parms.fnames_fparc{f};
    if parms.verbose==2
      fprintf('%s: resampling annotation file %s from fsaverage to %s...\n',...
        mfilename,fname_in,parms.subj);
    end;
    % call fs_annot2annot
    tmp_parms = [];
%    tmp_parms.outdir = [parms.outdir '/fparc_annot'];
    tmp_parms.outdir = [parms.cpath '/fparc_annot'];
    tmp_parms.source_subj = 'fsaverage';
    tmp_parms.subj = parms.subj;
    tmp_parms.subjdir = parms.subjdir;
    tmp_parms.verbose = (parms.verbose==2);
    tmp_parms.forceflag = parms.forceflag;
    args = mmil_parms2args(tmp_parms);
    fname_out = fs_annot2annot(fname_in,args{:});

    n = regexp(fname_out,'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
    if isempty(n)
      error('unexpected naming convention for aparc file %s\n',fname_out);
    end;
    parms.fnames_fparc{f} = fname_out;
    parms.fparc_hemis{f} = n.hemi;
    parms.fparc_names{f} = n.name;
  end;
  parms.num_fparcs = length(parms.fnames_fparc);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = run_analysis(parms)
  for i=1:length(parms.snum_sets)
    s = parms.snum_sets{i};
    fprintf('%s: running GLM analysis for',mfilename);
    if length(s)==1
      fprintf(' scan %d\n',s);
    else
      fprintf(' scans [%s]\n',strtrim(sprintf('%d ',s)));
    end;
    % run 3dDeconvolve
    parms = run_3dDeconv(parms,i);
    % resample output to T1 resolution
    if parms.resamp_flag
      fprintf('%s: resampling output to T1 resolution...\n',mfilename);
      parms = resamp_output(parms,i);
    end;
    % calculate averages for subcortical ROIs
    if parms.aseg_flag
      fprintf('%s: aseg ROI analysis...\n',mfilename);
      for tstat_flag=0:1
        aseg_analysis(parms,i,tstat_flag);
      end;
    end;
    % paint to surface
    if parms.paint_flag
      fprintf('%s: resampling output to cortical surface...\n',mfilename);
      parms = paint_output(parms,i);
      % calculate averages for surface ROIs
      if parms.aparc_flag
        fprintf('%s: aparc ROI analysis...\n',mfilename);
        for tstat_flag=0:1
          aparc_analysis(parms,i,tstat_flag);
        end;
      end;
      if parms.fparc_flag
        fprintf('%s: fparc ROI analysis...\n',mfilename);
        for tstat_flag=0:1
          fparc_analysis(parms,i,tstat_flag);
        end;
      end;
      if parms.smoothsteps
        fprintf('%s: smoothing output...\n',mfilename);
        parms = smooth_output(parms);
      end;
      if parms.sphere_flag
        fprintf('%s: resampling output to sphere...\n',mfilename);
        sphere_output(parms);
      end;
      % create output files for each condition
      fprintf('%s: splitting output...\n',mfilename);
      split_output(parms,i);
      % calculate averages for surface ROIs
    end;
    % threhold coefs with tstats
    if parms.pthresh~=0 && parms.pthresh~=1
      fprintf('%s: thresholding output...\n',mfilename);
      threshold_output(parms);
    end;
    % render surface map images
    if parms.render_flag && parms.paint_flag
      fprintf('%s: rendering output...\n',mfilename);
      render_output(parms,i);
    end;
    % save motion stats
    if parms.mc_flag
      fname_motion = sprintf('%s/%s_motion.mat',...
        parms.outdirs{i},parms.outstems{i});
      if ~exist(fname_motion,'file') || parms.forceflag
        fprintf('%s: saving motion stats to %s...\n',mfilename,fname_motion);
        motion_stats = calc_motion_stats(parms);
        if ~isempty(motion_stats)
          save(fname_motion,'-struct','motion_stats');
        end;
      end;
    end;
    % save scan info
    fname_info = sprintf('%s/%s_info.mat',parms.outdirs{i},parms.outstems{i});
    if ~exist(fname_info,'file') || parms.forceflag
      fprintf('%s: saving scan info to %s...\n',mfilename,fname_info);
      scan_info = get_scan_info(parms,s);
      if ~isempty(scan_info)
        save(fname_info,'scan_info');
      end;
    end;
    % save dof's
    fname_dof = sprintf('%s/%s_dof.mat',parms.outdirs{i},parms.outstems{i});
    if ~exist(fname_dof,'file') || parms.forceflag
      fprintf('%s: saving degrees of freedom to %s...\n',mfilename,fname_dof);
      stats_info = load_stats_info(parms);
      ind = select_coefs(parms);
      nconds = length(ind);
      dof1 = zeros(nconds,1);
      dof2 = zeros(nconds,1);
      for j=1:nconds
        k = ind(j);
        dof(j) = stats_info(k+1).dofs;
      end;
      condnames = parms.condnames;
      save(fname_dof,'nconds','condnames','dof');
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = run_3dDeconv(parms,i)
  % select input files
  s = parms.snum_sets{i};
  [tmp,s_ind] = intersect(parms.snums,s);
  if length(s)==1
    fname_in = parms.data_fnames{s_ind};
    if parms.mc_flag
      parms.fname_motion = parms.motion_fnames{s_ind};
    else
      parms.fname_motion = [];
    end;
    stim_fnames = parms.stim_fnames{s_ind};
  else
    fname_in = parms.data_fnames(s_ind);
    if parms.mc_flag
      parms.fname_motion = parms.motion_fnames(s_ind);
    else
      parms.fname_motion = [];
    end;
    stim_fnames = parms.stim_fnames(s_ind);
  end;
  % set output
  parms.outdir = parms.outdirs{i};
  parms.outstem = parms.outstems{i};
  % run mmil_3dDeconv
  args = mmil_parms2args(parms,parms.deconv_tags);
  [parms.fname_stats,parms.fname_info,parms.fnames_iresp] = ...
    mmil_3dDeconv(fname_in,stim_fnames,args{:});
  % check output
  if ~exist(parms.fname_stats,'file')
    error('3dDeconv output file %s not found',parms.fname_stats);
  end;
  if ~exist(parms.fname_info,'file')
    error('3dDeconv info file %s not found',parms.fname_info);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resamp_output(parms,i)
  % set parameters
  s = parms.snum_sets{i};
  [tmp,s_ind] = intersect(parms.snums,s);
  tp = [];
  tp.forceflag = parms.overwrite_flag;
  tp.fname_regdat = parms.reg_fnames{s_ind(1)};
  tp.fname_ref = sprintf('%s/mri/T1.mgz',parms.fspath);
  if ~exist(tp.fname_ref,'file')
    error('file %s not found',tp.fname_ref);
  end;

  % first frame of BOLD data
  fname_reg = parms.data_fnames{s_ind(1)};
  [fpath,fstem,fext] = fileparts(fname_reg);
  if ~isempty(parms.out_ext)
    fext = parms.out_ext;
  end;
  tp.fname_out = sprintf('%s/%s_resT1_f0%s',fpath,fstem,fext);
  tp.frames = 1;
  args = mmil_parms2args(tp);
  [M_ref2reg,subj,inplane,slicethick] = ...
    mmil_resample_by_regdat(fname_reg,args{:});

  % stats
  fname_reg = parms.fname_stats;
  [fpath,fstem,fext] = fileparts(fname_reg);
  if ~isempty(parms.out_ext)
    fext = parms.out_ext;
  end;
  tp.fname_out = sprintf('%s/%s_resT1%s',fpath,fstem,fext);
  tp.frames = [];
  args = mmil_parms2args(tp);
  [M_ref2reg,subj,inplane,slicethick] = ...
    mmil_resample_by_regdat(fname_reg,args{:});
  parms.fname_stats = tp.fname_out;

  % iresp
  if parms.iresp_flag
    for c=1:length(fnames_iresp)
      fname_reg = fnames_iresp{c};
      [fpath,fstem,fext] = fileparts(fname_in);
      tp.fname_out = sprintf('%s/%s_resT1%s',fpath,fstem,fext);
      tp.frames = [];
      args = mmil_parms2args(tp);
      [M_ref2reg,subj,inplane,slicethick] = ...
        mmil_resample_by_regdat(fname_reg,args{:});
      fnames_iresp{c} = tp.fname_out;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = paint_output(parms,i)
  s = parms.snum_sets{i};
  [tmp,s_ind] = intersect(parms.snums,s);
  if parms.resamp_flag
    fname_regdat = [];
  else
    fname_regdat = parms.reg_fnames{s_ind(1)};
  end;
  args = mmil_parms2args(parms,parms.paint_tags);

  % paint GLM stats
  fname_in = parms.fname_stats;
  parms.fname_stats = fs_paint(parms.subj,fname_in,...
    'regfile',fname_regdat,args{:});

  % paint impulse response functions
  if parms.iresp_flag
    [fpath,fstem,fext] = fileparts(fname_in);
    if ~isempty(parms.out_ext)
      fext = parms.out_ext;
    end;
    for c=1:length(fnames_iresp)
      fname_in = fnames_iresp{c};
      fname_out = fs_paint(parms.subj,fname_in,...
        'regfile',fname_regdat,args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = smooth_output(parms)
  for f=1:length(parms.fname_stats)
    fname_in = parms.fname_stats{f};
    [fpath,fstem,fext] = fileparts(fname_in);
    n = regexp(fstem,'(?<stem>.+)-(?<hemi>\wh)','names');
    if ~isempty(n)
      hemi = n.hemi;
      fstem = n.stem;
    else
      hemi = [];
    end;
    fname_out = sprintf('%s/%s-sm%d-%s%s',...
      fpath,fstem,parms.smoothsteps,hemi,fext);
    parms.fname_stats{f} = fs_surf2surf(fname_in,parms.subj,...
      'fname_out',fname_out,...
      'smooth_out',parms.smoothsteps,...
      'forceflag',parms.forceflag);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sphere_output(parms)
  for f=1:length(parms.fname_stats)
    fname_in = parms.fname_stats{f};
    [fpath,fstem,fext] = fileparts(fname_in);
    n = regexp(fstem,'(?<stem>.+)-(?<hemi>\wh)','names');
    if ~isempty(n)
      hemi = n.hemi;
      fstem = n.stem;
    else
      hemi = [];
    end;
    if parms.sphsmoothsteps
      fname_out = sprintf('%s/%s-sphere-sm%d-%s%s',...
        fpath,fstem,parms.sphsmoothsteps,hemi,fext);
    else
      fname_out = sprintf('%s/%s-sphere-%s%s',...
        fpath,fstem,hemi,fext);
    end;
    fname_out = fs_surf2surf(fname_in,parms.subj,...
      'fname_out',fname_out,...
      'trgsubj','ico',...
      'icolevel',7,...
      'smooth_out',parms.sphsmoothsteps,...
      'forceflag',parms.forceflag);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stats_info = load_stats_info(parms)
  stats_info = [];
  load(parms.fname_info);
  if isempty(stats_info)
    error('file %s did not contain expected stats info',parms.fname_info);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshold_output(parms)
  stats_info = load_stats_info(parms);
  if ~iscell(parms.fname_stats), parms.fname_stats = {parms.fname_stats}; end;
  for f=1:length(parms.fname_stats)
    fname_in = parms.fname_stats{f};
    [fpath,fstem,fext] = fileparts(fname_in);
    n = regexp(fstem,'(?<stem>.+)-(?<hemi>\wh)','names');
    if ~isempty(n)
      hemi = n.hemi;
      fstem = n.stem;
    else
      hemi = [];
    end;
    if ~isempty(hemi)
      fname_in = sprintf('%s/%s-%s%s',...
        fpath,fstem,hemi,fext);
      fname_out = sprintf('%s/%s-pthresh%0.1e-%s%s',...
        fpath,fstem,parms.pthresh,hemi,fext);
    else
      fname_in = sprintf('%s/%s%s',...
        fpath,fstem,fext);
      fname_out = sprintf('%s/%s-pthresh%0.1e%s',...
        fpath,fstem,parms.pthresh,fext);
    end;
    if ~exist(fname_out,'file') || parms.forceflag
      [vol,M] = fs_load_mgh(fname_in,[],[],0,1);
      nframes = size(vol,4);
      % select coef frames (all even)
      ind = select_coefs(parms,fname_in,nframes);
      vol_out = zeros([size(vol,1),size(vol,2),size(vol,3),length(ind)]);
      for j=1:length(ind)
        k = ind(j);
        if k+1>nframes
          error('expected more frames in %s',fname_in);
        end;
        tmp_vol_coef = vol(:,:,:,k);
        % can make assumptions about where tstats are relative to coefs
        tmp_vol_tstat = vol(:,:,:,k+1);
        dofs = stats_info(k+1).dofs;
        tthresh = abs(tinv(1-parms.pthresh,dofs));
        tmp_vol_coef(abs(tmp_vol_tstat)<tthresh)=0;
        vol_out(:,:,:,j) = tmp_vol_coef;
      end;
      % save thresholded coefs in multi-frame file
      fs_save_mgh(vol_out,fname_out,M);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_coef = select_coefs(parms,fname_in,nframes)
  if ~exist('fname_in','var'), fname_in = []; end;
  if ~exist('nframes','var'), nframes = []; end;
  % check that number of frames matches twice the number of contrasts plus 1
  if ~isempty(fname_in) && ~isempty(nframes)
    nexpect = 1 + 2*parms.nconds;
    if nframes ~= nexpect
      error('mismatch between number of frames in %s (%d) and expected (%d)',...
        fname_in,nframes,nexpect);
    end;
  end;
  % select coef frames (all even)
  ind_coef = 1 + [1:2:2*parms.nconds];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [contdir,fstem] = split_output(parms,i,subj_flag,post_thresh_flag)
  if ~exist('subj_flag','var') || isempty(subj_flag), subj_flag = 0; end;
  if ~exist('post_thresh_flag','var') || isempty(post_thresh_flag), post_thresh_flag = 0; end;
  
  indir = parms.outdirs{i};
  contdir = [parms.outdir '/contrasts'];
  fstem = set_fstem(parms,i,subj_flag,post_thresh_flag);
  mmil_mkdir(contdir);
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    fname_in = sprintf('%s/%s-%s%s',indir,fstem,hemi,parms.out_ext);
    if ~exist(fname_in,'file')
      error('file %s not found',fname_in);
    end;
    vol = [];
    for k=1:parms.nconds
      condname = parms.condnames{k};
      fname = sprintf('%s/%s_%s-%s.mgz',contdir,fstem,condname,hemi);
      if ~exist(fname,'file') || parms.forceflag
        if isempty(vol)
          vol = fs_load_mgh(fname_in);
          if ~post_thresh_flag || ismember(parms.pthresh,[0,1])
            ind = select_coefs(parms,fname_in,size(vol,4));
            vol = vol(:,:,:,ind);
          end;
          nframes = size(vol,4);
          if nframes~=parms.nconds
            error('mismatch between number of stimuli (%d) plus contrasts (%d) and number of frames (%d) in %s',...
              parms.nstims,parms.ncontrasts,nframes,fname_in);
          end;
        end;
        fs_save_mgh(vol(:,:,:,k),fname);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = set_fstem(parms,i,subj_flag,post_thresh_flag)
  fstem = [parms.outstems{i} '_3dDeconv'];
  if parms.resamp_flag
    fstem = [fstem '_resT1'];
  end;
  if parms.smoothsteps
    fstem = sprintf('%s-sm%d',fstem,parms.smoothsteps);
  end;
  if ~subj_flag && parms.sphere_flag
    fstem = sprintf('%s-sphere',fstem);
    if parms.sphsmoothsteps
      fstem = sprintf('%s-sm%d',fstem,parms.sphsmoothsteps);
    end;    
  end;
  if post_thresh_flag && ~ismember(parms.pthresh,[0,1])
    fstem = sprintf('%s-pthresh%0.1e',fstem,parms.pthresh);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = aseg_analysis(parms,i,tstat_flag)
  errcode = 0;
  fname_in = parms.fname_stats;  
  [M,volsz] = mmil_load_mgh_info(fname_in,parms.forceflag,parms.outdirs{i});
  nframes = volsz(4);
  parms.outdir = parms.outdirs{i};
  parms.outstem = [parms.outstems{i} '_3dDeconv'];
  if tstat_flag
    parms.outstem = [parms.outstem '_tstat'];
  else
    parms.outstem = [parms.outstem '_beta'];
  end;
  parms.frames = select_coefs(parms,fname_in,nframes);
  if tstat_flag
    parms.frames = parms.frames + 1;
  end;
  parms.verbose = (parms.verbose>0);
  args = mmil_parms2args(parms,parms.aseg_tags);
%  try
    mmil_aseg_analysis(fname_in,parms.fspath,args{:});
%  catch
%    fprintf('\n%s: WARNING: aseg analysis failed:\n%s\n\n',...
%      mfilename,lasterr);
%    errcode = 1;
%  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fparc_analysis(parms,i,tstat_flag)
  aparc_analysis(parms,i,tstat_flag,'fparc')
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aparc_analysis(parms,i,tstat_flag,annot_name)
  errcode = 0;
  if ~exist('annot_name','var') || isempty(annot_name)
    annot_name = 'aparc';
  end;
  outstem = sprintf('%s/%s_3dDeconv',parms.outdirs{i},parms.outstems{i});
  if tstat_flag
    outstem = [outstem '_tstat'];
  else
    outstem = [outstem '_beta'];
  end;
  fname_out = [outstem '_' annot_name '_roi_data.mat'];
  if ~exist(fname_out,'file') || parms.forceflag
    roi_data = [];
    num_annot = parms.(['num_' annot_name 's']);
    for k=1:num_annot
      hemi = parms.([annot_name '_hemis']){k};
      h = find(strcmp(hemi,parms.hemilist));
      fname_in = parms.fname_stats{h};
      [M,volsz] = mmil_load_mgh_info(fname_in,parms.forceflag,parms.outdirs{i});
      nframes = volsz(4);
      parms.hemi = hemi;
      parms.annot_name = annot_name;
      parms.fname_aparc = parms.(['fnames_' annot_name]){k};
      parms.frames = select_coefs(parms,fname_in,nframes);
      parms.verbose = (parms.verbose > 0);
      if tstat_flag
        parms.frames = parms.frames + 1;
      end;
      args = mmil_parms2args(parms,parms.surf_roi_tags);
%      try
        tmp_roi_data = mmil_surf_roi(fname_in,args{:});
%      catch
%        fprintf('\n%s: WARNING: %s analysis failed:\n%s\n\n',...
%          mfilename,annot_name,lasterr);
%        errcode = 1;
%        return;
%      end;
      roi_data = [roi_data,tmp_roi_data];
    end;
    save(fname_out,'roi_data');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function render_output(parms,i)
  [contdir,fstem] = split_output(parms,i,1,1);
  parms.surfdir = [parms.outdir '/surfs'];
  parms.outdir = [parms.outdir '/plots'];
  parms.surfname = parms.view_surfname;
  mmil_mkdir(parms.outdir);
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    for k=1:parms.nconds
      condname = parms.condnames{k};
      fname = sprintf('%s/%s_%s-%s.mgz',contdir,fstem,condname,hemi);
      [fpath,instem,fext] = fileparts(fname);
      for v=1:length(parms.viewlist)
        parms.view = parms.viewlist{v};
        parms.outstem = sprintf('%s-%s-%s',...
          instem,parms.surfname,parms.view);
        fname_out = [parms.outdir '/' parms.outstem '.tif'];
        if ~exist(fname_out,'file') || parms.forceflag
          args = mmil_parms2args(parms,parms.sview_tags);
          sv_surf_view(parms.subj,fname,args{:});
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_stats = calc_motion_stats(parms)
  motion_stats = [];
  motion_data = zeros(0,6);
  max_dz = 0; max_dy = 0; max_dx = 0;
  max_rz = 0; max_rx = 0; max_ry = 0;
  mean_trans = 0; max_trans = 0; mean_rot = 0; max_rot = 0;
  mean_motion = 0; mode_motion = 0; med_motion = 0;
  min_motion = Inf; max_motion = 0; mean_accel = 0;
  subthresh_nvols = 0; subthresh_perc = 0;
  thresholds = 0; nvols = 0;

  if ~iscell(parms.fname_motion)
    fnames = {parms.fname_motion};
  else
    fnames = parms.fname_motion;
  end;
  nscans = length(fnames);

  for i=1:nscans
    fname_motion = fnames{i};
    tmp_data = mmil_load_motion_1D(fname_motion,...
      'skipTRs',parms.skipTRs,'nframes',parms.numTRs);
    % calculate difference in head position or rotation
    motion_stats = mmil_motion_stats(tmp_data,...
      parms.motion_radius,parms.motion_absflag);
    max_dx = max(max_dx,motion_stats.max_dx);
    max_dy = max(max_dy,motion_stats.max_dy);
    max_dz = max(max_dz,motion_stats.max_dz);
    max_rx = max(max_rx,motion_stats.max_rx);
    max_ry = max(max_ry,motion_stats.max_ry);
    max_rz = max(max_rz,motion_stats.max_rz);
    max_trans = max(max_trans,motion_stats.max_trans);
    max_rot = max(max_rot,motion_stats.max_rot);
    max_motion = max(max_motion,motion_stats.max_motion);
    min_motion = min(min_motion,motion_stats.min_motion);
    mean_trans = mean_trans + motion_stats.mean_trans;
    mean_rot = mean_rot + motion_stats.mean_rot;
    mean_motion = mean_motion + motion_stats.mean_motion;
    mode_motion = mode_motion + motion_stats.mode_motion;
    med_motion = med_motion + motion_stats.med_motion;
    mean_accel = mean_accel + motion_stats.mean_accel;
    subthresh_nvols = subthresh_nvols + motion_stats.subthresh_nvols;
    subthresh_perc = subthresh_perc + motion_stats.subthresh_perc;
    thresholds = motion_stats.thresholds;
    nvols = nvols + motion_stats.nvols;
    motion_data = cat(1,motion_data,tmp_data);
  end;

  if parms.nscans > 1
    mean_trans = mean_trans / parms.nscans;
    mean_rot = mean_rot / parms.nscans;
    mean_motion = mean_motion / parms.nscans;
    mode_motion = mode_motion / parms.nscans;
    med_motion = med_motion / parms.nscans;
    mean_accel = mean_accel / parms.nscans;
  end;
  
  motion_stats = struct(...
    'motion_data',motion_data,...
    'max_dx',max_dx,...
    'max_dy',max_dy,...
    'max_dz',max_dz,...
    'max_rx',max_rx,...
    'max_ry',max_ry,...
    'max_rz',max_rz,...
    'mean_trans',mean_trans,...
    'max_trans',max_trans,...
    'mean_rot',mean_rot,...
    'max_rot',max_rot,...
    'mean_motion',mean_motion,...
    'mode_motion',mode_motion,...
    'med_motion',med_motion,...
    'min_motion',min_motion,...
    'max_motion',max_motion,...
    'mean_accel',mean_accel,...
    'subthresh_nvols',subthresh_nvols,...
    'subthresh_perc',subthresh_perc,...
    'thresholds',thresholds,...
    'nvols',nvols);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scan_info = get_scan_info(parms,snums)
  nscans = length(snums);
  scan_info = [];
  scan_info.nscans = nscans;
  scan_info.snums = snums;
  scan_info.ScanInfo = parms.ScanInfo;
  scan_info.SessInfo = parms.SessInfo;
  scan_info.fnames_motion = parms.fname_motion;
  scan_info.fnames_data = {};
  for i=1:nscans
    s = snums(i);
    scan_info.fnames_data{i} = parms.data_fnames{i};
    scan_info.TRs(i) = parms.ScanInfo(s).TR/1000;
    scan_info.nreps(i) = parms.ScanInfo(s).nreps;
    % NOTE: nreps and numTRs may be identical
    %       except for APE sequence, where numTRs = nreps/2          
    [M,volsz] = mmil_load_mgh_info(parms.data_fnames{i},parms.forceflag,parms.outdir);
    scan_info.numTRs(i) = volsz(4);
    % adjust TR if numTRs ~= nreps (i.e. for APE sequence)
    if scan_info.numTRs(i) ~= scan_info.nreps(i)
      scan_info.TRs(i) = scan_info.TRs(i) * scan_info.nreps(i) / scan_info.numTRs(i);
    end;
  end;
  if nscans==1
    scan_info.fnames_data = scan_info.fnames_data{1};
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

