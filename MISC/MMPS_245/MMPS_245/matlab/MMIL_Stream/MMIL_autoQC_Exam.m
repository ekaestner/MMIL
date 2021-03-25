function MMIL_autoQC_Exam(ContainerPath,varargin)
%function MMIL_autoQC_Exam(ContainerPath,[options])
%
% Purpose: Perform automated quality control procedures on raw data
%
% Usage:
%  MMIL_autoQC_Exam(ContainerPath,'key1', value1,...);
%
% Required Input:
%  ContainerPath: full path of MRIPROC, BOLDPROC, or DTIPROC directory
%
% Optional Paramters:
%  'outdir': output directory
%    provide full path or relative to ContainerPath
%    if empty, will write output to ContainerPath
%    {default = []}
%  'proc outdir': processing file output directory
%    relative to ContainerPath
%    if empty, will write processed output to ContainerPath
%    {default = []}
%  'infix': file name infix (e.g. [], 'corr')
%    {default = []}
%  'export_flag': [0|1] create summary images for manual review
%    {default = 1}
%  'cleanupflag': [0|1] remove temporary files
%    {default = 1}
%  'verbose': [0|1] display status messages and warnings
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional File Conversion Parameters:
%  'convert_flag': [0|1] convert files to other format
%    {default = 1}
%  'out_type': output file type
%    supported types: 'nii','mgh','mgz'
%    'nii' format can be used by FSL and AFNI
%    {default = 'nii'}
%  'out_orient': output slice orientation
%    if empty or omitted, keep original orientation
%      e.g. 'LPS', 'RAS', etc.
%    for use with FSL, 'LAS' may be preferred
%    {default = []}
%
% Created:  02/20/16 by Don Hagler
% Prev Mod: 05/02/17 by Don Hagler
% Last Mod: 09/26/17 by Don Hagler
%

%% todo: save qcinfo as json file
%%       need to extract from struct like in SummarizeQC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ContainerPath,varargin);

[parms,errcode] = get_info(parms);
if errcode || ~parms.ntypes, return; end;

if ~exist(parms.outdir,'dir')
  mmil_mkdir(parms.outdir);
  if parms.outdir_group_flag
    cmd = sprintf('chmod g+w %s',parms.outdir);
    [s,r] = unix(cmd);
  end;
end;
mmil_mkdir(parms.proc_outdir);

for i=1:parms.ntypes
  stype = parms.stypes{i};
  % create list of files, check for missing files
  [fnamelist,missing_files] = check_files(stype,parms);
  if isempty(fnamelist)
    fprintf('%s: WARNING: %s file list is empty\n',mfilename,stype);
  else
    % convert files to alternative format
    if parms.convert_flag
      convert_files(fnamelist,parms);
    end;
    % check images for each scan
    auto_qcinfo = check_scans(fnamelist,stype,parms);
    if isempty(auto_qcinfo)
      fprintf('%s: WARNING: %s auto_qcinfo is empty\n',mfilename,stype);
    else
      % report auto_qcinfo
      report_auto_qcinfo(auto_qcinfo,stype,parms);
    end;
  end;
  % report missing_files
  report_missing_files(missing_files,stype,parms);
  % generate images for manual review
  if parms.export_flag && ~isempty(fnamelist)
    export_images(fnamelist,stype,parms);
  end;
end;

if parms.cleanupflag
  cleanup_mgh_files(parms);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ContainerPath,options)
  parms = mmil_args2parms(options,{...
    'indir',ContainerPath,[],...
  ...
    'outdir',[],[],...
    'proc_outdir',[],[],...
    'infix',[],[],...
    'export_flag',true,[false true],...
    'cleanupflag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % file conversion
    'convert_flag',true,[false true],...
    'out_type','nii',{'nii','nii.gz'},...
    'out_orient',[],[],...
  ... % hidden
    'in_type','mgz',{'mgh','mgz'},...
    'fMRI_minTRs',10,[2,1000],...
    'dMRI_minTRs',10,[2,1000],...
    'motion_radius',50,[],... % for calculating distance from angle
    'dti_corr_suffix','ecc_mc',[],...
    'censor_min_fract',0.8,[],... % for calculating censor_min_ndirs
    'DTfit_nonlin_flag',false,[false true],...
    'DTI_conv_meas_list',{'FA','MD','b0','V0'},[],...
    'DTI_meas_list',{'FA','MD','b0'},[],...
    'DTI_err_list',{'rms_err'},{'err','rms_err','rms'},...
    'BOLD_meas_list',{'mean','std','snr'},[],...
    'bb_buffer',20,[0,100],... % mm
    'outdir_group_flag',false,[false true],...
  ... % hidden for epi_brainmask
    'log_flag',true,[false true],...
    'thresh',0.75,[],... %% todo: may need way to choose best value for
            ...          %%       different levels of intensity inhomogeneity
            ...          %%   or, need to correct for intensity inhomogeneity
    'fill1_smooth1',10,[],...
    'fill1_thresh1',0.95,[],...
    'fill1_smooth2',30,[],... % different from default
    'fill1_thresh2',0.1,[],...
    'fill1_smooth3',0,[],...
    'fill2_smooth1',20,[],...
    'fill2_thresh1',0.95,[],...
    'fill2_smooth2',30,[],... % different from default
    'fill2_thresh2',0.1,[],...
    'fill2_smooth3',20,[],...
    'fill2_thresh3',0.1,[],... % different from default
    'binary_flag',true,[false true],...
  ... % hidden for subplots
    'planestrings',{'HOR','SAG','COR'},[],...
    'multislice_nrows',[4],[],... % may be a vector
    'multislice_min',0.25,[],...
    'multislice_max',0.75,[],...
    'multislice_plim',[1,99],[],...
    'multislice_fig_size',[8 10],[],...
    'multiframe_plim',[5,99.9],[],...
    'multiframe_slice_fract',0.5,[],...
    'multiframe_frame_fracts',linspace(0,1,81),[],...
    'multiframe_fig_size',[10 10],[],...
    'threeviews_plim',[1,99],[],...
    'threeviews_slice_fract',[0.5,0.6,0.5],[],...
    'threeviews_fig_size',[10 4],[],...
    'motion_fig_size',[10 8],[],...
    'tif_dpi',300,[],...
  ...
    'mri_convert_tags',{'out_orient','options','verbose','forceflag'},[],...
    'fsl_brainmask_tags',{'outdir','verbose','forceflag',...
                          'out_orient','ext'},[],...
    'fsl_tSNR_tags',{'outdir','mcflag','hp_sigma','TR','verbose','forceflag',...
                     'out_orient','ext'},[],...
    'epi_tags',{'log_flag','thresh',...
                'fill1_smooth1','fill1_thresh1',...
                'fill1_smooth2','fill1_thresh2',...
                'fill1_smooth3','fill1_erode_flag',...
                'clip_edges_flag','clip_edges_width',...
                'fill2_smooth1','fill2_thresh1',...
                'fill2_smooth2','fill2_thresh2',...
                'fill2_smooth3','fill2_thresh3','fill2_erode_flag',...
                'binary_flag','forceflag'},[],...
    'dti_corr_tags',{'B0unwarp_flag','ecc_flag','censor_flag',...
                     'motion_B0uw_flag','mc_flag','gradunwarp_flag',...
                     'fname_out','fname_qmat','qmat','bvals','fname_censor',...
                     'censor_min_ndirs','censor_thresh','nonlin_flag',...
                     'b0_thresh','interpm','maskoutput','verbose','forceflag',...
                     'require_B0uw_flag','optimize_B0uw_flag',...
                     'inorm_B0uw_flag','fname_B0dx','fname_for','fname_rev',...
                     'fname_B0uw_ref','fname_B0uw_mask','fname_B0uw_reg',...
                     'PhaseDir','revflag','driftcorr',...
                     'motion_B0uw_iters','min_trans',...
                     'min_rot','fname_motion_B0uw',...
                     'gruw_type','gruw_unwarpflag','gruw_isoctrflag',...
                     'fname_ref','fname_mask','fname_reg',...
                     'mstep','scales','ecc_censor_niter','mc_censor_niter',...
                     'gruw_jacobian_flag','kernelWidthMax','lambda2',...
                     'kernelWidthMax_vec','lambda2_vec','multi_opt_flag',...
                     'outfix','outext','tmpext'},[],...%,'cleanupflag'},[],...
    'convert_DTmeas_tags',{'regT1flag','zpadflag','M','RegInfo',...
                           'measlist','mgz_flag','dat_flag','dat_measlist',...
                           'nii_flag','nii_measlist','nii_out_orient',...
                           'verbose_flag','forceflag','dat_revsliceflag',...
                           'dat_permvec'},[],...
  });

  if ~isempty(parms.infix)
    parms.infix = ['_' parms.infix];
  end;

  parms.in_ext = ['.' parms.in_type];
  parms.out_ext = ['.' parms.out_type];

  if isempty(parms.outdir)
    parms.outdir = parms.indir;
  elseif mmil_isrelative(parms.outdir)
    parms.outdir = [parms.indir '/' parms.outdir];
  end;
  
  if isempty(parms.proc_outdir)
    parms.proc_outdir = parms.indir;
  elseif mmil_isrelative(parms.proc_outdir)
    parms.proc_outdir = [parms.indir '/' parms.proc_outdir];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = get_info(parms)
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(parms.indir);
  if errcode || isempty(ContainerInfo), return; end;
  parms.ScanInfo = ContainerInfo.ScanInfo;
  parms.SeriesInfo = ContainerInfo.SeriesInfo;
  stypes = fieldnames(parms.ScanInfo);
  ind = find(~cellfun(@isempty,struct2cell(parms.ScanInfo)));
  parms.stypes = stypes(ind);
  parms.ntypes = length(ind);
  parms.DTI_ScanInfo = DTI_MMIL_Get_ScanInfo(parms.indir);
  parms.BOLD_ScanInfo = BOLD_MMIL_Get_ScanInfo(parms.indir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnamelist,missing_files] = check_files(stype,parms)
  fstemlist = create_fstemlist(stype,parms);
  nfiles = length(fstemlist);
  fnamelist = cell(nfiles,1);
  missing_files = [];
  for i=1:nfiles
    fstem = fstemlist{i};
    fname = sprintf('%s/%s%s%s',...
      parms.indir,fstem,parms.infix,parms.in_ext);
    if ~exist(fname,'file')
      if parms.verbose
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
      end;
      missing_files{end+1} = fname;
    else
      fnamelist{i} = fname;
    end;
  end;
  ind = find(~cellfun(@isempty,fnamelist));
  fnamelist = fnamelist(ind);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstemlist = create_fstemlist(stype,parms)
  fstemlist = [];
  tmpinfo = parms.ScanInfo.(stype);
  if isempty(tmpinfo), return; end;
  for i=1:length(tmpinfo)
    fstem = sprintf('%s%d',stype,i);
    switch stype
      case 'BOLD'
        switch tmpinfo(i).pepolar
          case 0
            fstemlist{end+1} = [fstem '_for'];
          case 1
            fstemlist{end+1} = [fstem '_rev'];
          case 2
            fstemlist{end+1} = [fstem '_rev'];
            fstemlist{end+1} = [fstem '_for'];
          case 3            
            fstemlist{end+1} = [fstem '_for'];
            fstemlist{end+1} = [fstem '_rev'];
        end;
      case 'DTI'
        switch tmpinfo(i).pepolar
          case 0
            fstemlist{end+1} = fstem;
          case 1
            fstemlist{end+1} = [fstem '_rev'];
          case 2
            fstemlist{end+1} = [fstem '_rev'];
            fstemlist{end+1} = fstem;
          case 3            
            fstemlist{end+1} = fstem;
            fstemlist{end+1} = [fstem '_rev'];
        end;
      otherwise
        fstemlist{end+1} = fstem;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_files(fnamelist,parms)
  for i=1:length(fnamelist)
    fname_in = fnamelist{i};
    fstem = get_fstem(fname_in);
    fname_out = sprintf('%s/%s%s%s',...
      parms.proc_outdir,fstem,parms.infix,parms.out_ext);
    if strcmp(fname_in,fname_out)
      if parms.verbose
        fprintf('%s: WARNING: input and output file names are identical: %s\n',...
          mfilename,fname_in);
      end;
      continue;
    end;
    args = mmil_parms2args(parms,parms.mri_convert_tags);
    fs_mri_convert(fname_in,fname_out,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_scans(fnamelist,stype,parms)
  auto_qcinfo = [];
  % choose function for checking image quality depending on scan type
  switch stype
    case {'MPR','FLASHhi','XetaT2'}
      checkfun = @(F) check_sMRI(F,stype,parms);
    case 'BOLD'
      checkfun = @(F) check_fMRI(F,stype,parms);
    case 'DTI'
      checkfun = @(F) check_dMRI(F,stype,parms);
    otherwise
      fprintf('%s: no quality checks implemented for type %s\n',...
        mfilename,stype);
      return;
  end;
  % check image quality for each file
  clear auto_qcinfo;
  for i=1:length(fnamelist)
    auto_qcinfo(i) = checkfun(fnamelist{i});
  end;
  if ~exist('auto_qcinfo','var')
    auto_qcinfo = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_images(fnamelist,stype,parms)
  % choose function for exporting images depending on scan type
  switch stype
    case {'MPR','FLASHhi','XetaT2'}
      expfun = @(F) export_sMRI(F,stype,parms);
    case 'BOLD'
      expfun = @(F) export_fMRI(F,stype,parms);
    case 'DTI'
      expfun = @(F) export_dMRI(F,stype,parms);
    otherwise
      fprintf('%s: no quality checks implemented for type %s\n',...
        mfilename,stype);
      return;
  end;
  % export images for each file
  for i=1:length(fnamelist), expfun(fnamelist{i}); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = init_auto_qcinfo(fname,stype,parms)
  % store stype
  auto_qcinfo.stype = stype;
  % store file stem
  auto_qcinfo.fstem = get_fstem(fname);
  % get scan information
  auto_qcinfo.scaninfo = get_scaninfo(fname,stype,parms);
  % get series information
  snum = auto_qcinfo.scaninfo.SeriesIndex;
  auto_qcinfo.seriesinfo = get_seriesinfo(fname,snum,parms);
  % get image information
  auto_qcinfo.imginfo = get_imginfo(fname,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_sMRI(fname,stype,parms)
  auto_qcinfo = [];
  if parms.verbose
    fprintf('%s: checking sMRI: %s...\n',mfilename,fname);
  end;

  % initialize auto_qcinfo
  auto_qcinfo = init_auto_qcinfo(fname,stype,parms);

  % create brain mask using FSL's BET
  fname_brainmask = create_fsl_brainmask(fname,parms);

  % calculate brain volume
  auto_qcinfo.brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms);

  % calculate SNR within brain mask
  auto_qcinfo.brain_SNR = calc_brain_SNR(fname,fname_brainmask,parms);

  %% todo: only do this for T1 using head mask?
  % calculate SNR outside of brainmask
  auto_qcinfo.nonbrain_SNR = calc_nonbrain_SNR(fname,fname_brainmask,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_fMRI(fname,stype,parms)
  auto_qcinfo = [];
  if parms.verbose
    fprintf('%s: checking fMRI: %s...\n',mfilename,fname);
  end;

  % initialize auto_qcinfo
  auto_qcinfo = init_auto_qcinfo(fname,stype,parms);

  % create brain mask using epi_brainmask
  fname_brainmask = create_epi_brainmask(fname,parms);

  % calculate brain volume
  auto_qcinfo.brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms);

  % calculate SNR within brain mask
  auto_qcinfo.brain_SNR = calc_brain_SNR(fname,fname_brainmask,parms);

  if auto_qcinfo.imginfo.nreps>=parms.fMRI_minTRs
    % calculate tSNR
    fname_snr = calc_fmri_tSNR(fname,auto_qcinfo.scaninfo,parms);

    % calculate average tSNR within brain mask
    auto_qcinfo.brain_tSNR = calc_brain_tSNR(fname_snr,fname_brainmask,parms);

    % calculate mean frame to frame motion (using mcflirt output)
    auto_qcinfo.motion = calc_fmri_motion(fname,parms);
  else
    auto_qcinfo.brain_tSNR = [];
    auto_qcinfo.motion = [];
  end;
  
  %% todo: calculate FWHM
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_dMRI(fname,stype,parms)
  auto_qcinfo = [];
  if parms.verbose
    fprintf('%s: checking dMRI: %s...\n',mfilename,fname);
  end;

  % initialize auto_qcinfo
  auto_qcinfo = init_auto_qcinfo(fname,stype,parms);

  % create brain mask using epi_brainmask
  fname_brainmask = create_epi_brainmask(fname,parms);

  % calculate brain volume
  auto_qcinfo.brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms);

  % calculate SNR within brain mask
  auto_qcinfo.brain_SNR = calc_brain_SNR(fname,fname_brainmask,parms);

  if auto_qcinfo.imginfo.nreps>=parms.dMRI_minTRs
    % perform eddy current and motion correction
    %   and calculate mean frame to frame motion
    [auto_qcinfo.motion,auto_qcinfo.censor,fname_mc,fname_qmat] =...
      calc_dmri_motion(fname,auto_qcinfo.scaninfo,parms);

    % calculate diffusion tensor fit
    fname_DTfit = calc_DTfit(fname_mc,fname_qmat,parms);

    if ~isempty(fname_DTfit)
      % calculate measures from diffusion tensor fit (FA, MD, b0)
      fname_DTmeas = calc_DTmeas(fname_DTfit,parms);

      % save DTmeas as mgz and convert to nii
      convert_DTmeas(fname_DTmeas,parms);

      % calculate average DTI measures within brain mask
      auto_qcinfo.brain_DTmeas = ...
        calc_brain_DTmeas(fname_DTmeas,fname_brainmask,parms);

      % calculate residual error from diffusion tensor fit
      fname_DTerr = calc_DTerr(fname,fname_DTfit,parms);

      % save DTerr as mgz and convert to nii
      convert_DTerr(fname_DTerr,parms);

      % calculate average DTI err within brain mask
      auto_qcinfo.brain_DTerr = ...
        calc_brain_DTerr(fname,fname_DTerr,fname_brainmask,parms);
    else
      auto_qcinfo.brain_DTmeas = [];
      auto_qcinfo.brain_DTerr = [];
    end;
  else
    auto_qcinfo.motion = [];
    auto_qcinfo.censor = [];
    auto_qcinfo.brain_DTmeas = [];
    auto_qcinfo.brain_DTerr = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_sMRI(fname,stype,parms)
  if parms.verbose
    fprintf('%s: exporting sMRI: %s...\n',mfilename,fname);
  end;
  % convert to mgh for faster loading each time
  fname = unzip_mgz(fname);
  plot_multislice(fname,parms);
  plot_threeviews(fname,parms);
  plot_threeviews_mask(fname,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_dMRI(fname,stype,parms)
  if parms.verbose
    fprintf('%s: exporting dMRI: %s...\n',mfilename,fname);
  end;
  % convert to mgh for faster loading each time
  fname_orig = fname;
  fname = unzip_mgz(fname);
  plot_multislice(fname,parms);
  % multi-frame plot with universal scaling
  plot_multiframe(fname,parms,stype,0);
  % multi-frame plot with grouped scaling
  plot_multiframe(fname,parms,stype,1);
  plot_threeviews(fname,parms);
  plot_threeviews_mask(fname,parms);
  % multi-slice plot of measures (b0, FA, MD)
  plot_multislice_measures(fname,parms,parms.DTI_meas_list,parms.dti_corr_suffix);
  % three-view plot of measures (b0, FA, MD)
  plot_threeviews_measures(fname,parms,parms.DTI_meas_list,parms.dti_corr_suffix);
  % three-view plot of rgb (FA * abs(V0))
  plot_threeviews_rgb(fname,parms,parms.dti_corr_suffix);
  % multislice plots of mean err
  plot_multislice_measures(fname,parms,parms.DTI_err_list,parms.dti_corr_suffix);
  % multiframe plots of err
  plot_multiframe_measures(fname,stype,parms,parms.DTI_err_list,parms.dti_corr_suffix);
  % three-view plot of err
  plot_threeviews_measures(fname,parms,parms.DTI_err_list,parms.dti_corr_suffix);
  % plot motion time series
  imginfo = get_imginfo(fname,parms);
  if imginfo.nreps>=parms.dMRI_minTRs
    [tmp,fstem] = fileparts(fname);
    scaninfo = get_scaninfo(fname,stype,parms);
    motion_data = calc_dmri_motion(fname_orig,scaninfo,parms);
    plot_motion(fstem,motion_data,parms);
    plot_motion_fd(fstem,motion_data,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_fMRI(fname,stype,parms)
  if parms.verbose
    fprintf('%s: exporting fMRI: %s...\n',mfilename,fname);
  end;
  % convert to mgh for faster loading each time
  fname_orig = fname;
  fname = unzip_mgz(fname);
  plot_multislice(fname,parms);
  plot_multiframe(fname,parms);
  plot_threeviews(fname,parms);
  plot_threeviews_mask(fname,parms);
  plot_multislice_measures(fname,parms,parms.BOLD_meas_list);
  plot_threeviews_measures(fname,parms,parms.BOLD_meas_list);
  % plot motion time series
  imginfo = get_imginfo(fname,parms);
  if imginfo.nreps>=parms.fMRI_minTRs
    [tmp,fstem] = fileparts(fname);
    motion_data = calc_fmri_motion(fname_orig,parms);
    plot_motion(fstem,motion_data,parms);
    plot_motion_fd(fstem,motion_data,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = unzip_mgz(fname)
  [fpath,fstem] = fileparts(fname);
  fname_out = sprintf('%s/%s.mgh',fpath,fstem);
  if ~exist(fname_out,'file')
    fs_copy_mgh(fname,fname_out); % note: this is faster than fs_mri_convert
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multislice(fname,parms)
  [tmp,fstem] = fileparts(fname);
  suffix = 'ms';
  tparms = [];
  tparms.outdir = parms.outdir;
  tparms.plim = parms.multislice_plim;
  tparms.fig_size = parms.multislice_fig_size;
  tparms.tif_dpi = parms.tif_dpi;
  tparms.forceflag = parms.forceflag;
  for j=1:length(parms.multislice_nrows)
    nrows = parms.multislice_nrows(j);
    nslices = nrows.^2;
    tparms.slice_fracts = ...
      linspace(parms.multislice_min,parms.multislice_max,nslices);
    for i=1:length(parms.planestrings)
      tparms.planestring = parms.planestrings{i};
      tparms.outstem = sprintf('%s_%s_%s%d',...
        fstem,tparms.planestring,suffix,nslices);
      args = mmil_parms2args(tparms);
      mmil_multislice_subplots(fname,args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multiframe(fname,parms,stype,groupflag)
  if ~exist('stype','var'), stype = []; end;
  if ~exist('groupflag','var') || isempty(groupflag), groupflag = 0; end;
  [tmp,fstem] = fileparts(fname);
  suffix = 'mf';
  tparms = [];
  tparms.outdir = parms.outdir;
  tparms.slice_fract = parms.multiframe_slice_fract;
  tparms.frame_fracts = parms.multiframe_frame_fracts;
  tparms.plim = parms.multiframe_plim;
  tparms.fig_size = parms.multiframe_fig_size;
  tparms.tif_dpi = parms.tif_dpi;
  tparms.forceflag = parms.forceflag;
  % skip this one if there are not multiple frames
  [M,volsz] = mmil_load_mgh_info(fname,0,parms.proc_outdir);
  nf = volsz(4);
  if nf==1, return; end;
  % get bvals for dMRI to scale each group separately
  if strcmp(stype,'DTI')
    tparms.frame_fracts = []; % plot all frames
    if groupflag
      s = get_scannum(fname,stype);
      qmat = parms.DTI_ScanInfo(s).qmat;
      qlength = max(sqrt(sum(qmat.^2,2)),eps);
      bval = parms.DTI_ScanInfo(s).bval;
      bvals = bval.*(qlength.^2);
      tparms.groups = round(bvals);
      suffix = [suffix '_bvals'];
      if nf~=length(bvals), return; end;
    end;
  end;
  % create multiframe plots for each slice plane
  for i=1:length(parms.planestrings)
    tparms.planestring = parms.planestrings{i};
    tparms.outstem = sprintf('%s_%s_%s',fstem,tparms.planestring,suffix);
    args = mmil_parms2args(tparms);
    mmil_multiframe_subplots(fname,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_threeviews(fname,parms)
  [tmp,fstem] = fileparts(fname);
  suffix = 'tv';
  tparms = [];
  tparms.outdir = parms.outdir;
  tparms.outstem = sprintf('%s_%s',fstem,suffix);
  tparms.plim = parms.threeviews_plim;
  tparms.fig_size = parms.threeviews_fig_size;
  tparms.tif_dpi = parms.tif_dpi;
  tparms.planestrings = parms.planestrings;
  tparms.slice_fract = parms.threeviews_slice_fract;
  tparms.forceflag = parms.forceflag;
  args = mmil_parms2args(tparms);
  mmil_subplots(fname,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_threeviews_mask(fname,parms,infix)
  if ~exist('infix','var') || isempty(infix), infix = 'brainmask'; end;
  [tmp,fstem] = fileparts(fname);
  fname_mask = sprintf('%s/%s_%s%s',...
    parms.indir,fstem,infix,parms.in_ext);
  if ~exist(fname_mask,'file'), return; end;
  fname_mask = unzip_mgz(fname_mask);
  plot_threeviews(fname_mask,parms)
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multislice_measures(fname,parms,measlist,infix)
  [tmp,fstem] = fileparts(fname);
  if ~exist('infix','var'), infix = []; end;
  if ~isempty(infix), infix = ['_' infix]; end;
  for i=1:length(measlist)
    meas = measlist{i};
    fname_meas = sprintf('%s/%s%s_%s%s',...
      parms.indir,fstem,infix,meas,parms.in_ext);
    if ~exist(fname_meas,'file'), continue; end;
    fname_meas = unzip_mgz(fname_meas);
    plot_multislice(fname_meas,parms)
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multiframe_measures(fname,stype,parms,measlist,infix)
  [tmp,fstem] = fileparts(fname);
  if ~exist('infix','var'), infix = []; end;
  if ~isempty(infix), infix = ['_' infix]; end;
  for i=1:length(measlist)
    meas = measlist{i};
    fname_meas = sprintf('%s/%s%s_%s%s',...
      parms.indir,fstem,infix,meas,parms.in_ext);
    if ~exist(fname_meas,'file'), continue; end;
    fname_meas = unzip_mgz(fname_meas);
    plot_multiframe(fname_meas,parms,stype,0)
    plot_multiframe(fname_meas,parms,stype,1)
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_threeviews_measures(fname,parms,measlist,infix)
  [tmp,fstem] = fileparts(fname);
  if ~exist('infix','var'), infix = []; end;
  if ~isempty(infix), infix = ['_' infix]; end;
  for i=1:length(measlist)
    meas = measlist{i};
    fname_meas = sprintf('%s/%s%s_%s%s',...
      parms.indir,fstem,infix,meas,parms.in_ext);
    if ~exist(fname_meas,'file'), continue; end;
    fname_meas = unzip_mgz(fname_meas);
    plot_threeviews(fname_meas,parms)
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_threeviews_rgb(fname,parms,infix)
  [tmp,fstem] = fileparts(fname);
  suffix = 'tv';
  fname_FA = sprintf('%s/%s_%s_FA%s',parms.indir,fstem,infix,parms.in_ext);
  fname_V0 = sprintf('%s/%s_%s_V0%s',parms.indir,fstem,infix,parms.in_ext);
  if ~exist(fname_FA,'file') || ~exist(fname_V0,'file'), return; end;
  tparms = [];
  tparms.outdir = parms.outdir;
  tparms.outstem = sprintf('%s_%s_rgb_%s',fstem,infix,suffix);
  tparms.plim = parms.threeviews_plim;
  tparms.fig_size = parms.threeviews_fig_size;
  tparms.tif_dpi = parms.tif_dpi;
  tparms.planestrings = parms.planestrings;
  tparms.slice_fract = parms.threeviews_slice_fract;
  tparms.forceflag = parms.forceflag;
  args = mmil_parms2args(tparms);
  mmil_dti_rgb_subplots(fname_FA,fname_V0,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_motion(fstem,motion_data,parms)
  suffix = 'motion';
  fname_out = sprintf('%s/%s_%s.tif',parms.outdir,fstem,suffix);
  if ~exist(fname_out,'file') || parms.forceflag
    figure;
    set(gcf,'Visible','off');
    subplot(2,1,1);
    plot(motion_data.motion_tseries(:,1:3));
    axis tight;
    xlabel('time points');
    ylabel('translation (mm)');
    legend({'L-R (dx)','A-P (dy)','I-S (dz)'},'location','EastOutside');
    subplot(2,1,2);
    plot(motion_data.motion_tseries(:,4:6));
    axis tight;
    xlabel('time points');
    ylabel('rotation (deg)');
    legend({'pitch (rx)','yaw (ry)','roll (rz)'},'location','EastOutside');
    if ~isempty(parms.motion_fig_size)
      fig_size = parms.motion_fig_size;
      if length(fig_size)==1
        fig_size = [fig_size fig_size];
      end;
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperSize', [fig_size(1) fig_size(2)]);
      set(gcf, 'PaperPositionMode', 'manual');
      set(gcf, 'PaperPosition', [0 0 fig_size(1) fig_size(2)]);
    end;
    print(gcf,'-dtiff',fname_out,sprintf('-r %d',parms.tif_dpi));
    close(gcf);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_motion_fd(fstem,motion_data,parms)
  suffix = 'motion_fd';
  fname_out = sprintf('%s/%s_%s.tif',parms.outdir,fstem,suffix);
  if ~exist(fname_out,'file') || parms.forceflag
    figure;
    set(gcf,'Visible','off');
    plot(motion_data.motion_fd);
    axis tight;
    xlabel('time points');
    ylabel('frame-wise displacement (mm)');
    if ~isempty(parms.motion_fig_size)
      fig_size = parms.motion_fig_size;
      if length(fig_size)==1
        fig_size = [fig_size fig_size];
      end;
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperSize', [fig_size(1) fig_size(2)]);
      set(gcf, 'PaperPositionMode', 'manual');
      set(gcf, 'PaperPosition', [0 0 fig_size(1) fig_size(2)]);
    end;
    print(gcf,'-dtiff',fname_out,sprintf('-r %d',parms.tif_dpi));
    close(gcf);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = get_fstem(fname)
  [tmp,fstem] = fileparts(fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scaninfo = get_scaninfo(fname,stype,parms)
  scaninfo = [];
  if parms.verbose
    fprintf('%s: getting scan info...\n',mfilename);
  end;
  s = get_scannum(fname,stype);
  switch stype
    case 'DTI'
      scaninfo = parms.DTI_ScanInfo(s);
    case 'BOLD'
      scaninfo = parms.BOLD_ScanInfo(s);
    otherwise
      scaninfo = parms.ScanInfo.(stype)(s);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = get_scannum(fname,stype)
  fstem = get_fstem(fname);
  n = regexp(fstem,[stype '(?<snum>\d+)'],'names');
  if isempty(n)
    error('failed to get scan number from file name %s',fname);
  end;
  s = str2num(n.snum);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seriesinfo = get_seriesinfo(fname,snum,parms);
  seriesinfo = [];
  if parms.verbose
    fprintf('%s: getting series info...\n',mfilename);
  end;
  seriesinfo = parms.SeriesInfo(snum);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imginfo = get_imginfo(fname,parms)
  imginfo = [];
  if parms.verbose
    fprintf('%s: getting image info...\n',mfilename);
  end;
  % load mgh header info
  [imginfo.M,imginfo.volsz] =...
    mmil_load_mgh_info(fname,parms.forceflag,parms.proc_outdir);
  % determine slice orientation
  imginfo.orient = fs_read_orient([],imginfo.M);
  % calculate voxel size
  imginfo.voxsize = fs_voxsize(imginfo.M);
  % calculate voxel volume
  imginfo.voxvol = prod(imginfo.voxsize);
  % calculate number of voxels
  imginfo.nvox = prod(imginfo.volsz(1:3));
  % get number of repetitions
  imginfo.nreps = imginfo.volsz(4);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_brainmask = create_fsl_brainmask(fname,parms)
  fname_brainmask = [];
  tparms = parms;
  tparms.outdir = parms.proc_outdir;
  args = mmil_parms2args(tparms,parms.fsl_brainmask_tags);
  fname_brainmask = mmil_fsl_brainmask(fname,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_brainmask = create_epi_brainmask(fname,parms)
  [fpath,fstem,fext] = fileparts(fname);
  fname_brainmask = sprintf('%s/%s_brainmask%s',parms.proc_outdir,fstem,fext);
  args = mmil_parms2args(parms,parms.epi_tags);
  epi_brainmask(fname,fname_brainmask,args{:});
  fname_brainmask_nii = sprintf('%s/%s_brainmask.nii',parms.proc_outdir,fstem);
  args = mmil_parms2args(parms,parms.mri_convert_tags);
  fs_mri_convert(fname_brainmask,fname_brainmask_nii,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms)
  brain_vol = [];
  fstem = get_fstem(fname_brainmask);
  fname_out = sprintf('%s/%s_vol.mat',parms.proc_outdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating brain volume...\n',mfilename);
    end;
    [vol,M] = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol>0);
    nvox_brain = length(ind_brain);
    brain_vol = nvox_brain * auto_qcinfo.imginfo.voxvol;
    save(fname_out,'brain_vol','nvox_brain');
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_SNR = calc_brain_SNR(fname,fname_brainmask,parms)
  brain_SNR = [];
  fstem = get_fstem(fname);
  fname_out = sprintf('%s/%s_brain_SNR.mat',parms.proc_outdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating SNR within brain...\n',mfilename);
    end;
    % load first frame only
    [vol,M] = fs_load_mgh(fname,[],1);
    vol_brain = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol_brain>0);
    vals = vol(ind_brain);
    % calculate SNR as ratio between mean and standard deviation
    brain_SNR.mean = mean(vals);
    brain_SNR.std = std(vals);
    brain_SNR.SNR = brain_SNR.mean / brain_SNR.std;
    % calculate additional stats
    brain_SNR.min = min(vals);
    brain_SNR.max = max(vals);
    brain_SNR.median = median(vals);
    % save results
    save(fname_out,'brain_SNR');
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: create head mask
%% todo: dilate head mask
%% todo: separate ROIs for left/right, front/back, top
%%       depends on orientation -- so should probably reslice to standard
%%       and choose ROIs based on phase-encode direction

function nonbrain_SNR = calc_nonbrain_SNR(fname,fname_brainmask,parms)
  nonbrain_SNR = [];
  fstem = get_fstem(fname);
  fname_out = sprintf('%s/%s_nonbrain_SNR.mat',parms.proc_outdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating SNR outside brain...\n',mfilename);
    end;
    % load image volume
    [vol,M,tmp,volsz] = fs_load_mgh(fname);
    voxsize = fs_voxsize(M);
    % load brain mask
    vol_brain = fs_load_mgh(fname_brainmask);
    % exclude zero values
    vol_image = 1.0*(vol>0);

    % get bounds of brain mask
    bb_v = ceil(voxsize * parms.bb_buffer);
    ind_brain = find(vol_brain>0);
    [x,y,z] = ind2sub(volsz,ind_brain);
    bb_x = minmax(x');
    bb_y = minmax(y');
    bb_z = minmax(z');

    % create bounding box
    vol_box = ones(volsz);
    vol_box(1:max(1,bb_x(1)-bb_v(1)),:,:) = 0;
    vol_box(min(volsz(1),bb_x(2)+bb_v(1)):volsz(1),:,:) = 0;
    vol_box(:,1:max(1,bb_y(1)-bb_v(2)),:) = 0;
    vol_box(:,min(volsz(2),bb_y(2)+bb_v(2)):volsz(2),:) = 0;
    vol_box(:,:,1:max(1,bb_z(1)-bb_v(3))) = 0;
    vol_box(:,:,min(volsz(3),bb_z(2)+bb_v(3)):volsz(3)) = 0;

    % create negative of bounding box
    vol_nonbrain = 1.0*(vol_box==0) .* vol_image;

    ind_nonbrain = find(vol_nonbrain);
    vals = vol(ind_nonbrain);
    nonbrain_SNR.mean = mean(vals);
    nonbrain_SNR.std = std(vals);
    nonbrain_SNR.SNR = nonbrain_SNR.mean / nonbrain_SNR.std;

    save(fname_out,'nonbrain_SNR');
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_snr = calc_fmri_tSNR(fname,scaninfo,parms)
  fname_snr = [];
  tparms = parms;
  tparms.TR = scaninfo.TR/1000;
  tparms.outdir = parms.proc_outdir;
  args = mmil_parms2args(tparms,parms.fsl_tSNR_tags);
  fname_snr = mmil_fsl_calc_tSNR(fname,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_data = calc_fmri_motion(fname,parms)
  motion_data = [];
  fstem = get_fstem(fname);
  fname_par = sprintf('%s/%s_mcf.par',parms.proc_outdir,fstem);
  fname_motion = sprintf('%s/%s_motion.mat',parms.proc_outdir,fstem);
  if ~exist(fname_motion,'file') || parms.forceflag
    if ~exist(fname_par,'file'), return; end;
    motion_tseries =  mmil_load_motion_par(fname_par);
    % calculate motion stats
    [motion_stats,motion_fd] = mmil_motion_stats(motion_tseries);
    [motion_stats_nody,motion_fd_nody] = mmil_motion_stats(motion_tseries,...
                                                           [],[],[],1);
    save(fname_motion,'motion_tseries','motion_fd','motion_stats',...
                      'motion_fd_nody','motion_stats_nody');
  end;
  motion_data = load(fname_motion);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_tSNR = calc_brain_tSNR(fname_snr,fname_brainmask,parms)
  brain_tSNR = [];
  fstem = get_fstem(fname_snr);
  fname_out = sprintf('%s/%s_brain_tSNR.mat',parms.proc_outdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating tSNR within brain...\n',mfilename);
    end;
    vol_snr = fs_load_mgh(fname_snr);
    vol_brain = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol_brain>0);
    vals = vol_snr(ind_brain);
    brain_tSNR.mean = mean(vals);
    brain_tSNR.median = median(vals);
    brain_tSNR.std = std(vals);
    save(fname_out,'brain_tSNR');
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [motion_data,censor_info,fname_mc,fname_qmat] = ...
                                        calc_dmri_motion(fname,scaninfo,parms)
  fname_mc = [];
  motion_data = [];
  censor_info = [];

  % set output file name
  [fpath,fstem,fext] = fileparts(fname);
  fname_mc = sprintf('%s/%s_%s%s',...
    parms.proc_outdir,fstem,parms.dti_corr_suffix,fext);
  fname_qmat = sprintf('%s/%s_%s_qmat.mat',...
    parms.proc_outdir,fstem,parms.dti_corr_suffix);
  fname_censor = sprintf('%s/%s_censor.mat',parms.proc_outdir,fstem);

  % get diffusion directions and b-values from scaninfo
  qmat = scaninfo.qmat;
  bvals = scaninfo.bval*sum(qmat.^2,2);

  % set parameters for dti_correct_data
  parms.B0unwarp_flag = 0;
  parms.ecc_flag = 1;
  parms.censor_flag = 1;
  parms.censor_min_ndirs = parms.censor_min_fract*length(bvals);
  parms.motion_B0uw_flag = 0;
  parms.mc_flag = 1;
  parms.gradunwarp_flag = 0;
  parms.interpm = 1; % linear interpolation for better speed
  parms.qmat = qmat;
  parms.bvals = bvals;
  parms.fname_out = fname_mc;
  parms.fname_qmat = fname_qmat;

  % use dti_correct_data for eddy current distortion and motion
  if parms.verbose
   fprintf('%s: performing eddy current and motion correction using dti_correct_data...\n',...
     mfilename);
   tic;
  end;
  args = mmil_parms2args(parms,parms.dti_corr_tags);
  [fname_mc,fname_B0dx,fname_qmat] = dti_correct_data(fname,args{:});
  if parms.verbose, toc; end;

  fname_motion = sprintf('%s/%s_motion.mat',parms.proc_outdir,fstem);
  if ~exist(fname_motion,'file') || parms.forceflag
    % load motion file
    if exist(fname_qmat,'file')
      tmp = load(fname_qmat);
      motion_tseries = tmp.motion_tseries;
      % calculate motion stats
      [motion_stats,motion_fd] = mmil_motion_stats(motion_tseries);
      [motion_stats_nody,motion_fd_nody] = mmil_motion_stats(motion_tseries,...
                                                            [],[],[],1);
      save(fname_motion,'motion_tseries','motion_fd','motion_stats',...
                        'motion_fd_nody','motion_stats_nody');
    else
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_qmat);
    end;
  end;
  if exist(fname_motion,'file')
    motion_data = load(fname_motion);
  end;

  % load censoring information and summarize
  if exist(fname_censor,'file')
    censor_info = load(fname_censor);
    censor_info.nbad_frame_slices = sum(sum(censor_info.censor_mat));
    censor_info.nbad_frames = sum(sum(censor_info.censor_mat,1) > 0);
    censor_info.nbad_slices = sum(sum(censor_info.censor_mat,2) > 0);
    censor_info.nbad_frames_per_slice = sum(censor_info.censor_mat,2);
    censor_info.max_nbad_frames_per_slice = ...
      max(censor_info.nbad_frames_per_slice);
    censor_info.nbad_slices_per_frame = sum(censor_info.censor_mat,1);
    censor_info.max_nbad_slices_per_frame = ...
      max(censor_info.nbad_slices_per_frame);
  else
    fprintf('%s: WARNING: file %s not found\n',mfilename,fname_censor);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_DTfit = calc_DTfit(fname,fname_qmat,parms)
  fstem = get_fstem(fname);
  fname_DTfit = sprintf('%s/%s_DTfit.mat',parms.proc_outdir,fstem);
  if ~exist(fname_DTfit,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: loading %s...\n',mfilename,fname);
    end;
    [vol,M] = fs_load_mgh(fname);
    load(fname_qmat);
    nframes = size(qmat,1);
    if parms.verbose
      fprintf('%s: running tensor calculations with %d frames...\n',...
        mfilename,nframes);
      tic
    end;
    DTfit = dti_fit_tensor(vol,qmat,...
      'bvals',bvals,'M',M,'nonlin_flag',parms.DTfit_nonlin_flag);
    if parms.verbose, toc; end;
    if isempty(DTfit)
      fprintf('%s: WARNING: tensor calculations failed for %s\n',...
        mfilename,fname);
      fname_DTfit = [];
      return;
    end;
    % save output structure in matfile
    save(fname_DTfit,'DTfit','M'); % include M for saving mgh files
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_DTmeas = calc_DTmeas(fname_DTfit,parms)
  fname_DTmeas = regexprep(fname_DTfit,'DTfit','DTmeas');
  if ~exist(fname_DTmeas,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating DT measures...\n',mfilename);
      tic
    end;
    % load DT fit
    load(fname_DTfit);
    % calculate Eigenvectors, FA, etc.
    DTmeas = dti_calc_DTmeas(DTfit,0);
    % save output structure in matfile
    save(fname_DTmeas,'DTmeas','M'); % include M for saving mgh files
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_DTmeas(fname_DTmeas,parms)
  fstem = get_fstem(fname_DTmeas);
  fstem = regexprep(fstem,'_DTmeas','');
  fstem = sprintf('%s/%s',parms.proc_outdir,fstem);
  all_exist = 1;
  % check if output files exist
  ext_list = {'.mgz','.nii'};
  for m=1:length(parms.DTI_conv_meas_list)
    meas = parms.DTI_conv_meas_list{m};
    for k=1:length(ext_list)
      ext = ext_list{k};
      fname = sprintf('%s%s',fstem,ext);
      if ~exist(fname,'file')
        all_exist = 0;
        break;
      end;
    end;
    if ~all_exist, break; end;
  end;
  if ~all_exist || parms.forceflag
    if parms.verbose
      fprintf('%s: converting DT measures...\n',mfilename);
      tic
    end;
    load(fname_DTmeas);
    parms.M = M;
    parms.measlist = parms.DTI_conv_meas_list;
    parms.nii_flag = 1;
    parms.nii_measlist = parms.measlist;
    args = mmil_parms2args(parms,parms.convert_DTmeas_tags);
    dti_convert_DTmeas(DTmeas,fstem,args{:});
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_DTmeas = calc_brain_DTmeas(fname_DTmeas,fname_brainmask,parms)
  brain_DTmeas = [];
  fname_brain_DTmeas = regexprep(fname_DTmeas,'DTmeas','brain_DTmeas');
  if ~exist(fname_brain_DTmeas,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating DTI measures within brain...\n',mfilename);
    end;
    load(fname_DTmeas);
    vol_brain = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol_brain>0);
    for m=1:length(parms.DTI_meas_list)
      meas = parms.DTI_meas_list{m};
      if strcmp(meas,'MD')
        vol = 1000*mean(DTmeas.volE,4);
        vals = vol(ind_brain);
      else
        vals = DTmeas.(['vol' meas])(ind_brain);
      end;
      brain_DTmeas.(meas) = [];
      brain_DTmeas.(meas).mean = mean(vals);
      brain_DTmeas.(meas).median = median(vals);
      brain_DTmeas.(meas).std = std(vals);
    end;
    save(fname_brain_DTmeas,'brain_DTmeas');
  else
    load(fname_brain_DTmeas);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_DTerr = calc_DTerr(fname,fname_DTfit,parms)
  fname_DTerr = regexprep(fname_DTfit,'DTfit','DTerr');
  if ~exist(fname_DTerr,'file') || parms.forceflag
    % calculate residual error
    if parms.verbose
      fprintf('%s: calculating DT error...\n',mfilename);
      tic
    end;

    % load DT fit
    load(fname_DTfit);
    % load data
    [vol,M] = fs_load_mgh(fname);

    % calculate synthesized volume
    vol_synth = dti_synth_vol(DTfit);

    % calculate residual error for each frame
    DTerr = [];
    DTerr.vol_err = vol - vol_synth;

    % calculate RMS signal and error across frames
    % one slice at a time to save memory
    volsz = DTfit.volsz(1:3);
    nslices = volsz(3);
    DTerr.vol_rms = zeros(volsz);
    DTerr.vol_rms_err = zeros(volsz);
    for z = 1:nslices
      slc = vol(:,:,z,:);
      slc_synth = vol_synth(:,:,z,:);
      DTerr.vol_rms(:,:,z) = sqrt(mean(slc.^2,4));
      DTerr.vol_rms_err(:,:,z) = sqrt(mean((slc_synth-slc).^2,4));
    end;
    DTerr.vol_rms_err_ratio = DTerr.vol_rms_err ./ max(DTerr.vol_rms(:));

    % save output structure in matfile
    save(fname_DTerr,'DTerr','M'); % include M for saving mgh files

    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_DTerr(fname_DTerr,parms)
  fstem = get_fstem(fname_DTerr);
  fstem = regexprep(fstem,'_DTerr','');
  fstem = sprintf('%s/%s',parms.proc_outdir,fstem);
  DTerr = [];
  for m=1:length(parms.DTI_err_list)
    errname = parms.DTI_err_list{m};
    fname_mgz = sprintf('%s_%s.mgz',fstem,errname);
    fname_nii = sprintf('%s_%s.nii',fstem,errname);
    if ~exist(fname_mgz,'file') || parms.forceflag
      if isempty(DTerr)
        if parms.verbose
          fprintf('%s: converting DT err...\n',mfilename);
          tic
        end;
        load(fname_DTerr);
      end;
      vol = DTerr.(['vol_' errname]);
      fs_save_mgh(vol,fname_mgz,M);    
    end;
    % convert to nii
    args = mmil_parms2args(parms,parms.mri_convert_tags);
    fs_mri_convert(fname_mgz,fname_nii,args{:});
    if parms.verbose && ~isempty(DTerr), toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_DTerr = calc_brain_DTerr(fname,fname_DTerr,fname_brainmask,parms)
  brain_DTerr = [];
  fname_brain_DTerr = regexprep(fname_DTerr,'DTerr','brain_DTerr');
  if ~exist(fname_brain_DTerr,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating DTI residual error within brain...\n',mfilename);
    end;

    % load data
    [vol,M] = fs_load_mgh(fname);
    % load DT err
    load(fname_DTerr);
    % load brainmask
    vol_brain = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol_brain>0);

    % calc RMS signal and err within brain for each frame
    nframes = size(vol,4);
    brain_rms_vals = zeros(1,nframes);
    brain_rms_err = zeros(1,nframes);
    for f=1:nframes
      vol_tmp = vol(:,:,:,f);
      brain_rms_vals(f) = sqrt(mean(vol_tmp(ind_brain).^2));
      vol_tmp = DTerr.vol_err(:,:,:,f);
      brain_rms_err(f) = sqrt(mean(vol_tmp(ind_brain).^2));
    end;
    brain_DTerr.brain_rms_err = brain_rms_err;
    brain_DTerr.brain_rms_vals = brain_rms_vals;
    brain_DTerr.brain_rms_err_ratio = brain_rms_err ./ brain_rms_vals;

    brain_DTerr.brain_rms_err_ratio_mean = mean(brain_DTerr.brain_rms_err_ratio);
    brain_DTerr.brain_rms_err_ratio_median = median(brain_DTerr.brain_rms_err_ratio);
    brain_DTerr.brain_rms_err_ratio_std = std(brain_DTerr.brain_rms_err_ratio);

    save(fname_brain_DTerr,'brain_DTerr');
  else
    load(fname_brain_DTerr);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function report_auto_qcinfo(auto_qcinfo,stype,parms)
  % save qc info to mat file
  fname_auto_qcinfo = sprintf('%s/%s_auto_qcinfo.mat',parms.outdir,stype);
  save(fname_auto_qcinfo,'auto_qcinfo');

  % write summary of auto_qcinfo to text file
  fname_auto_qcinfo = sprintf('%s/%s_auto_qcinfo.txt',parms.outdir,stype);
  fid = fopen(fname_auto_qcinfo,'wt');
  if fid<0, error('failed to open %s for writing',fname_auto_qcinfo); end;
  for i=1:length(auto_qcinfo)
    fprintf(fid,'%s\n',auto_qcinfo(i).fstem);
    fprintf(fid,'  nvox = %s\n',sprintf('%d ',auto_qcinfo(i).imginfo.volsz(1:3)));
    fprintf(fid,'  nreps = %d\n',auto_qcinfo(i).imginfo.nreps);
    fprintf(fid,'  orient = %s\n',auto_qcinfo(i).imginfo.orient);
    fprintf(fid,'  voxel volume = %0.1f\n',auto_qcinfo(i).imginfo.voxvol);
    fprintf(fid,'  brain volume = %0.1f\n',auto_qcinfo(i).brain_vol);

    % write SNR inside brain
    if ~isempty(auto_qcinfo(i).brain_SNR)
      fprintf(fid,'  brain mean intensity = %0.3f\n',...
        auto_qcinfo(i).brain_SNR.mean);
      fprintf(fid,'  brain standard deviation of intensity = %0.3f\n',...
        auto_qcinfo(i).brain_SNR.std);
      fprintf(fid,'  brain SNR = %0.3f\n',...
        auto_qcinfo(i).brain_SNR.SNR);
      fprintf(fid,'  brain min intensity = %0.3f\n',...
        auto_qcinfo(i).brain_SNR.min);
      fprintf(fid,'  brain max intensity = %0.3f\n',...
        auto_qcinfo(i).brain_SNR.max);
      fprintf(fid,'  brain median intensity = %0.3f\n',...
        auto_qcinfo(i).brain_SNR.median);
    end;
    
    % write stype-specific metrics
    switch stype
      case {'MPR','FLASHhi','XetaT2'}
        % write SNR outside brain
        if ~isempty(auto_qcinfo(i).nonbrain_SNR)
          fprintf(fid,'  non-brain mean intensity = %0.3f\n',...
            auto_qcinfo(i).nonbrain_SNR.mean);
          fprintf(fid,'  non-brain standard deviation of intensity = %0.3f\n',...
            auto_qcinfo(i).nonbrain_SNR.std);
          fprintf(fid,'  non-brain SNR = %0.3f\n',...
            auto_qcinfo(i).nonbrain_SNR.SNR);
        end;
      case 'BOLD'
        % write TR
        TR = auto_qcinfo(i).scaninfo.TR/1000;
        fprintf(fid,'  TR = %0.3f s\n',TR);
        % write motion stats, etc.
        if ~isempty(auto_qcinfo(i).motion)
          fprintf(fid,'  mean frame-wise displacement = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.mean_motion);
          fprintf(fid,'  standard deviation of frame-wise displacement = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.std_motion);
          fprintf(fid,'  mean frame-to-frame translation = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.mean_trans);
          fprintf(fid,'  mean frame-to-frame rotation = %0.3f deg\n',...
            auto_qcinfo(i).motion.motion_stats.mean_rot);
          fprintf(fid,'  max displamcent in x relative to baseline = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.max_dx);
          fprintf(fid,'  max displamcent in y relative to baseline = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.max_dy);
          fprintf(fid,'  max displamcent in z relative to baseline = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.max_dz);
          fprintf(fid,'  max rotation about x relative to baseline = %0.3f deg\n',...
            auto_qcinfo(i).motion.motion_stats.max_rx);
          fprintf(fid,'  max rotation about y relative to baseline = %0.3f deg\n',...
            auto_qcinfo(i).motion.motion_stats.max_ry);
          fprintf(fid,'  max rotation about z relative to baseline = %0.3f deg\n',...
            auto_qcinfo(i).motion.motion_stats.max_rz);
          thresholds = auto_qcinfo(i).motion.motion_stats.thresholds;
          nthresh = length(thresholds);
          for t=1:nthresh
            fprintf(fid,'  scan time with motion < %0.1f mm = %0.3f s\n',...
              thresholds(t),TR*auto_qcinfo(i).motion.motion_stats.subthresh_nvols(t));
            fprintf(fid,'  percent of scan with motion < %0.1f mm = %0.3f\n',...
              thresholds(t),auto_qcinfo(i).motion.motion_stats.subthresh_perc(t));
          end;
        end;
        % write tSNR inside brain
        if ~isempty(auto_qcinfo(i).brain_tSNR)
          fprintf(fid,'  brain tSNR mean = %0.3f\n',...
            auto_qcinfo(i).brain_tSNR.mean);
          fprintf(fid,'  brain tSNR median = %0.3f\n',...
            auto_qcinfo(i).brain_tSNR.median);
          fprintf(fid,'  brain tSNR std = %0.3f\n',...
            auto_qcinfo(i).brain_tSNR.std);
        end;
      case 'DTI'
        % write TR
        TR = auto_qcinfo(i).scaninfo.TR/1000;
        fprintf(fid,'  TR = %0.3f s\n',TR);
        % write motion stats, etc.
        if ~isempty(auto_qcinfo(i).motion)
          fprintf(fid,'  mean frame-wise displacement = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.mean_motion);
          fprintf(fid,'  standard deviation of frame-wise displacement = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.std_motion);
          fprintf(fid,'  mean frame-to-frame translation = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.mean_trans);
          fprintf(fid,'  mean frame-to-frame rotation = %0.3f deg\n',...
            auto_qcinfo(i).motion.motion_stats.mean_rot);
          fprintf(fid,'  max displamcent in x relative to baseline = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.max_dx);
          fprintf(fid,'  max displamcent in y relative to baseline = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.max_dy);
          fprintf(fid,'  max displamcent in z relative to baseline = %0.3f mm\n',...
            auto_qcinfo(i).motion.motion_stats.max_dz);
          fprintf(fid,'  max rotation about x relative to baseline = %0.3f deg\n',...
            auto_qcinfo(i).motion.motion_stats.max_rx);
          fprintf(fid,'  max rotation about y relative to baseline = %0.3f deg\n',...
            auto_qcinfo(i).motion.motion_stats.max_ry);
          fprintf(fid,'  max rotation about z relative to baseline = %0.3f deg\n',...
            auto_qcinfo(i).motion.motion_stats.max_rz);
          thresholds = auto_qcinfo(i).motion.motion_stats.thresholds;
          nthresh = length(thresholds);
          for t=1:nthresh
            fprintf(fid,'  scan time with motion < %0.1f mm = %0.3f s\n',...
              thresholds(t),TR*auto_qcinfo(i).motion.motion_stats.subthresh_nvols(t));
            fprintf(fid,'  percent of scan with motion < %0.1f mm = %0.3f\n',...
              thresholds(t),auto_qcinfo(i).motion.motion_stats.subthresh_perc(t));
          end;
        end;
        % write number of bad slices
        if ~isempty(auto_qcinfo(i).censor)
          fprintf(fid,'  number of censored frame-slices = %d\n',...
            auto_qcinfo(i).censor.nbad_frame_slices);
          fprintf(fid,'  number of frames with censored slices = %d\n',...
            auto_qcinfo(i).censor.nbad_frames);
          fprintf(fid,'  number of slices with censored frames = %d\n',...
            auto_qcinfo(i).censor.nbad_slices);
          fprintf(fid,'  max number of censored frames per slice = %d\n',...
            auto_qcinfo(i).censor.max_nbad_frames_per_slice);
          fprintf(fid,'  max number of censored slices per frame = %d\n',...
            auto_qcinfo(i).censor.max_nbad_slices_per_frame);
        end;
        % write about average brain DTmeas (e.g. FA, MD, b0)
        if ~isempty(auto_qcinfo(i).brain_DTmeas)
          for m=1:length(parms.DTI_meas_list)
            meas = parms.DTI_meas_list{m};
            fprintf(fid,'  brain DTI %s mean = %0.3f\n',...
              meas,auto_qcinfo(i).brain_DTmeas.(meas).mean);
            fprintf(fid,'  brain DTI %s median = %0.3f\n',...
              meas,auto_qcinfo(i).brain_DTmeas.(meas).median);
            fprintf(fid,'  brain DTI %s std = %0.3f\n',...
              meas,auto_qcinfo(i).brain_DTmeas.(meas).std);
          end;
        end;
        % write about average brain DTerr
        if ~isempty(auto_qcinfo(i).brain_DTerr)
          fprintf(fid,'  brain DTI err mean = %0.3f\n',...
            auto_qcinfo(i).brain_DTerr.brain_rms_err_ratio_mean);
          fprintf(fid,'  brain DTI err median = %0.3f\n',...
            auto_qcinfo(i).brain_DTerr.brain_rms_err_ratio_median);
          fprintf(fid,'  brain DTI err std = %0.3f\n',...
            auto_qcinfo(i).brain_DTerr.brain_rms_err_ratio_std);
        end;
      otherwise
        fprintf('%s: no quality checks implemented for type %s\n',...
          mfilename,stype);
        return;
    end;
  end;
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function report_missing_files(missing_files,stype,parms)
  if ~isempty(missing_files)
    fname_out = sprintf('%s/%s_missing_files.txt',parms.outdir,stype);
    fid = fopen(fname_out,'wt');
    if fid<0, error('failed to open %s for writing',fname_out); end;
    if parms.verbose
      fprintf('%s: WARNING: %d missing %s files\n',...
        mfilename,length(missing_files),stype);
    end;
    for i=1:length(missing_files)
      fprintf(fid,'%s\n',missing_files{i});
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_mgh_files(parms)
  dlist = dir(sprintf('%s/*.mgh',parms.indir));
  if ~isempty(dlist)
    cmd = sprintf('rm %s/*.mgh',parms.indir);
    [s,r] = unix(cmd);
    if s
      error('cmd %s failed:\n%s',cmd,r);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
