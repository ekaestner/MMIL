function errcode = ABCD_Analyze_rsBOLD_Exam(ContainerPath,FSContainerPath,varargin)
%function errcode = ABCD_Analyze_rsBOLD_Exam(ContainerPath,FSContainerPath,[options])
%
% Purpose: Correlation analysis for resting-state fMRI data
%   Samples BOLD data to surface, extracts average time series
%     for each aparc cortical ROI, calculates cross-correlation between ROIs,
%     creates cortical surface correlation maps
%
% Usage:
%  ABCD_Analyze_rsBOLD_Exam(ContainerPath,FSContainerPath,'key1', value1,...)
%  e.g. ABCD_Analyze_rsBOLD_Exam(ContainerPath,FSContainerPath,...
%         'snums',[3,4],'infix','corr');
%
% Required input:
%  ContainerPath: Full path of BOLDPROC directory containing BOLD scans
%  FSContainerPath: full path of FSURF directory
%
% Optional parameters:
%  'outdir': output directory
%    provide full path or will be relative to ContainerPath
%    {default = 'rsBOLD_analysis'}
%  'outstem': file stem for ROI and corr output file names
%    {default = 'rsBOLD'}
%  'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%  'forceflag': [0|1] whether to run calculations even if output files exist
%    {default = 0}
%
% Optional parameters for data selection:
%  'snums': vector of scan numbers to analyze
%    if empty, will use all valid BOLD scans in ContainerPath
%    {default = []}
%  'snums_valid': vector of scan numbers that were processed
%    if empty, will assume all were processed
%    {default = []}
%  'infix': BOLD file name infix (e.g. '', 'corr')
%    {default = []}
%  'revflag': [0|1|2] specify whether to use non-rev or rev data
%    0: process only forward phase-encode polarity data
%    1: process only reverse phase-encode polarity data
%    2: use both forward and reverse data
%    {default = 2}
%  'min_numTRs': minimum number of TRs required to include data file
%    {default = 20}
%  'max_nscans': maximum number of scans allowed (use lowest scan numbers)
%    if empty, use all available scans
%    {defalut = []}
%
% Optional parameters for preprocessing:
%  'mc_flag': [0|1] whether within-scan motion correction was done
%    0: motion correction not done or ignore motion correction
%    1: collect motion estimates but do not regress out
%    2: collect motion estimates and regress out
%    {default = 2}
%  'mc_resp_flag': [0|1] notch filter motion estimates for aliased respiration signals
%    ignored if mc_flag = 0
%    note: this is done after nuisance regression, but before censoring
%    {default = 1}
%  'mc_bandpass_flag': [0|1] bandpass filter motion estimates
%    ignored if mc_flag = 0 or bandpass_flag = 0
%    note: this is done after nuisance regression and censoring
%    {default = 1}
%  'motion_radius': approximate radius of typical head (mm)
%    for calculating distance from an angle of rotation
%    {default = 50}
%  'motion_absflag': calculate frame-to-frame motion as sum of abs values
%    otherwise, calculate square root of sum of squares
%    {default = 1}
%  'motion_nodyflag': calculate motion ignoring dy (translation in y)
%    {default = 0}
%  'mask_flag': [0|1] mask input volumes based on mean intensities
%    {default = 0}
%  'mask_thresh': relative threshold applied to mean volume
%    cumulative probability of intensity (max is 1)
%    {default = 0.7}
%  'norm_flag': [0|1] normalize input timeseries by mean for each voxel
%    NOTE: if detrend>0, new mean will be 0, but values will be % change
%    {default = 1}
%  'normval':  value to be assigned to mean
%    {defalut = 100}
%  'thresh': when normalizing to mean, set voxels with original values
%    less than this to 0
%    {default = 10}
%  'detrend': [0|1|2] whether and how to detrend input timeseries
%    0: no detrend (not recommended)
%    1: linear detrend
%    2: quadratic detrend
%    3: cubic detrend
%    {default = 2}
%  'skipTRs': number of TRs at beginning of scan to be removed
%    {default = 0}
%  'resp_low': respiration rate low frequency cut-off (respirations/minute)
%    {default = 18.6}
%  'resp_high': respiration rate high frequency cut-off (respirations/minute)
%    {default = 25.7}
%  'bandpass_flag': [0|1] bandpass filter time series
%    {default = 1}
%  'bandpass_low': low frequency cut-off (Hz)
%    {default = 0.01}
%  'bandpass_high': high frequency cut-off (Hz)
%    {default = 0.08}
%  'bandpass_tb': bandpass filter transition band (Hz)
%    {default = 0.001}
%  'bandpass_zpad_flag': [0|1] pad time series with zeros on each end to
%     supress edge artifacts
%     {default = 1}
%  'wm_flag': [0|1|2] extract and regress out time course of white matter
%    0: do nothing
%    1: extract white matter time course but do not regress out
%    2: extract white matter time course and regress out
%    {default = 2}
%  'wm_codes': aseg ROI codes for white matter
%    {default = [2,41]}
%  'wm_erode_flag': [0|1] "erode" white matter mask by smoothing and thresholding
%    after resampling to BOLD resolution
%    {default = 1}
%  'wm_erode_nvoxels': smoothing kernel sigma for wm erosion (# of voxels)
%    {default = 1}
%  'csf_flag': [0|1|2] extract and regress out time course of CSF
%    0: do nothing
%    1: extract CSF time course but do not regress out
%    2: extract CSF time course and regress out
%    {default = 0}
%  'csf_codes': aseg ROI codes for CSF (i.e. ventricles)
%    {default = [4,43]}
%  'csf_erode_flag': [0|1] "erode" CSF mask by smoothing and thresholding
%    after resampling to BOLD resolution
%    {default = 1}
%  'csf_erode_nvoxels': smoothing kernel sigma for csf erosion (# of voxels)
%    {default = 1}
%  'brain_flag': [0|1|2] extract and regress out time course of whole brain
%    0: do nothing
%    1: extract brain time course but do not regress out
%    2: extract brain time course and regress out
%    {default = 0}
%  'brain_codes': aseg ROI codes for whole brain
%    {default = [1:255]}
%  'brain_erode_flag': [0|1] "erode" brain mask by smoothing and thresholding
%    after resampling to BOLD resolution
%    {default = 0}
%  'brain_erode_nvoxels': smoothing kernel sigma for brain erosion (# of voxels)
%    {default = 1}
%  'deriv_flag': [0|1] include derivatives of motion, white matter, CSF,
%    and brain as regressors (if mc_flag, wm_flag, csf_flag, or brain_flag)
%    {default = 1}
%  'concat_flag': [0|1] concatenate data across scans before calculating correlations
%    otherwise, calculate correlations for each scan and then average
%    {default = 0}
%
% Optional parameters for sampling to surface:
%  'projfrac_flag': [0|1] whether to use projdist (0) or projfract (1)
%    {default = 0}
%  'projdist': distance (mm) to project along surface vertex normal
%    {default = 1}
%  'projfrac': fractional distance to project along surface vertex normal
%    relative to cortical thickness
%    {default = 0.5}
%  'projfrac_avg': vector of [min max del] for averaging multiple samples
%     with mri_surf2surf projfract-avg option
%    If empty, use projfrac instead if projfrac_flag=1
%    {default = []}
%  'projdist_avg': vector of [min max del] for averaging multiple samples
%     with mri_surf2surf projdist-avg option
%    If empty, use projdist instead
%    {default = []}
%  'surfname': name of surface onto which to sample volume data
%    {default = white}
%  'ico_flag': [0|1] resample data to icosahedral sphere
%    {default = 1}
%  'ico_order': icosahedral order (0-7)
%      Order  Number of Vertices
%        0              12
%        1              42
%        2             162
%        3             642
%        4            2562
%        5           10242
%        6           40962
%        7          163842
%    {default = 3}
%  'ico_trunc_flag': [0|1] for ico_order<7, truncate extra values
%     otherwise, actually resample to ico_order
%    {default = 0}
%  'ico_presmooth': number of smoothing steps applied on native surface
%    before resampling to ico
%    NOTE: FWHM ~ 1.25 * sqrt(N)
%    {default = 64}
%
% Optional parameters for extracting ROI time courses:
%  'roi_outfix': additional string added to ROI output files
%    e.g. to use existing preprocessed data and surface time courses
%         but rerun ROI extraction with different settings
%    {default = []}
%  'aparc_flag': extract time series for aparc cortical surface ROIs
%    {default = 1}
%  'fnames_aparc': cell array of annotation files
%     if empty, will use ?h.aparc.annot files in subjdir/subj/label
%     otherwise, must use ?h.<name>.annot naming convention
%    {default = []}
%  'fparc_flag': extract time series for fparc cortical surface ROIs
%    ignored if fnames_fparc and fname_points are empty
%    {default = 1}
%  'fnames_fparc': cell array of annotation files in fsaverage space
%    will be resampled to to individual subject space before use
%    {default = []}
%  'fname_points': name of csv file containing talairach points
%     will be sampled onto fsaverage surface, resampled to individual subject
%     and then combined to generate annot files (see fs_tal2annot)
%    {default = []}
%  'aseg_flag': extract time series for aseg volume ROIs
%    {default = 1}
%  'aseg_erode_flag': [0|1] "erode" aseg ROIs by smoothing and thresholding
%    after resampling to BOLD resolution
%    {default = 0}
%  'aseg_erode_nvoxels': smoothing kernel sigma for aseg erosion (# of voxels)
%    {default = 1}
%  'fname_aseg': name of aseg file
%     if empty, will use aseg.mgz in subjdir/subj/mri
%    {default = []}
%  'csv_tseries_flag': [0|1] write ROI times series to csv spreadsheet
%    {default = 1}
%
% Optional parameters for correlation calculations:
%  'corr_outfix': additional string added to correlation output files
%    e.g. to use existing data and ROI time courses
%         but rerun correlations with different settings
%    {default = []}
%  'corr_seed_roinames': cell array of ROI names to be used as seeds
%    if empty, will use each ROI as a seed for corrleation calculations
%    if names are partial match to more than one ROI,
%      time courses will be averaged
%    {default = []}
%  'corr_roi_flag': [0|1] calculate cross-correlation for ROIs
%    {default = 1}
%  'corr_roi_ico_flag': [0|1] calculate correlation between
%    ROIs and ico vertices
%    {default = 1}
%  'corr_roi_ico_maps_flag': [0|1] save ico correlation maps
%    separately for each ROI seed
%    ignored if corr_roi_ico_flag = 0
%    {default = 0}
%  'corr_ico_flag': [0|1] calculate cross-correlation for ico vertices
%    {default = 1}
%  'censor_thresh': framewise displacement censoring threshold (mm)
%     for correlation calculations
%     motion censoring not done if mc_flag = 0 or censor_thresh = Inf
%     {default = 0.2}
%  'contig_thresh': mininum number of contiguous frames required
%     for motion censoring for correlation calculations
%     {default = 5}
%  'resid_censor_thresh': framewise displacement censoring threshold (mm)
%     for pre-residualization
%     {default = 0.3}
%  'resid_contig_thresh': mininum number of contiguous frames required
%     for motion censoring for pre-residualization
%     {default = 5}
%
% Created:  07/03/17 by Don Hagler
% Prev Mod: 10/23/17 by Don Hagler
% Last Mod: 11/02/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% based on BOLD_MMIL_Resting_Analysis
% Created:  03/06/12 by Don Hagler
% Last Mod: 06/30/17 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: remove spikes? -- AFNI has function 3dDespike

%% todo: plots of motion, wm_tseries

%% todo: accept weighted ROIs (multi-frame mgz) instead of aparc

%% todo: bar plots of correlation, option bar_plot_roinames

%% todo: mc_inter_flag -- if 0, different regfiles for each scan

%% todo: residualize motion but do partial correlation with GLM for wm?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters and check for problems
[parms,errcode] = check_input(ContainerPath,FSContainerPath,varargin);
if errcode, return; end;

% find processed BOLD data
if parms.verbose==2
%  fprintf('%s: locate_files...\n',mfilename);
end;
[parms,errcode] = locate_files(parms);
if errcode, return; end;

% resample aparc ROIs to subject
if parms.fparc_flag
  if ~isempty(parms.fnames_fparc)
    if parms.verbose==2
%      fprintf('%s: resample_fparc_annot...\n',mfilename);
    end;
    parms = resample_fparc_annot(parms);
  end;
  if ~isempty(parms.fname_points)
    if parms.verbose==2
%      fprintf('%s: create_points_annot...\n',mfilename);
    end;
    parms = create_points_annot(parms);
  end;
end;

% extract white matter, csf, and brain time series
if parms.verbose==2
%  fprintf('%s: extract_group_tseries...\n',mfilename);
end;
for g=1:length(parms.group_names)
  group_name = parms.group_names{g};
  if parms.([group_name '_flag'])>0
    parms = extract_group_tseries(parms,group_name);
  end;
end;

% calculate mean for each voxel
if parms.verbose==2
%  fprintf('%s: calc_mean...\n',mfilename);
end;
parms = calc_mean(parms);

% preprocess BOLD data
if parms.verbose==2
%  fprintf('%s: prep_data...\n',mfilename);
end;
parms = prep_data(parms);

% save info about scans
if parms.verbose==2
%  fprintf('%s: save_info...\n',mfilename);
end;
parms = save_info(parms);

% calculate motion statistics
if parms.mc_flag
  if parms.verbose==2
%    fprintf('%s: calc_motion_stats...\n',mfilename);
  end;
  parms = calc_motion_stats(parms);
end;

% concatenate data from multiple BOLD scans
if parms.concat_flag
  if parms.verbose==2
%    fprintf('%s: concat_data...\n',mfilename);
  end;
  parms = concat_data(parms);
  nruns = 1;
else
  nruns = parms.nscans;
end;

% concatenate motion data from multiple BOLD scans
if parms.verbose==2
%  fprintf('%s: concat_motion_data...\n',mfilename);
end;
parms = concat_motion_data(parms);

% loop over runs
for i=1:nruns
  if ~parms.concat_flag
    parms.fname_BOLD = parms.fnames_BOLD{i};
    parms.fname_motion = parms.fnames_motion_stats{i};
    parms.run_outfix = sprintf('_scan%d',i);
  end;
  
  % sample BOLD data to cortical surface
  if parms.aparc_flag || parms.fparc_flag || parms.ico_flag
    if parms.verbose==2
%      fprintf('%s: paint_data...\n',mfilename);
    end;
    parms = paint_data(parms);
  end;

  % extract time series for surface and/or volume ROIs
  if parms.aparc_flag || parms.fparc_flag || parms.aseg_flag
    if parms.verbose==2
%      fprintf('%s: extract_tseries...\n',mfilename);
    end;
    parms = extract_tseries(parms);
  end;

  % resample to ico
  if parms.ico_flag
    if parms.verbose==2
%      fprintf('%s: resample_ico...\n',mfilename);
    end;
    parms = resample_ico(parms);
  end;

  % calculate signal variance
  if parms.verbose==2
%    fprintf('%s: calc_variance...\n',mfilename);
  end;
  parms = calc_variance(parms);
  
  % calculate correlations
  if parms.verbose==2
%    fprintf('%s: calc_correlations...\n',mfilename);
  end;
  parms = calc_correlations(parms);
end;

if ~parms.concat_flag
  % calculate average variance
  if parms.verbose==2
%    fprintf('%s: calc_avg_variance...\n',mfilename);
  end;
  parms = calc_avg_variance(parms);

  % calculate average correlation
  if parms.verbose==2
%    fprintf('%s: calc_avg_correlations...\n',mfilename);
  end;
  parms = calc_avg_correlations(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_input(ContainerPath,FSContainerPath,args)
  errcode = 0;
  parms_filter = {...
    'ContainerPath',ContainerPath,[],...
    'FSContainerPath',FSContainerPath,[],...
...
    'outdir','rsBOLD_analysis',[],...
    'outstem','rsBOLD',[],...
    'verbose',1,[0:2],...
    'forceflag',false,[false true],...
... % for data selection
    'snums',[],[],...
    'snums_valid',[],[],...
    'infix',[],[],...
    'revflag',2,[0,1,2],...
    'min_numTRs',20,[1,1000],...
    'max_nscans',[],[1,Inf],...
... % for preprocessing
    'mc_flag',2,[0:2],...
    'mc_resp_flag',true,[false true],...
    'mc_bandpass_flag',true,[false true],...
    'mask_flag',false,[false true],...
    'mask_thresh',0.7,[],...
    'norm_flag',true,[false true],...
    'normval',100,[1 1e10],...
    'thresh',10,[0 Inf],...
    'detrend',2,[0:2],...
    'skipTRs',0,[0,100],...
    'resp_low',18.6,[],...
    'resp_high',25.7,[],...
    'bandpass_flag',true,[false true],...
    'bandpass_low',0.01,[],...
    'bandpass_high',0.08,[],...
    'bandpass_tb',0.001,[],...
    'bandpass_zpad_flag',true,[false true],...
    'wm_flag',2,[0:2],...
    'wm_codes',[2,41],[1,Inf],...
    'wm_erode_flag',true,[false true],...
    'wm_erode_nvoxels',1,[1:100],...
    'csf_flag',0,[0:2],...
    'csf_codes',[4,43],[1,Inf],...
    'csf_erode_flag',true,[false true],...
    'csf_erode_nvoxels',1,[1:100],...
    'brain_flag',0,[0:2],...
    'brain_codes',[1:255],[1,Inf],...
    'brain_erode_flag',false,[false true],...
    'brain_erode_nvoxels',1,[1:100],...
    'deriv_flag',true,[false true],...
    'concat_flag',false,[false true],...
... % for sampling to surface
    'projfrac_flag',false,[false true],...
    'projdist',1,[-10,10],...
    'projfrac',0.5,[-2,2],...
    'projdist_avg',[],[],...
    'projfrac_avg',[],[],...
    'mask_midbrain_flag',false,[false true],...
    'surfname','white',[],...
    'ico_flag',true,[false true],...
    'ico_order',3,[0 7],...
    'ico_trunc_flag',false,[false true],...
    'ico_presmooth',64,[0,1000],...
... % for fs_tal2annot
    'surf_sm',100,[10,1000],...
    'tal_thresh',0.02,[0.01,0.5],...
    'tal_verbose',false,[false true],...
... % for extracting ROI time courses
    'roi_outfix',[],[],...
    'aparc_flag',true,[false true],...
    'fnames_aparc',[],[],...
    'fparc_flag',true,[false true],...
    'fnames_fparc',[],[],...
    'fname_points',[],[],...
    'aseg_flag',true,[false true],...
    'aseg_erode_flag',false,[false true],...
    'aseg_erode_nvoxels',1,[1:100],...
    'fname_aseg',[],[],...
    'csv_tseries_flag',true,[false true],...
    'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79],[1,Inf],...
... % for correlation calculations
    'corr_outfix',[],[],...
    'corr_seed_roinames',[],[],...
    'corr_roi_flag',true,[false true],...
    'corr_roi_ico_flag',true,[false true],...
    'corr_roi_ico_maps_flag',false,[false true],...
    'corr_ico_flag',true,[false true],...
    'censor_thresh',0.2,[0.001,Inf],...
    'contig_thresh',5,[0,20],...
    'resid_censor_thresh',0.3,[0.001,Inf],...
    'resid_contig_thresh',5,[0,20],...
...
    'run_outfix',[],[],...
    'motion_radius',50,[],... % for calculating distance from angle
    'motion_absflag',true,[false true],...
    'motion_nodyflag',false,[false true],...
    'fnamestem','BOLD',[],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'outtype','mgz',{'mgh','mgz'},...
    'ico_nverts',[42 162 642 2562 10242 40962 163842],[],...
...
    'group_names',{'wm','csf','brain'},[],...
    'detrend_tags',{'skipTRs','norm_flag','detrend','normval','thresh',...
                    'fname_motion','regressors','ind_valid'},[],...
    'bandpass_tags',{'TR','bandpass_low','bandpass_high','bandpass_tb',...
                     'bandpass_zpad_flag'},[],...
    'extract_group_tags',{'outdir','outstem','fname_aseg','aseg_codes',...
                          'mask_name','regfile','subjdir',... % 'skipTRs'
                          'erode_flag','erode_nvoxels','forceflag'},[],...
    'paint_tags',{'regfile' 'projfrac_flag' 'projdist' 'projfrac'...
                  'projdist_avg' 'projfrac_avg' 'subjdir' 'surfname'...
                  'mask_midbrain_flag' 'outtype' 'forceflag'},[],...
    'extract_vol_tags',{'outdir' 'regfile' 'subjdir' 'forceflag',...
                        'fname_aseg','aseg_roilist'},[],... % 'skipTRs'
    'extract_surf_tags',{'outdir' 'regfile' 'subjdir' 'forceflag'...
                         'fnames_aparc' 'hemilist' 'fnames_surf' 'verbose'...
                         'projfrac_flag' 'projdist' 'projfrac' 'projdist_avg'...
                         'projfrac_avg' 'surfname'},[],... % 'skipTRs'
    'tal2annot_tags',{'labeldir','MNIflag','vol_sm','surf_sm','thresh',...
                      'outdir','source_subj','fname_ctab','subj','subjdir',...
                      'annotname','ico','verbose','forceflag','hemilist'},[],...
  };
  parms = mmil_args2parms(args,parms_filter);

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

  % check fname_aseg
  if parms.aseg_flag
    if isempty(parms.fname_aseg)
      parms.fname_aseg = sprintf('%s/mri/aseg.mgz',FSContainerPath);
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
          FSContainerPath,hemi);
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
  if parms.fparc_flag &&...
    (~isempty(parms.fnames_fparc) || ~isempty(parms.fname_points))
    if ~isempty(parms.fnames_fparc)
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
    if ~isempty(parms.fname_points)
      if ~exist(parms.fname_points,'file')
        fprintf('%s: ERROR: points file %s not found\n',...
          mfilename,parms.fname_points);
        errcode = 1;
        return;
      end;
    end;
  end;

  % set subj name
  [parms.subjdir,parms.subj,tmp] = fileparts(FSContainerPath);
  parms.subj = [parms.subj tmp];

  % set output file stem
  if mmil_isrelative(parms.outdir)
    parms.outdir = [parms.ContainerPath '/' parms.outdir];
  end;
  parms.outstem = [parms.outdir '/' parms.outstem];
  if ~isempty(parms.roi_outfix)
    parms.roi_outstem = [parms.outstem '_' parms.roi_outfix];
  else
    parms.roi_outstem = parms.outstem;
  end;
  if ~isempty(parms.corr_outfix)
    parms.corr_outstem = [parms.roi_outstem '_' parms.corr_outfix];
  else
    parms.corr_outstem = parms.roi_outstem;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = locate_files(parms)
  errcode = 0;
  % load scan info
  [parms.ScanInfo,parms.SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(...
    parms.ContainerPath,'snums',parms.snums_valid,...
    'revflag',parms.revflag,'fnamestem',parms.fnamestem);
  if errcode, return; end;
  if isempty(parms.ScanInfo)
    fprintf('%s: ERROR: no valid BOLD data found in %s\n',...
      mfilename,parms.ContainerPath);
    errcode = 1;
    return;
  end;
  if isempty(parms.snums)
    if isempty(parms.snums_valid)
      parms.snums = [1:parms.SessInfo.nscans];
    else
      parms.snums = parms.snums_valid;
    end;
  else
    ind_bad_snums = find(parms.snums<1 | parms.snums>parms.SessInfo.nscans);
    if ~isempty(ind_bad_snums)
      bad_snums = parms.snums(ind_bad_snums);
      fprintf('%s: ERROR: bad %s snums ( %s)\n',...
        mfilename,parms.fnamestem,sprintf('%d ',bad_snums));
      errcode = 1;
      return;
    end;  
  end;
  parms.nscans = length(parms.snums);

  % set regfile
  if parms.verbose==2
    fprintf('%s: using scan %d as reference\n',...
      mfilename,parms.SessInfo.regT1_ref);
  end;
  % set fstem_ref
  fstem_ref = [parms.ContainerPath '/' parms.SessInfo.fstem_regT1_ref];
  if ~isempty(parms.infix), fstem_ref = [fstem_ref '_' parms.infix]; end;
  parms.regfile = [fstem_ref '_register.dat'];
  if ~exist(parms.regfile,'file')
    fprintf('%s: ERROR: registration file %s not found\n',...
      mfilename,parms.regfile);
    errcode = 1;
    return;
  end;
  
  % get file names and number of TRs for each BOLD scan
  parms.fnames_BOLD = cell(1,parms.nscans);
  parms.fnames_motion = cell(1,parms.nscans);
  parms.TRs = zeros(1,parms.nscans);
  parms.nreps = zeros(1,parms.nscans);
  for i=1:parms.nscans
    s = parms.snums(i);
    fstem = [parms.ContainerPath '/' parms.ScanInfo(s).fstem];
    if ~isempty(parms.infix)
      fstem = [fstem '_' parms.infix];
    end;

    parms.nreps(i) = parms.ScanInfo(s).nreps;
    if parms.nreps(i)<parms.skipTRs+parms.min_numTRs
      if parms.verbose
        fprintf('%s: WARNING: skipping %s scan %d with only %d TRs\n',...
          mfilename,parms.fnamestem,s,parms.nreps(i));
      end;
      continue;  
    end;
    parms.TRs(i) = parms.ScanInfo(s).TR/1000;

    % set fname_BOLD
    fname_BOLD = [fstem '.mgz'];
    if ~exist(fname_BOLD,'file')
      fprintf('%s: ERROR: file %s not found\n',mfilename,fname_BOLD);
      errcode = 1;
      return;
    end;
    parms.fnames_BOLD{i} = fname_BOLD;

    % set fname_motion
    if parms.mc_flag
      fname_motion = [fstem '_motion.1D'];
      if ~exist(fname_motion,'file')
        fprintf('%s: ERROR: motion 1D file %s not found\n',...
          mfilename,fname_motion);
        errcode = 1;
        return;
      end;
    else
      fname_motion = [];
    end;
    parms.fnames_motion{i} = fname_motion;
  end;

  % exclude if file name is empty
  ind_keep = find(~cellfun(@isempty,parms.fnames_BOLD));
  if isempty(ind_keep)  
    % quit if no valid scans
    fprintf('%s: ERROR: no valid BOLD scans found in %s\n',...
      mfilename,parms.ContainerPath);
    errcode = 1;
    return;
  elseif length(ind_keep) > parms.max_nscans
    % limit to first max_nscans
    fprintf('%s: WARNING: using first %d of %d scans...\n',...
      mfilename,parms.max_nscans,length(ind_keep));
    ind_keep = ind_keep(1:parms.max_nscans);
  end;
  parms.fnames_BOLD = parms.fnames_BOLD(ind_keep);
  parms.fnames_motion = parms.fnames_motion(ind_keep);
  parms.TRs = parms.TRs(ind_keep);
  parms.nreps = parms.nreps(ind_keep);
  parms.snums = parms.snums(ind_keep);
  parms.nscans = length(parms.snums);

  % get number of TRs for each scan
  parms.numTRs = zeros(size(parms.TRs));
  for i=1:parms.nscans
    mmil_mkdir(parms.outdir);
    [M,volsz] = mmil_load_mgh_info(parms.fnames_BOLD{i},parms.forceflag,parms.outdir);
    parms.numTRs(i) = volsz(4);
    % adjust TR if numTRs ~= nreps (i.e. for APE sequence)
    if parms.numTRs(i) ~= parms.nreps(i)
      parms.TRs(i) = parms.TRs(i) * parms.nreps(i) / parms.numTRs(i);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_fparc_annot(parms)
  i = parms.num_aparcs;
  for f=1:length(parms.fnames_fparc)
    fname_in = parms.fnames_fparc{f};
    if parms.verbose==2
      fprintf('%s: resampling annotation file %s from fsaverage to %s...\n',...
        mfilename,fname_in,parms.subj);
    end;
    % call fs_annot2annot
    tmp_parms = [];
%    tmp_parms.outdir = [parms.outdir '/fparc_annot'];
    tmp_parms.outdir = [parms.ContainerPath '/fparc_annot'];
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
    i = i + 1;
    parms.fnames_aparc{i} = fname_out;
    parms.aparc_hemis{i} = n.hemi;
    parms.aparc_names{i} = n.name;
  end;
  parms.num_aparcs = length(parms.fnames_aparc);
  parms.aparc_flag = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_points_annot(parms)
  if parms.verbose==2
    fprintf('%s: creating annotation files from talairach points in %s...\n',...
      mfilename,parms.fname_points);
  end;
  % create annot files from talairach points
  tmp_parms = parms;
  tmp_parms.thresh = parms.tal_thresh;
  tmp_parms.verbose = parms.tal_verbose;
  [tmp,outfix] = fileparts(parms.fname_points);
  tmp_parms.outdir = [parms.outdir '/talrc_annot'];
  tmp_parms.labeldir = sprintf('%s/fsaverage_labels',tmp_parms.outdir);
  tmp_parms.fname_ctab = sprintf('%s/%s.ctab',tmp_parms.outdir,outfix);
  args = mmil_parms2args(tmp_parms,parms.tal2annot_tags);
  [fnames_annot,fnames_label] = fs_tal2annot(parms.fname_points,args{:});  
  % add new annot files to parms.fnames_aparc
  i = parms.num_aparcs;
  for f=1:length(fnames_annot)
    fname = fnames_annot{f};
    n = regexp(fname,'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
    if isempty(n)
      error('unexpected naming convention for aparc file %s\n',fname);
    end;
    i = i + 1;
    parms.fnames_aparc{i} = fname;
    parms.aparc_hemis{i} = n.hemi;
    parms.aparc_names{i} = n.name;
  end;
  parms.num_aparcs = length(parms.fnames_aparc);
  parms.aparc_flag = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = extract_group_tseries(parms,group_name)
  fname_tag = ['fnames_' group_name '_tseries'];
  parms.(fname_tag) = cell(parms.nscans,1);
  for i=1:parms.nscans
    tparms = parms;
    tparms.outstem = sprintf('%s_data_scan%d',parms.outstem,i);
    fname_tseries = sprintf('%s_%s_tseries.mat',...
      tparms.outstem,group_name);
    if ~exist(fname_tseries,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: extracting %s time series for scan %d...\n',...
          mfilename,group_name,i);
      end;
      tparms.aseg_codes = parms.([group_name '_codes']);
      tparms.mask_name = [group_name '_mask'];
      tparms.erode_flag = parms.([group_name '_erode_flag']);
      tparms.erode_nvoxels = parms.([group_name '_erode_nvoxels']);
      args = mmil_parms2args(tparms,parms.extract_group_tags);
      % extract tseries
      tseries = mmil_extract_aseg_group_tseries(parms.subj,...
        parms.fnames_BOLD{i},args{:});
      % normalize and remove mean
      if parms.norm_flag
        meanval = mean(tseries);
        tseries = parms.normval*(tseries - meanval)/(meanval+eps);
      end;
      var_name = [group_name '_tseries'];
      eval([var_name ' = tseries;']);
      save(fname_tseries,var_name,'-v7.3');
    end;
    parms.(fname_tag){i} = fname_tseries;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_mean(parms)
  parms.fnames_BOLD_mean = cell(size(parms.fnames_BOLD));
  for i=1:parms.nscans
    fname_BOLD = parms.fnames_BOLD{i};
    [tmp,fstem] = fileparts(fname_BOLD);
    fname_BOLD_mean = sprintf('%s_data_scan%d_mean.mgz',...
      parms.outstem,i);
    if ~exist(fname_BOLD_mean,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: calculating mean of %s...\n',mfilename,fname_BOLD);
      end;
      [vol,M,tmp,volsz] = fs_load_mgh(fname_BOLD);
      vol = mean(vol,4);
      fs_save_mgh(vol,fname_BOLD_mean,M);
    end;
    parms.fnames_BOLD_mean{i} = fname_BOLD_mean;
  end;
  % average mean images
  parms.fname_BOLD_mean = [parms.outstem '_mean.mgz'];
  if ~exist(parms.fname_BOLD_mean,'file') || parms.forceflag
    if parms.nscans==1
      if exist(parms.fname_BOLD_mean,'file'), delete(parms.fname_BOLD_mean); end;
      cmd = sprintf('ln -s %s %s',...
        parms.fnames_BOLD_mean{1},parms.fname_BOLD_mean);
      [s,r] = mmil_unix(cmd);
      if s, error('command "%s" failed:\n%s',cmd,r); end;
    else
      if parms.verbose==2
        fprintf('%s: average %s mean volume...\n',mfilename,parms.fnamestem);
      end;
      mmil_concat(parms.fnames_BOLD_mean,parms.fname_BOLD_mean,...
        'meanflag',1,'forceflag',parms.forceflag);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_data(parms)
  parms.fnames_BOLD_orig = parms.fnames_BOLD;
  if parms.mask_flag ||...
     parms.mc_flag==2 || parms.wm_flag==2 ||...
     parms.csf_flag==2 || parms.brain_flag==2 ||...
     parms.norm_flag || parms.detrend || ...
     parms.skipTRs>0 || parms.bandpass_flag
    if parms.mask_flag
      % create brain mask from mean image based on intensity
      parms.fname_BOLD_mask = [parms.outstem '_mask.mgz'];
      if ~exist(parms.fname_BOLD_mask) || parms.forceflag
        if parms.verbose==2
          fprintf('%s: creating mask from %s mean volume...\n',...
            mfilename,parms.fnamestem);
        end;
        mmil_quick_brainmask([],...
          'fname_in',parms.fname_BOLD_mean,...
          'fname_out',parms.fname_BOLD_mask,...
          'thresh',parms.mask_thresh,...
          'fill1_smooth1',0,...
          'clip_edges_flag',0,...
          'fill2_smooth1',0,...
          'forceflag',parms.forceflag);
      end;
    end;
    vol_mask = [];
    for i=1:parms.nscans
      fname_BOLD = parms.fnames_BOLD{i};
      [tmp,fstem] = fileparts(fname_BOLD);
      fname_BOLD_prep = sprintf('%s_data_scan%d.mgz',...
        parms.outstem,i);
      if ~exist(fname_BOLD_prep,'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: preprocessing %s...\n',mfilename,fname_BOLD);
        end;
        [vol,M,tmp,volsz] = fs_load_mgh(fname_BOLD);
        if parms.mask_flag
          % apply mask to timeseries volume
          if isempty(vol_mask)
            vol_mask = fs_load_mgh(parms.fname_BOLD_mask);
          end;
          for f=1:volsz(4)
            vol(:,:,:,f) = vol_mask .* vol(:,:,:,f);
          end;
        end;
        tmp_parms = parms;
        tmp_parms.regressors = [];
        tmp_parms.fname_motion = [];
        % add motion estimates as regressors
        if parms.mc_flag==2
          motion_data = mmil_load_motion_1D(parms.fnames_motion{i},...
            'nframes',parms.numTRs(i));
          if parms.mc_resp_flag
            motion_data = abcd_resp_filter_motion(motion_data,parms.TRs(i),...
              parms.resp_low,parms.resp_high);
          end;
          tmp_parms.regressors = cat(2,tmp_parms.regressors,motion_data);
          if parms.deriv_flag
            motion_deriv = cat(1,zeros(1,6),diff(motion_data));
            tmp_parms.regressors = cat(2,tmp_parms.regressors,motion_deriv);
          end;
          % set valid frames based on framewise displacement
          censor_tseries = mmil_censor_motion(motion_data,...
            parms.resid_censor_thresh,parms.resid_contig_thresh,...
            parms.motion_radius,parms.motion_absflag,parms.motion_nodyflag);
          tmp_parms.ind_valid = find(~censor_tseries);
        end;
        % add aseg group ROI tseries as regressors
        for g=1:length(parms.group_names)
          group_name = parms.group_names{g};
          if parms.([group_name '_flag']) == 2
            % load tseries file
            tmp = load(parms.(['fnames_' group_name '_tseries']){i});
            tseries = tmp.([group_name '_tseries']);  
            % add to regressors
            tmp_parms.regressors = cat(2,tmp_parms.regressors,tseries);
            if parms.deriv_flag
              % calculate derivative
              deriv = [0;diff(tseries)];
              % add to regressors
              tmp_parms.regressors = cat(2,tmp_parms.regressors,deriv);
            end;
          end;
        end;
        args = mmil_parms2args(tmp_parms,parms.detrend_tags);
        vol = mmil_detrend_vol(vol,args{:});
        if parms.bandpass_flag
          parms.TR = parms.TRs(i);
          args = mmil_parms2args(parms,parms.bandpass_tags);
          vol = mmil_bandpass_vol(vol,args{:});
        end;
        fs_save_mgh(vol,fname_BOLD_prep,M);
      end;
      parms.fnames_BOLD{i} = fname_BOLD_prep;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = save_info(parms)
  parms.fname_info = sprintf('%s_info.mat',parms.outstem);
  if ~exist(parms.fname_info,'file') || parms.forceflag
    info = [];
    info.ScanInfo = parms.ScanInfo;
    info.SessInfo = parms.SessInfo;
    info.fnames_motion = parms.fnames_motion;
    info.fnames_BOLD_orig = parms.fnames_BOLD_orig;
    info.fnames_BOLD = parms.fnames_BOLD;
    info.TRs = parms.TRs;
    info.nreps = parms.nreps;   
    info.numTRs = parms.numTRs;
    % NOTE: nreps and numTRs may be identical
    %       except for APE sequence, where numTRs = nreps/2    
    save(parms.fname_info,'info');
  end;
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_motion_stats(parms)
  if parms.mc_flag
    parms.fnames_motion_stats = cell(1,parms.nscans);
    for i=1:parms.nscans
      parms.fnames_motion_stats{i} = ...
        sprintf('%s_motion_scan%d.mat',parms.outstem,i);
      if ~exist(parms.fnames_motion_stats{i},'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: calculating motion stats for scan %d...\n',mfilename,i);
        end;
        % load motion time courses
        motion_data_orig = mmil_load_motion_1D(parms.fnames_motion{i},...
          'skipTRs',parms.skipTRs,'nframes',parms.numTRs(i));
        % calculate difference in head position or rotation
        [motion_stats_orig,motion_fd_orig] = mmil_motion_stats(motion_data_orig,...
          parms.motion_radius,parms.motion_absflag,...
          parms.censor_thresh,parms.motion_nodyflag);
        % resp filter motion estimates to remove respiration noise
        if parms.mc_resp_flag
          motion_data = abcd_resp_filter_motion(motion_data_orig,parms.TRs(i),...
            parms.resp_low,parms.resp_high);
          % calculate difference in head position or rotation
          [motion_stats,motion_fd] = mmil_motion_stats(motion_data,...
            parms.motion_radius,parms.motion_absflag,...
            parms.censor_thresh,parms.motion_nodyflag);
          % normalize motion_fd
          motion_fd = motion_fd * (sum(motion_fd_orig)/sum(motion_fd));
        else
          motion_data = motion_data_orig;
          motion_stats = motion_stats_orig;
          motion_fid = motion_stats_fd;
        end;
        motion_stats.motion_data_raw = motion_data;
        % create time series of censor flags
        censor_tseries = mmil_censor_motion(motion_data,...
          parms.censor_thresh,parms.contig_thresh,...
          parms.motion_radius,parms.motion_absflag,parms.motion_nodyflag);
        % bandpass filter motion estimates to match data
        if parms.mc_bandpass_flag && parms.bandpass_flag
          motion_data = ts_freq_filt(motion_data,1/parms.TRs(i),...
            [parms.bandpass_low,parms.bandpass_high],...
            [parms.bandpass_tb,parms.bandpass_tb],'bandpass',...
            parms.bandpass_zpad_flag);
        end;
        % save to mat file
        motion_stats.motion_data = motion_data;
        motion_stats.censor_tseries = censor_tseries;
        save(parms.fnames_motion_stats{i},'-struct','motion_stats','-v7.3');
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = concat_data(parms)
  % concatenate BOLD data files
  parms.fname_BOLD = [parms.outstem '_data.mgz'];
  if ~exist(parms.fname_BOLD,'file') || parms.forceflag
    if parms.nscans==1
      if exist(parms.fname_BOLD,'file'), delete(parms.fname_BOLD); end;
      cmd = sprintf('ln -s %s %s',parms.fnames_BOLD{1},parms.fname_BOLD);
      [s,r] = mmil_unix(cmd);
      if s, error('command "%s" failed:\n%s',cmd,r); end;
    else
      if parms.verbose==2
        fprintf('%s: concatenating %s data...\n',mfilename,parms.fnamestem);
      end;
      mmil_concat(parms.fnames_BOLD,parms.fname_BOLD,'forceflag',parms.forceflag);
    end;
  end;

  % concat ROI group tseries, bandpass filter
  for g=1:length(parms.group_names)
    group_name = parms.group_names{g};
    if parms.([group_name '_flag'])
      fname_tag = ['fname_' group_name '_tseries'];
      parms.(fname_tag) =...
        [parms.outstem '_' group_name '_tseries.mat'];
      if ~exist(parms.(fname_tag),'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: concatenating %s time series data...\n',...
            mfilename,group_name);
        end;
        tseries_concat = [];
        for i=1:parms.nscans
          % load tseries file
          tmp = load(parms.(['fnames_' group_name '_tseries']){i});
          tseries = tmp.([group_name '_tseries']);  
          % bandpass filter to match data
          if parms.bandpass_flag
            tseries = ts_freq_filt(tseries,1/parms.TRs(i),...
              [parms.bandpass_low,parms.bandpass_high],...
              [parms.bandpass_tb,parms.bandpass_tb],'bandpass',...
              parms.bandpass_zpad_flag);
          end;
          % concatenate
          tseries_concat = cat(1,tseries_concat,tseries);
        end;
        % save file  
        var_name = [group_name '_tseries'];
        eval([var_name ' = tseries_concat;']);
        save(parms.(fname_tag),var_name);
      end;  
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = concat_motion_data(parms)
  % concatenate motion data
  if parms.mc_flag
    parms.fname_motion = [parms.outstem '_motion.mat'];
    if ~exist(parms.fname_motion,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: concatenating motion data...\n',mfilename);
      end;

      motion_data_raw = zeros(0,6);
      motion_data = zeros(0,6);
      censor_tseries = zeros(0,1);
      max_rz = 0; max_rx = 0; max_ry = 0;
      max_dz = 0; max_dy = 0; max_dx = 0;
      mean_trans = 0; max_trans = 0; mean_rot = 0; max_rot = 0;
      mean_motion = 0; mode_motion = 0; med_motion = 0;
      min_motion = Inf; max_motion = 0; mean_accel = 0;
      subthresh_nvols = 0; subthresh_perc = 0;
      thresholds = 0; nvols = 0;

      for i=1:parms.nscans
        % load pre-saved motion stats mat file
        motion_stats = load(parms.fnames_motion_stats{i});
        motion_data_raw = cat(1,motion_data_raw,motion_stats.motion_data_raw);
        motion_data = cat(1,motion_data,motion_stats.motion_data);
        censor_tseries = cat(1,censor_tseries,motion_stats.censor_tseries);
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
      end;

      if parms.nscans > 1
        mean_trans = mean_trans / parms.nscans;
        mean_rot = mean_rot / parms.nscans;
        mean_motion = mean_motion / parms.nscans;
        mode_motion = mode_motion / parms.nscans;
        med_motion = med_motion / parms.nscans;
        mean_accel = mean_accel / parms.nscans;
      end;

      save(parms.fname_motion,...
        'motion_data_raw',...
        'motion_data','censor_tseries',...
        'max_dx','max_dy','max_dz',...
        'max_rx','max_ry','max_rz',...
        'mean_trans','max_trans',...
        'mean_rot','max_rot',...
        'mean_motion','mode_motion','med_motion',...
        'min_motion','max_motion','mean_accel',...
        'subthresh_nvols','subthresh_perc','thresholds','nvols','-v7.3');
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = paint_data(parms)
  % sample 4D data to surface
  [tmp,fstem] = fileparts(parms.fname_BOLD);
  outstem = sprintf('%s/%s',parms.outdir,fstem);
  args = mmil_parms2args(parms,parms.paint_tags);
  parms.fname_BOLD_surf = fs_paint(parms.subj,parms.fname_BOLD,...
    'outstem',outstem,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = extract_surf_tseries(parms)
  parms.fname_surf_data = sprintf('%s%s_surf_roi_data.mat',...
    parms.roi_outstem,parms.run_outfix);
  if ~exist(parms.fname_surf_data,'file') || parms.forceflag
    % extract time series for aparc ROIs
    data = [];
    roinames = {};
    for i=1:parms.num_aparcs
      hemi = parms.aparc_hemis{i};
      h = find(strcmp(hemi,parms.hemilist));
      tmp_parms = parms;
      tmp_parms.hemilist = parms.hemilist(h);
      tmp_parms.fnames_surf = parms.fname_BOLD_surf(h);
      tmp_parms.fnames_aparc = parms.fnames_aparc(i);
%      tmp_parms.verbose = (parms.verbose==2);
      tmp_parms.verbose = 0;
      if parms.verbose==2
        fprintf('%s: extracting time series from %s for surface ROIs from %s...\n',...
          mfilename,parms.fname_BOLD_surf{h},parms.fnames_aparc{i});
      end;
      args = mmil_parms2args(tmp_parms,parms.extract_surf_tags);
      [tmp_data,tmp_roinames] = mmil_extract_aparc_tseries(parms.subj,[],args{:});
      % exclude ROIs with name 'unknown'
      ind_known = find(cellfun(@isempty,regexp(tmp_roinames,'unknown')));
      tmp_data = tmp_data(:,ind_known);
      tmp_roinames = tmp_roinames(ind_known);
      % combine data for multiple parcellations
      data = cat(2,data,tmp_data);
      roinames = cat(2,roinames,tmp_roinames);
    end;
    save(parms.fname_surf_data,'data','roinames','-v7.3');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = extract_vol_tseries(parms)
  parms.fname_vol_data = sprintf('%s%s_vol_roi_data.mat',...
    parms.roi_outstem,parms.run_outfix);
  if ~exist(parms.fname_vol_data,'file') || parms.forceflag
    % extract time series for aseg ROIs
    args = mmil_parms2args(parms,parms.extract_vol_tags);
    if parms.verbose==2
      fprintf('%s: extracting volume ROI time series from %s...\n',...
        mfilename,parms.fname_BOLD);
    end;
    [data,roinames] = mmil_extract_aseg_tseries(parms.subj,...
      parms.fname_BOLD,'erode_flag',parms.aseg_erode_flag,...
      'erode_nvoxels',parms.aseg_erode_nvoxels,args{:});
    save(parms.fname_vol_data,'data','roinames','-v7.3');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = combine_roi_tseries(parms)
  parms.fname_roi_data = sprintf('%s%s_roi_data.mat',...
    parms.roi_outstem,parms.run_outfix);
  if ~exist(parms.fname_roi_data,'file') || parms.forceflag
    % write mat file containing time series for each ROI  
    if parms.aparc_flag
      parms.fname_surf_data = sprintf('%s%s_surf_roi_data.mat',...
        parms.roi_outstem,parms.run_outfix);
      surf_data = load(parms.fname_surf_data);
    end;
    if parms.aseg_flag
      parms.fname_vol_data = sprintf('%s%s_vol_roi_data.mat',...
        parms.roi_outstem,parms.run_outfix);
      vol_data = load(parms.fname_vol_data);
    end;
    if parms.aparc_flag && parms.aseg_flag
      % combine lh and rh aparc and aseg data and create column labels
      data = cat(2,surf_data.data,vol_data.data);
      roinames = cat(2,surf_data.roinames,vol_data.roinames);
    elseif parms.aparc_flag
      data = surf_data.data;
      roinames = surf_data.roinames;
    else
      data = vol_data.data;
      roinames = vol_data.roinames;
    end;
    save(parms.fname_roi_data,'data','roinames','-v7.3');
  else
    load(parms.fname_roi_data);
  end;
  parms.roi_data = data;
  parms.roinames = roinames;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [timevec,timestrs] = set_time(parms,ntpoints)
  timevec = []; timestrs = [];
  % check if currently data file is for one scan or all (concatenated)
  k = find(strcmp(parms.fname_BOLD,parms.fnames_BOLD));
  if ~isempty(k)
    nscans = 1;
    numTRs = parms.numTRs(k);
    TRs = parms.TRs(k);
  else
    nscans = parms.nscans;
    numTRs = parms.numTRs;
    TRs = parms.TRs;
  end;
  % check number of data points
  ntpoints_data = sum(numTRs) - (nscans * parms.skipTRs);
  if ntpoints_data ~= ntpoints
    error('numTRs (%d) minus total number of skipped TRs (%d) does not match number of time points (%d)',...
      sum(numTRs),nscans * parms.skipTRs,ntpoints);
  end;
  % construct time vector
  timestrs = cell(1,ntpoints);
  time = zeros(1,ntpoints);
  t = 1;
  time = 0;
  for i=1:nscans
    TR = TRs(i);
    n = numTRs(i) - parms.skipTRs;
    for r=1:n
      timestrs{t} = sprintf('%0.2f',time);
      timevec(t) = time;
      time = time + TR;
      t = t + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = write_csv_tseries(parms,timestrs)
  parms.fname_csv_tseries = sprintf('%s%s_roi_tseries.csv',...
    parms.roi_outstem,parms.run_outfix);
  if ~exist(parms.fname_csv_tseries,'file') || parms.forceflag
    % write output file
    mmil_write_csv(parms.fname_csv_tseries,parms.roi_data,...
      'col_labels',parms.roinames,...
      'row_labels',timestrs','firstcol_label','time (sec)');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = extract_tseries(parms)
  % extract time series for surface ROIs
  if parms.aparc_flag
    parms = extract_surf_tseries(parms);
  end;
  % extract time series for volume ROIs
  if parms.aseg_flag
    parms = extract_vol_tseries(parms);
  end;
  % combine ROI time series data and save
  parms = combine_roi_tseries(parms);
  ntpoints = size(parms.roi_data,1);
  [timevec,timestrs] = set_time(parms,ntpoints);
  % save as csv file
  if parms.csv_tseries_flag
    parms = write_csv_tseries(parms,timestrs);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_ico(parms)
  [tmp,fstem] = fileparts(parms.fname_BOLD);
  outstem = sprintf('%s/%s',parms.outdir,fstem);
  tmp_parms = [];
  tmp_parms.trgsubj = 'ico';
  if ~parms.ico_trunc_flag
    tmp_parms.icolevel = parms.ico_order;
  else
    tmp_parms.icolevel = 7;
  end;
  tmp_parms.smooth_in = parms.ico_presmooth;
  tmp_parms.subjdir = parms.subjdir;
  tmp_parms.intype = parms.outtype;
  tmp_parms.forceflag = parms.forceflag;
%  tmp_parms.verbose = (parms.verbose==2);
  tmp_parms.verbose = 0;
  for h=1:parms.nhemi
    tmp_parms.hemi = parms.hemilist{h};
    tmp_parms.fname_out = sprintf('%s-sm%d-ico%d-%s.%s',...
      outstem,parms.ico_presmooth,tmp_parms.icolevel,...
      tmp_parms.hemi,parms.outtype);
    parms.fname_BOLD_ico{h} = tmp_parms.fname_out;
    if ~exist(tmp_parms.fname_out,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: resampling %s to ico...\n',...
          mfilename,parms.fname_BOLD_surf{h});
      end;
      args = mmil_parms2args(tmp_parms);
      fs_surf2surf(parms.fname_BOLD_surf{h},parms.subj,args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_variance(parms)
  % initialize vector of frame indices to use in analysis
  ind_valid = []; ntpoints = -1;

  % calculate variance across time for ROIs
  if parms.aparc_flag || parms.aseg_flag
    fname_out = sprintf('%s%s_roi_data_var.mat',...
      parms.roi_outstem,parms.run_outfix);
    if ~exist(fname_out,'file') || parms.forceflag
      if ntpoints<0
        % use censor_tseries to identify valid frames
        [ind_valid,ntpoints] = get_valid_frames(parms);
      end;
      [roi_data,roinames] = get_roi_data(parms);
      roi_data = roi_data(ind_valid,:);
      roi_var = var(roi_data,1,1);
      save(fname_out,'roi_var','roinames','ntpoints','-v7.3');
    end;
  end;
  
  % calculate variance across time for ico vertices
  if parms.ico_flag
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      fname_out = sprintf('%s%s-sm%d-ico%d-var-%s.%s',...
        parms.outstem,parms.run_outfix,...
        parms.ico_presmooth,parms.ico_order,hemi,parms.outtype);
      if ~exist(fname_out,'file') || parms.forceflag
        if ntpoints<0
          % use censor_tseries to identify valid frames
          [ind_valid,ntpoints] = get_valid_frames(parms);
        end;
        fname_in = parms.fname_BOLD_ico{h};
        vals = fs_load_mgh(fname_in,[],ind_valid);
        vals_var = var(vals,1,4);
        fs_save_mgh(vals_var,fname_out);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_valid,ntpoints] = get_valid_frames(parms)
  motion_stats = load(parms.fname_motion);
  censor_tseries = motion_stats.censor_tseries;
  ind_valid = find(~censor_tseries);
  ntpoints = length(ind_valid);
return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_correlations(parms)
  % initialize vector of frame indices to use in analysis
  ind_valid = []; ntpoints = -1;
  % initialize data
  seed_data = []; roi_data = []; ico_data = [];
  % calculate correlation matrix between seed ROIs and all ROIs
  if parms.corr_roi_flag && (parms.aparc_flag || parms.aseg_flag)
    parms.fname_roi_corr = sprintf('%s%s_roi_corr.mat',...
      parms.corr_outstem,parms.run_outfix);
    if ~exist(parms.fname_roi_corr,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: calculating correlation between ROI time courses...\n',...
          mfilename);
      end;
      if ntpoints<0
        % use censor_tseries to identify valid frames
        [ind_valid,ntpoints] = get_valid_frames(parms);
      end;
      [roi_data,roinames,time] = get_roi_data(parms);
      roi_data = roi_data(ind_valid,:);
      time = time(ind_valid);
      [seed_data,seed_roinames,seed_roinames_full] =...
                                   select_seed_data(parms,roi_data,roinames);
      if isempty(ind_valid)
        R = []; Z = [];
      else
        % calculate cross-correlation matrix
        R = corr(roi_data,seed_data);
        % Fisher's r to z transformation
        Z = atanh(R);
      end;
      % save results
      save(parms.fname_roi_corr,...
        'seed_data','seed_roinames','seed_roinames_full',...
        'roi_data','roinames','time','R','Z','ntpoints','-v7.3');
    end;
    % plot correlation image
    plot_corr_image(parms);
    % save correlation matrix to csv file
    save_corr_csv(parms);
  end;
  
  % calculate cross-correlation matrix between ico vertices
  if parms.corr_ico_flag && parms.ico_flag
    parms.fname_ico_corr = sprintf('%s%s_ico%d_corr.mat',...
      parms.corr_outstem,parms.run_outfix,parms.ico_order);
    if ~exist(parms.fname_ico_corr,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: calculating correlation between ico vertices...\n',...
          mfilename);
      end;
      if ntpoints<0
        % use censor_tseries to identify valid frames
        [ind_valid,ntpoints] = get_valid_frames(parms);
      end;
      [ico_data,hemis,vertices,time] = get_ico_data(parms);
      ico_data = ico_data(ind_valid,:);
      time = time(ind_valid);
      if isempty(ind_valid)
        R = []; Z = [];
      else
        % calculate cross-correlation matrix
        R = corr(ico_data);
        % Fisher's r to z transformation
        Z = atanh(R);
      end;
      % save results
      save(parms.fname_ico_corr,...
        'ico_data','hemis','vertices','time','R','Z','ntpoints','-v7.3');
    end;  
  end;
  
  % calculate correlation matrix between seed ROIs and ico vertices
  if parms.corr_roi_ico_flag && parms.ico_flag &&...
     (parms.aparc_flag || parms.aseg_flag)
     parms.fname_roi_ico_corr = sprintf('%s%s_roi_ico%d_corr.mat',...
      parms.corr_outstem,parms.run_outfix,parms.ico_order);
    if ~exist(parms.fname_roi_ico_corr,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: calculating correlation between ROI and ico vertex time courses...\n',...
          mfilename);
      end;
      if ntpoints<0
        % use censor_tseries to identify valid frames
        [ind_valid,ntpoints] = get_valid_frames(parms);
      end;
      if isempty(ico_data)
        [ico_data,hemis,vertices,time] = get_ico_data(parms);
        ico_data = ico_data(ind_valid,:);
        time = time(ind_valid);
      end;
      if isempty(seed_data)
        [roi_data,roinames,time] = get_roi_data(parms);
        roi_data = roi_data(ind_valid,:);
        time = time(ind_valid);
        [seed_data,seed_roinames,seed_roinames_full] =...
                                     select_seed_data(parms,roi_data,roinames);
      end;
      if isempty(ind_valid)
        R = []; Z = [];
      else
        % calculate correlations
        R = corr(ico_data,seed_data);
        % Fisher's r to z transformation
        Z = atanh(R);
      end;
      save(parms.fname_roi_ico_corr,...
        'seed_data','seed_roinames','seed_roinames_full',...
        'ico_data','hemis','vertices','time','R','Z','ntpoints','-v7.3');
    end;
    % save maps of seed ROIs to ico vertices as individual maps for each seed
    if parms.corr_roi_ico_maps_flag
      parms = save_corr_roi_ico_maps(parms);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_corr_image(parms)
  fname_plot = sprintf('%s%s_roi_corr.tif',...
    parms.corr_outstem,parms.run_outfix);
  if ~exist(fname_plot,'file') || parms.forceflag
    load(parms.fname_roi_corr);
    figure;
    imagesc(R,[-1,1]);
    colormap mmil_cmap_blueblackred;
    axis off;
    set(gcf,'visible','off');
    print(gcf,'-dtiff',fname_plot);
    close(gcf);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_corr_csv(parms)
  fname_out = sprintf('%s%s_roi_corr.csv',...
    parms.corr_outstem,parms.run_outfix);
  if ~exist(fname_out,'file') || parms.forceflag
    load(parms.fname_roi_corr);
    data = num2cell(R);
    data = cat(1,roinames,data);
    data = cat(2,cat(2,{'roiname'},seed_roinames)',data);
    mmil_write_csv(fname_out,data);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_avg_variance(parms)
  % calculate pooled variance across runs for ROIs
  if parms.aparc_flag || parms.aseg_flag
    parms.run_outfix = [];
    fname_out = sprintf('%s%s_roi_data_var.mat',...
      parms.roi_outstem,parms.run_outfix);
    if ~exist(fname_out,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: calculating average ROI variance...\n',...
          mfilename);
      end;
      wt = 0; roi_var = 0; roinames = []; ntpoints = 0;
      for i=1:parms.nscans
        % load file for individual scan
        parms.run_outfix = sprintf('_scan%d',i);
        fname_in = sprintf('%s%s_roi_data_var.mat',...
          parms.roi_outstem,parms.run_outfix);
        tmp = load(fname_in);
        roinames = tmp.roinames;
        % weight by number of TRs used in calculation
        n = tmp.ntpoints;
        if n>1
          wt = wt + n-1;
          roi_var = roi_var + (n-1)*tmp.roi_var;
          ntpoints = ntpoints + n;
        end;
      end;
      % average across runs, weighted by total numTRs
      roi_var = roi_var / wt;      
      % save results
      save(fname_out,'roi_var','roinames','ntpoints','-v7.3');
    end;
  end;

  % calculate variance across time for ico vertices
  if parms.ico_flag
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      parms.run_outfix = [];
      fname_out = sprintf('%s%s-sm%d-ico%d-var-%s.%s',...
        parms.outstem,parms.run_outfix,...
        parms.ico_presmooth,parms.ico_order,hemi,parms.outtype);
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: calculating average ico variance...\n',...
            mfilename);
        end;
        wt = 0; vals_var = 0; ntpoints = 0;
        for i=1:parms.nscans
          % load file for individual scan
          parms.run_outfix = sprintf('_scan%d',i);
          fname_in = sprintf('%s%s-sm%d-ico%d-var-%s.%s',...
            parms.outstem,parms.run_outfix,...
            parms.ico_presmooth,parms.ico_order,hemi,parms.outtype);
          tmp_var = fs_load_mgh(fname_in);
          % weight by number of TRs used in calculation
          parms.motion = parms.fnames_motion_stats{i};
          [ind_valid,ntpoints] = get_valid_frames(parms);
          n = ntpoints;
          if n>1
            wt = wt + n-1;
            vals_var = vals_var + (n-1)*tmp_var;
            ntpoints = ntpoints + n;
          end;
        end;
        % average across runs, weighted by total numTRs
        vals_var = vals_var / wt;      
        fs_save_mgh(vals_var,fname_out);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calc_avg_correlations(parms)
  % calculate correlation matrix between seed ROIs and all ROIs
  if parms.corr_roi_flag && (parms.aparc_flag || parms.aseg_flag)
    parms.run_outfix = [];
    parms.fname_roi_corr = sprintf('%s%s_roi_corr.mat',...
      parms.corr_outstem,parms.run_outfix);
    if ~exist(parms.fname_roi_corr,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: calculating average correlation between ROI time courses...\n',...
          mfilename);
      end;
      seed_data = []; roi_data = []; time = []; Z = 0; ntpoints = 0;
      for i=1:parms.nscans
        % load file for individual scan
        parms.run_outfix = sprintf('_scan%d',i);
        fname_in = sprintf('%s%s_roi_corr.mat',...
          parms.corr_outstem,parms.run_outfix);
        tmp = load(fname_in);
        seed_roinames = tmp.seed_roinames;
        seed_roinames_full = tmp.seed_roinames_full;
        roinames = tmp.roinames;
        % concatenate time series data
        seed_data = cat(1,seed_data,tmp.seed_data);
        roi_data = cat(1,roi_data,tmp.roi_data);
        % concat time vector
        if i==1
          t_offset = 0;
        else
          t_offset = 2*time(end) - time(end-1);
        end;
        time = cat(2,time,tmp.time + t_offset);
        % weight by number of TRs used in calculation
        n = tmp.ntpoints;
        if n>1
          ntpoints = ntpoints + n;
          Z = Z + n*tmp.Z;
        end;
      end;
      % average across runs
      Z = Z / ntpoints;
      % z to r transformation
      R = tanh(Z);
      % save results
      save(parms.fname_roi_corr,...
        'seed_data','seed_roinames','seed_roinames_full',...
        'roi_data','roinames','time','R','Z','ntpoints','-v7.3');
    end;
    % plot correlation image
    plot_corr_image(parms);
    % save correlation matrix to text file
    save_corr_csv(parms);
  end;  
  
  % calculate cross-correlation matrix between ico vertices
  if parms.corr_ico_flag && parms.ico_flag
    parms.run_outfix = [];
    parms.fname_ico_corr = sprintf('%s%s_ico%d_corr.mat',...
      parms.corr_outstem,parms.run_outfix,parms.ico_order);
    if ~exist(parms.fname_ico_corr,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: calculating average correlation between ico vertices...\n',...
          mfilename);
      end;
      ico_data = []; time = []; Z = 0; ntpoints = 0;
      for i=1:parms.nscans
        % load file for individual scan
        parms.run_outfix = sprintf('_scan%d',i);
        fname_in = sprintf('%s%s_ico%d_corr.mat',...
          parms.corr_outstem,parms.run_outfix,parms.ico_order);
        tmp = load(fname_in);
        hemis = tmp.hemis;
        vertices = tmp.vertices;
        % concatenate time series data
        ico_data = cat(2,ico_data,tmp.ico_data);
        % concat time vector
        if i==1
          t_offset = 0;
        else
          t_offset = 2*time(end) - time(end-1);
        end;
        time = cat(2,time,tmp.time + t_offset);
        % weight by numTRs
        n = tmp.ntpoints;
        if n>1
          ntpoints = ntpoints + n;
          Z = Z + n*tmp.Z;
        end;
      end;
      % average across runs
      Z = Z / ntpoints;
      % z to r transformation
      R = tanh(Z);
      % save results
      save(parms.fname_ico_corr,...
        'ico_data','hemis','vertices','time','R','Z','ntpoints','-v7.3');
    end;
  end;
  
  % calculate correlation matrix between seed ROIs and ico vertices
  if parms.corr_roi_ico_flag && parms.ico_flag &&...
     (parms.aparc_flag || parms.aseg_flag)
    parms.run_outfix = [];
    parms.fname_roi_ico_corr = sprintf('%s%s_roi_ico%d_corr.mat',...
      parms.corr_outstem,parms.run_outfix,parms.ico_order);
    if ~exist(parms.fname_roi_ico_corr,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: calculating correlation between ROI and ico vertex time courses...\n',...
          mfilename);
      end;
      seed_data = []; ico_data = []; time = []; Z = 0; ntpoints = 0;
      for i=1:parms.nscans
        % load file for individual scan
        parms.run_outfix = sprintf('_scan%d',i);
        fname_in = sprintf('%s%s_roi_ico%d_corr.mat',...
          parms.corr_outstem,parms.run_outfix,parms.ico_order);
        tmp = load(fname_in);
        seed_roinames = tmp.seed_roinames;
        seed_roinames_full = tmp.seed_roinames_full;
        hemis = tmp.hemis;
        vertices = tmp.vertices;
        % concatenate time series data
        seed_data = cat(1,seed_data,tmp.seed_data);
        ico_data = cat(2,ico_data,tmp.ico_data);
        % concat time vector
        if i==1
          t_offset = 0;
        else
          t_offset = 2*time(end) - time(end-1);
        end;
        time = cat(2,time,tmp.time + t_offset);
        % weight by numTRs
        n = tmp.ntpoints;
        if n>1
          ntpoints = ntpoints + n;
          Z = Z + n*tmp.Z;
        end;
      end;
      % average across runs
      Z = Z / ntpoints;
      % z to r transformation
      R = tanh(Z);
      % save results
      save(parms.fname_roi_ico_corr,...
        'seed_data','seed_roinames','seed_roinames_full',...
        'ico_data','hemis','vertices','time','R','Z','ntpoints','-v7.3');
    end;
    % save maps of seed ROIs to ico vertices as individual maps for each seed
    if parms.corr_roi_ico_maps_flag
      parms.run_outfix = [];
      parms = save_corr_roi_ico_maps(parms);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = save_corr_roi_ico_maps(parms)
  % check if output files exist
  all_exist = 1;
  seed_roinames = get_seed_roinames(parms);
  for r=1:length(seed_roinames)
    seed_roiname = regexprep(seed_roinames{r},' ','_');
    for h=1:length(parms.hemilist)
      hemi = parms.hemilist{h};
      fname = sprintf('%s%s_%s_ico%d_corr-%s.%s',...
        parms.corr_outstem,parms.run_outfix,...
        seed_roiname,parms.ico_order,hemi,parms.outtype);
      if ~exist(fname,'file') || parms.forceflag
        all_exist = 0;
        break;
      end;
    end;
    if ~all_exist, break; end;
  end;
  if ~all_exist
    load(parms.fname_roi_ico_corr);
    for r=1:length(seed_roinames)
      seed_roiname = regexprep(seed_roinames{r},' ','_');
      for h=1:length(parms.hemilist)
        hemi = parms.hemilist{h};
        fname = sprintf('%s%s_%s_ico%d_corr-%s.%s',...
          parms.corr_outstem,parms.run_outfix,...
          seed_roiname,parms.ico_order,hemi,parms.outtype);
        if ~exist(fname,'file') || parms.forceflag
          ind = find(hemis==h);
          vals = Z(ind,r);
          fs_save_mgh(vals,fname);
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,roinames,time] = get_roi_data(parms)
  data = parms.roi_data;
  roinames = parms.roinames;
  ntpoints = size(parms.roi_data,1);
  time = set_time(parms,ntpoints);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [seed_data,seed_roinames,seed_roinames_full] = ...
                                    select_seed_data(parms,roi_data,roinames)
  % extract time courses for seed ROIs
  %   if more than one match, average
  if isempty(parms.corr_seed_roinames)
    seed_roinames = roinames;
    seed_roinames_full = roinames;
    seed_data = roi_data;
  else
    seed_roinames = parms.corr_seed_roinames;
    if ~iscell(seed_roinames), seed_roinames = {seed_roinames}; end;
    nseeds = length(seed_roinames);
    seed_data = zeros(size(roi_data,1),nseeds);
    seed_roinames_full = cell(1,nseeds);
    for s=1:nseeds
      seed_roiname = seed_roinames{s};
      ind_seed = find(~cellfun(@isempty,regexp(roinames,seed_roiname)));
      if isempty(ind_seed)
        error('seed_roiname %s not found in data roinames',seed_roiname);
      end;
      seed_roinames_full{s} = roinames(ind_seed);
      seed_data(:,s) = mean(roi_data(:,ind_seed),2);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seed_roinames = get_seed_roinames(parms)
  if isempty(parms.corr_seed_roinames)
    seed_roinames = parms.roinames;
  else
    seed_roinames = parms.corr_seed_roinames;
  end;
  if ~iscell(seed_roinames), seed_roinames = {seed_roinames}; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,hemis,vertices,time] = get_ico_data(parms)
  data = [];
  hemis = [];
  vertices = [];
  for h=1:parms.nhemi
    hemi_data = squeeze(fs_load_mgh(parms.fname_BOLD_ico{h}))';
    % select subset of vertices
    if parms.ico_trunc_flag && parms.ico_order~=7
      nverts = parms.ico_nverts(parms.ico_order);
      hemi_data = hemi_data(:,1:nverts);
    end;
    nverts = size(hemi_data,2);
    data = cat(2,data,hemi_data);
    hemis = cat(2,hemis,h*ones(1,nverts));
    vertices = cat(2,vertices,[1:nverts]);
  end;
  ntpoints = size(data,1);
  time = set_time(parms,ntpoints);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

