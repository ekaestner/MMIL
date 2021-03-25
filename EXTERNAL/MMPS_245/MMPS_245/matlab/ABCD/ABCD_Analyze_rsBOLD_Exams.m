function ABCD_Analyze_rsBOLD_Exams(ProjID,varargin)
% function ABCD_Analyze_rsBOLD_Exams(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject including fields:
%    SubjID, VisitID, STRUCT_VisitID
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: proc_bold, fsurf
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'fname_info': full pathname of procotol compliance info file
%     used to determine which BOLD series to include in analysis
%     if empty, will be ~/MetaData/{ProjID}/{ProjID}_pcinfo.csv
%     if not found, will include all available BOLD series in analysis
%     {default = []}
%  'SeriesType': series type according to column in fname_info
%     {default = 'rsfMRI'}
%  'qcflag': [0|1] whether to exclude subjects with StudyInfo.QC=0
%    {default = 1}
%  'batchname': name of output batchdir
%    {default = 'ABCD_Analyze_rsBOLD_Exams'}
%  'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional parameters for specifying output:
%  'outdir': output directory (relative to ContainerPath)
%    {default = 'rsBOLD_analysis'}
%  'oustem': file stem for output file names
%    {default = 'rsBOLD'}
%
% Optional parameters for data selection:
%  'snums': vector of scan numbers to analyze
%    if empty, will use all valid BOLD scans in ContainerPath
%    {default = []}
%  'snums_valid': vector of scan numbers that were processed
%    if empty, will assume all were processed
%    {default = []}
%  'infix': input file will be sprintf('BOLD%d_%s.mgz',snum,infix)
%    to specify no infix (i.e. raw data), use 'none'
%    {default = 'corr_resBOLD'}
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
%    {default: 2}
%  'mc_resp_flag': [0|1] notch filter motion estimates for aliased respiration signals
%    ignored if mc_flag = 0
%    note: this is done after nuisance regression, but before censoring
%    {default = 1}
%  'mc_bandpass_flag': [0|1] bandpass filter motion estimates
%    ignored if mc_flag = 0 or bandpass_flag = 0
%    note: this is done after being used as nuisance regressors
%    {default: 1}
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
%     cumulative probability of intensity (max is 1)
%     {default = 0.7}
%  'norm_flag': [0|1] whether to normalize input timeseries
%    by mean for each voxel (new mean = 100)
%    NOTE: if detrend>0, new mean will be 0, but values will be % change
%    {default: 1}
%  'thresh': when normalize to mean, set voxels with original values
%    less than this to 0
%    {default: 10}
%  'detrend': [0|1|2] whether and how to detrend input timeseries
%    0: no detrend (not recommended)
%    1: linear detrend
%    2: quadratic detrend
%    3: cubic detrend
%    {default: 2}
%  'skipTRs': number of TRs at beginning of scan to be ignored in time shifting
%    {default = 0}
%  'resp_low': low frequency cut-off (Hz)
%    {default = 18.6}
%  'resp_high': high frequency cut-off (Hz)
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
%  'ico_flag': resample data to icosahedral sphere
%    {default = 1}
%  'ico_order': icosahedral order (0-7)
%    {default = 3}
%  'ico_trunc_flag': [0|1] for ico_order<7, truncate extra values
%     otherwise, actually resample to ico_order
%    {default = 0}
%  'ico_presmooth': number of smoothing steps applied on native surface
%    before resampling to ico
%    NOTE: FWHM ~ 1.25 * sqrt(N)
%    {default = 64}
%
% Optional parameters for extracting ROI time courses
%  'roi_outfix': additional string added to ROI output files
%    e.g. to use existing preprocessed data and surface time courses
%         but rerun ROI extraction with different settings
%    {default = []}
%  'aparc_flag': extract time series for aparc cortical surface ROIs
%    {default = 1}
%  'fnames_aparc': cell array of annotation files (one for each hemisphere)
%     (relative to subjdir/subj/label)
%     if empty, will use ?h.aparc.annot files in subjdir/subj/label
%    {default = []}
%  'fparc_flag': extract time series for fparc cortical surface ROIs
%    ignored if fnames_fparc and fname_points are empty
%    {default = 1}
%  'fnames_fparc': cell array of annotation files in fsaverage space
%    will be resampled to individual subject space before use
%    {default = []}
%  'fname_points': name of csv file containing talairach points
%     will be sampled onto fsaverage surface and used to generate annot files
%     see fs_tal2annot
%    {default = []}
%  'aseg_flag': extract time series for aseg volume ROIs
%    {default = 1}
%  'aseg_erode_flag': [0|1] "erode" aseg ROIs by smoothing and thresholding
%    after resampling to BOLD resolution
%    {default = 0}
%  'aseg_erode_nvoxels': smoothing kernel sigma for aseg erosion (# of voxels)
%    {default = 1}
%  'fname_aseg': name of aseg file (relative to subjdir/subj/mri)
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
%     {default = 0.2}
%  'contig_thresh': mininum number of contiguous frames required (motion censoring)
%     {default = 5}
%
% Created:  07/03/17 by Don Hagler
% Prev Mod: 10/18/17 by Don Hagler
% Last Mod: 10/23/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% based on MMIL_Analyze_rsBOLD_Exams
% Created:  07/17/12 by Don Hagler
% Last Mod: 10/05/16 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchname','ABCD_Analyze_rsBOLD_Exams',[],...
  'fname_info',[],[],...
  'SeriesType','rsfMRI',[],...
  'qcflag',true,[false true],...
  'verbose',1,[0:2],...
  'forceflag',false,[false true],...
...
  'outdir','rsBOLD_analysis',[],...
  'outstem','rsBOLD',[],...
... % for data selection
  'snums',[],[],...
  'snums_valid',[],[],...
  'infix','corr_resBOLD',[],...
  'revflag',2,[0,1,2],...  
  'min_numTRs',20,[1,1000],...
  'max_nscans',[],[1,Inf],...
... % for preprocessing
  'mc_flag',2,[0:2],...
  'mc_resp_flag',true,[false true],...
  'mc_bandpass_flag',true,[false true],...
  'motion_radius',50,[],... % for calculating distance from angle
  'motion_absflag',true,[false true],...
  'motion_nodyflag',false,[false true],...
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
...
  'required_containers',{'proc_bold','fsurf'},[],...
  'modality','MRI',[],...
  'QC_raw',true,[false true],... % only applies if manual raw QC exists
  'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists
  'QC_recon',true,[false true],...
...
  'skipTRs_GE',5,[],... %% todo: must check GE software version for DV26 (skipTRs_GEdv26 = 16)
  'skipTRs_Philips',8,[],...
  'skipTRs_Siemens',8,[],...
...
  'info_tags',{'VisitIDs','SubjIDs','StudyInfo','RootDirs','ignore_VisitInfo_flag',...
               'user','numvec_tags'},[],...
  'analy_tags',{'outdir' 'outstem' 'verbose' 'forceflag' ... % for data selection
                'snums' 'snums_valid' 'infix' 'revflag' 'min_numTRs' 'max_nscans' ... % for preprocessing
                'mc_flag' 'mc_bandpass_flag' 'mc_resp_flag'...
                'mask_flag' 'mask_thresh' 'norm_flag' 'normval' 'thresh' 'detrend' 'skipTRs'...
                'resp_low' 'resp_high'...
                'bandpass_flag' 'bandpass_low' 'bandpass_high' 'bandpass_tb' 'bandpass_zpad_flag'...
                'wm_flag' 'wm_codes' 'wm_erode_flag' 'wm_erode_nvoxels' 'csf_flag' 'csf_codes' 'csf_erode_flag' 'csf_erode_nvoxels' 'brain_flag' 'brain_codes' 'brain_erode_flag' 'brain_erode_nvoxels' 'deriv_flag' ... % for sampling to surface
                'projfrac_flag' 'projdist' 'projfrac' 'projdist_avg' 'projfrac_avg' 'mask_midbrain_flag' 'surfname' 'ico_flag' 'ico_order' 'ico_trunc_flag' 'ico_presmooth' ... % for fs_tal2annot
                'surf_sm' 'tal_thresh' 'tal_verbose' ... % for extracting ROI time courses
                'roi_outfix' 'aparc_flag' 'fnames_aparc'...
                'fparc_flag' 'fnames_fparc' 'fname_points'...
                'aseg_flag' 'aseg_erode_flag' 'aseg_erode_nvoxels'...
                'fname_aseg' 'csv_tseries_flag' 'aseg_roilist' ... % for correlation calculations
                'corr_outfix' 'corr_seed_roinames' 'corr_roi_flag' 'corr_roi_ico_flag' 'corr_roi_ico_maps_flag' 'corr_ico_flag' ...
                'censor_thresh' 'contig_thresh'...
                'motion_radius' 'motion_absflag' 'motion_nodyflag'...
                'fnamestem' 'hemilist' 'outtype' 'ico_nverts' },[],...
};
parms = mmil_args2parms(varargin,parms_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check fnames_aparc has two elements
if ~isempty(parms.fnames_aparc)
  if ~iscell(parms.fnames_aparc)
    parms.fnames_aparc = {parms.fnames_aparc};
  end;
end;

if strcmp(parms.infix,'none'), parms.infix = []; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = mmil_parms2args(parms,parms.info_tags);
[ProjInfo,StudyInfo,RootDirs] = MMIL_Quick_Check_ProjID(ProjID,args{:});

if isempty(StudyInfo), error('empty StudyInfo'); end;
if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load info file to find series with the right type
if isempty(parms.fname_info)
  parms.fname_info = sprintf('%s/MetaData/%s/%s_pcinfo.csv',...
    RootDirs.home,ProjID,ProjID);
end;
if ~exist(parms.fname_info,'file')
  fprintf('%s: WARNING: file %s not found, will analyze all available BOLD series\n',...
    mfilename,parms.fname_info);
  pcinfo = [];
else
  all_pcinfo = abcd_load_csv(parms.fname_info);
  ind_type = find(strcmp({all_pcinfo.SeriesType},parms.SeriesType));
  pcinfo = all_pcinfo(ind_type);
  pcinfo_VisitIDs = {pcinfo.VisitID};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s/*\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;
for i=1:length(StudyInfo)
  VisitID = StudyInfo(i).VisitID;

  % find containers for this VisitID
  ContainerPath = MMIL_Get_Container(RootDirs,VisitID,'proc_bold');
%  RawContainerPath = MMIL_Get_Container(RootDirs,VisitID,'raw');
  FSContainerPath = MMIL_Get_Container(RootDirs,VisitID,'fsurf');

  if isempty(ContainerPath)
    fprintf('%s: WARNING: missing proc_bold container for %s\n',...
      mfilename,VisitID);
    continue;
  end;
%  if isempty(RawContainerPath)
%    fprintf('%s: WARNING: missing raw container for %s\n',...
%      mfilename,VisitID);
%    continue;
%  end;
  if isempty(FSContainerPath)
    fprintf('%s: WARNING: missing fsurf container for %s\n',...
      mfilename,VisitID);
    continue;
  end;

  % find pcinfo entries for this VisitID
  if ~isempty(pcinfo)
    ind_pcinfo = find(strcmp(VisitID,pcinfo_VisitIDs));
    if isempty(ind_pcinfo)
      fprintf('%s: skipping %s because no series with SeriesType %s\n',...
        VisitID,parms.SeriesType);
      continue;
    end;
    pcinfo_SeUIDs = {pcinfo(ind_pcinfo).SeriesInstanceUID};
  else
    pcinfo_SeUIDs = [];
  end;

  targs = {ContainerPath,FSContainerPath};
  tparms = parms;

  % set fnames_aparc and fname_aseg using relative locations
  if ~isempty(parms.fnames_aparc)
    for h=1:length(parms.fnames_aparc)
      tmp = parms.fnames_aparc{h};
      if mmil_isrelative(tmp)
        tparms.fnames_aparc{h} = sprintf('%s/label/%s',FSContainerPath,tmp);
      end;
    end;
  end;
  if ~isempty(parms.fname_aseg)
    tmp = parms.fname_aseg;
    if mmil_isrelative(tmp)
      tparms.fname_aseg = sprintf('%s/mri/%s',FSContainerPath,tmp);
    end;
  end;

  % load ContainerInfo from raw and proc containers
%  [RawContainerInfo,errcode] = MMIL_Load_ContainerInfo(RawContainerPath);
%  if errcode, continue; end;
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
  if errcode, continue; end;

  % set skipTRs according to Scanner Manufacturer
  switch ContainerInfo.Manufacturer
    case 'GE MEDICAL SYSTEMS'
      tparms.skipTRs = parms.skipTRs_GE;
    case 'Philips Medical Systems'
      tparms.skipTRs = parms.skipTRs_Philips;
    case 'SIEMENS'
      tparms.skipTRs = parms.skipTRs_Siemens;
    otherwise
      fprintf('%s: WARNING: unsupported manufacturer "%s", setting skipTRs = 0\n',...
        mfilename,ContainerInfo.Manufacturer);
      tparms.skipTRs = 0;
  end;

  if ~isempty(pcinfo_SeUIDs)
    % match scans based on SeUID
%    SeUIDs = {RawContainerInfo.SeriesInfo.SeriesInstanceUID};
    SeUIDs = {ContainerInfo.SeriesInfo.SeriesInstanceUID};
    [tmp,ind_pcinfo,ind_series] = intersect(pcinfo_SeUIDs,SeUIDs);
    if isempty(ind_series)
      fprintf('%s: WARNING: mismatch in SeriesInstanceUIDs for VisitID %s... skipping\n',...
        mfilename,VisitID);
      %% todo: change this message? just means no resting state scans, right?
%      keyboard
      continue;
    end;
    % find corresponding BOLD snums
    ind_BOLD = [ContainerInfo.ScanInfo.BOLD.SeriesIndex];
    [tmp,ind_s,ind_B] = intersect(ind_series,ind_BOLD);
    tparms.snums = ind_B;
  end;

  % set valid snums (from BOLDScanNums in VisitInfo file)
  if isfield(StudyInfo(i),'BOLDScanNums') && ~isempty(StudyInfo(i).BOLDScanNums)
    tparms.snums_valid = StudyInfo(i).BOLDScanNums;
    tparms.snums = intersect(tparms.snums,tparms.snums_valid);
  end;

  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem);
  jobfname = [batchdir '/' jobID '.m'];
  mmil_write_script(jobfname,'ABCD_Analyze_rsBOLD_Exam',...
    targs,parms.analy_tags,tparms);

  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
  j=j+1;
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmajjobs15 %s\n',parms.batchname);

