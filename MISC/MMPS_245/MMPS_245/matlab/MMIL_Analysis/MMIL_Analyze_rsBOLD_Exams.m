function MMIL_Analyze_rsBOLD_Exams(ProjID,varargin)
% function MMIL_Analyze_rsBOLD_Exams(ProjID,[options])
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
%  'qcflag': [0|1] whether to exclude subjects with StudyInfo.QC=0
%    {default = 1}
%  'batchname': name of output batchdir
%    {default = 'MMIL_Analyze_rsBOLD_Exams'}
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
%  'infix': BOLD file name infix (e.g. '', 'corr')
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
%  'mc_bandpass_flag': [0|1] bandpass filter motion estimates
%    ignored if mc_flag = 0 or bandpass_flag = 0
%    note: this is done after being used as nuisance regressors
%    {default: 1}
%  'motion_absflag': calculate frame-to-frame motion as sum of abs values
%     otherwise, calculate square root of sum of squares
%     {default = 1}
%  'mask_flag': [0|1] mask input volumes based on mean intensities
%    {default = 1}
%  'mask_thresh': relative threshold applied to mean volume
%     cumulative probability of intensity (max is 1)
%     {default = 0.9}
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
%  'bandpass_flag': [0|1] bandpass filter time series
%    {default = 1}
%  'bandpass_low': low frequency cut-off (Hz)
%    {default = 0.01}
%  'bandpass_high': high frequency cut-off (Hz)
%    {default = 0.08}
%  'bandpass_tb': bandpass filter transition band (Hz)
%    {default = 0.001}
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
% Optional parameters for correlation calculatiosn:
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
%
% Created:  07/17/12 by Don Hagler
% Last Mod: 10/05/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchname','MMIL_Analyze_rsBOLD_Exams',[],...
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
  'mc_bandpass_flag',true,[false true],...
  'motion_absflag',true,[false true],...
  'mask_flag',true,[false true],...
  'mask_thresh',0.9,[],...
  'norm_flag',true,[false true],...
  'normval',100,[1 1e10],...
  'thresh',10,[0 Inf],...
  'detrend',2,[0:2],...
  'skipTRs',0,[0,100],...
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
...
  'required_containers',{'proc_bold','fsurf'},[],...
  'modality','MRI',[],...
  'QC_raw',true,[false true],... % only applies if manual raw QC exists
  'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists
  'QC_recon',true,[false true],...
};
parms = mmil_args2parms(varargin,parms_filter);

% check fnames_aparc has two elements
if ~isempty(parms.fnames_aparc)
%  if ~iscell(parms.fnames_aparc) || length(parms.fnames_aparc)~=2
%    error('fnames_aparc must be cell array with two elements');
%  end;
  if ~iscell(parms.fnames_aparc)
    parms.fnames_aparc = {parms.fnames_aparc};
  end;
end;

% excl_tags are fields that should not be passed to BOLD_MMIL_Resting_Analysis
excl_tags = {'StudyInfo' 'RootDirs' 'batchname' 'qcflag'...
  'required_containers' 'modality' 'QC_raw' 'QC_BOLD','QC_recon'};
tags = setdiff(fieldnames(parms),excl_tags);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

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
  tmpInfo = StudyInfo(i);
  VisitID = tmpInfo.VisitID;
  ContainerPath = fullfile(RootDirs.proc_bold,StudyInfo(i).proc_bold);
  FSContainerPath = fullfile(RootDirs.fsurf,StudyInfo(i).fsurf);
  tmp_args = {ContainerPath,FSContainerPath};
  tmp_parms = parms;

  % set fnames_aparc and fname_aseg using relative locations
  if ~isempty(parms.fnames_aparc)
    for h=1:length(parms.fnames_aparc)
      tmp = parms.fnames_aparc{h};
      if mmil_isrelative(tmp)
        tmp_parms.fnames_aparc{h} = sprintf('%s/label/%s',FSContainerPath,tmp);
      end;
    end;
  end;
  if ~isempty(parms.fname_aseg)
    tmp = parms.fname_aseg;
    if mmil_isrelative(tmp)
      tmp_parms.fname_aseg = sprintf('%s/mri/%s',FSContainerPath,tmp);
    end;
  end;

  if isfield(tmpInfo,'BOLDScanNums') && ~isempty(tmpInfo.BOLDScanNums)
    tmp_parms.snums_valid = tmpInfo.BOLDScanNums;
  end;

  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem);
  jobfname = [batchdir '/' jobID '.m'];
  mmil_write_script(jobfname,'BOLD_MMIL_Resting_Analysis',...
    tmp_args,tags,tmp_parms);

  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
  j=j+1;
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs2 %s\n',parms.batchname);

