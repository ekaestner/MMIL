function results = MMIL_Summarize_MRI_DTI_Analysis(ProjID,varargin)
% function results = MMIL_Summarize_MRI_DTI_Analysis(ProjID,[options])
%
% Purpose: create summary spreadsheet with MRI and DTI ROI measures
%
% Usage:
%  results = MMIL_Summarize_MRI_DTI_Analysis(ProjID,'key1', value1,...);
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%      proc_dti, fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'QC_raw' : [0|1] use good raw QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_RawQC.mat
%     {default = 1}
%  'QC_recon' : [0|1] use recon QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_FSReconQC.csv
%     {default = 1}
%  'QC_DTI' : [0|1] use DTIreg QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_DTIQC.csv
%     {default = 1}
%  'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
%     {default = 1}
%
% Optional Parameters Controlling which MRI Stats to Compile:
%  'thick_flag': cortical thickness from aparc stats
%    {default = 1}
%  'area_flag': cortical area from aparc stats
%    {default = 1}
%  'cortT1w_flag': average T1-weighted values for aparc ROIs
%    and calculate gray-white contrast
%    {default = 1}
%  'subvol_flag': subcortical volume from aseg stats
%    {default = 1}
%  'subT1w_flag': average T1-weighted values for aseg ROIs 
%    {default = 1}  
%  'fuzzy_flag': [0|1] use fuzzy cluster weighted surface ROIs
%      for thickness, area, and T1w
%    applies if one or more of thick_flag, area_flag, or cortT1w_flag = 1
%    {default = 1}
%  'fuzzy_order': [2|4|12|18] number of fuzzy cluster ROIs
%    note: set of 18 includes combined sets of 2, 4, and 12
%    {default = 18}
%
% Optional Parameters that specify which DTI analysis results to compile:
%   'measlist': list of DT "measures" to extract for each fiber
%      e.g. 'FA', 'MD', 'LD', 'TD', 'b0', 'b0N', 'T2w', 'T1w'
%      'b0' is the b=0 volume collected with diffusion scan
%      'b0N' is b=0 volume normalized by mean signal in ventricles
%      'T2w' is b=0 volume normalized by linear fit to MD
%      'T1w' is the T1-weighted FreeSurfer nu.mgz volume
%       use 'none' to extract none of these measures
%     {default = {'FA','MD','LD','TD','T2w','T1w'}}
%   'scalefacts': scaling factors applied to each measure in measlist
%     if empty, will use 1 for each
%     if not empty, length must match length of measlist
%     {default = []}
%   'inputlist': list of input files in addition to or in place of measlist
%     relative to each ContainerPath
%     {default = []}
%   'motion_flag': [0|1] compile DTI head motion measures
%     {default = 1}
%   'fiber_flag': [0|1] compile white matter fiber ROI results
%     {default = 1}
%   'fiber_vol_flag': [0|1] calculate fiber volumes from number of
%     values extracted for first element of measlist
%     ignored if fiber_flag = 0
%     NOTE: fiber volume is only meaningful for manual fibers
%       or atlas fibers with thresh_prob > 0
%     {default = 0}
%   'aseg_flag': [0|1] compile aseg ROI results
%     {default = 1}
%   'wmparc_flag': [0|1] compile wmparc ROI results
%     {default = 0}
%   'cortsurf_flag': [0|1] compile cortical surface ROI results
%     {default = 1}
%
% Optional Parameters that determine input diffusion data:
%  'snums_flag': [0|1|2|3] which set of "snums" to use
%     0: use all available scans
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%    {default = 1}
%  'snum_index': index specifying which scan number of DTIScanNums
%     (or DTIScanNums2) to use in spread sheet (must have run DT
%     calculations separately for each scan)
%     If empty, use DT measures calculated from all DTIScanNums
%     {default = []}
%  'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix: 'corr_resDTI'
%     {default = []}
%  'auto_infix_flag': [0|1] set infix automatically based on typical
%     processing and settings in ProjInfo
%     ignored if infix is not empty
%     {default = 1}
%  'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%  'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%  'min_ndirs': minimum number of diffusion directions to be included
%     {default = 6}
%  'full_fstem_flag': [0|1] full stem included in DT measures analysis output
%     otherwise, use only meas name (e.g. 'FA' instead of 'DTI_scans_1_2_...FA')
%     if 0, options above (snums_flag, snum_index, infix, etc.) are ignored
%    {default = 0}
%  'DT_analdir': name of DTI analysis folder relative to ContainerPath
%    {default = 'DTanalysis'}
%  'DT_outfix': string attached to DT calculation and analysis file names
%    {default = []}
%
% Optional Parameters specific to DTI fiber analysis:
%   'fibers': vector of fiber numbers
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'atlas_flag': whether to use atlas fibers and if so,
%      what type of atlas fibers
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'weighted_avg_flag': [0|1] use weighted average fiber ROI measures
%     {default = 1}
%   'thresh_FA': FA threshold used to calculate fiber ROI measures
%     {default = 0}
%   'thresh_prob': fiber probability threshold
%     {default = 0}
%   'xcg_flag': [0|1] exclude CSF and gray-mattter from fiber ROIs
%     {default = 0}
%   'masksf_flag': [0|1] voxels with multiple fibers excluded
%     {default = 0}
%   'fiber_atlasname': name of the atlas used in AtlasTrack
%     and attached to output analysis mat files (empty = default atlas)
%     {default = []}
%   'fname_fiber_legend': name of fiber legend csv file
%     if empty, use default ($MMPS_PARMS/DTI_Fiber/DTI_Fiber_Legend.csv)
%     {default = []}
%
% Optional Parameters for fiber volume calculations:
%   'fiber_vol_norm_flag': [0|1] normalize fiber volumes by ICV or total brain
%     volume obtained from FreeSurfer aseg stats (expressed as percentage)
%     {default = 0}
%   'fiber_vol_norm_code': aseg ROI code used for normalizing fiber volumes
%     e.g. WholeBrain = 10001, BrainMaskVolume = 20001
%          BrainSegmentationVolume = 20002, ICV = 20003
%     {default = 20003}
%   'resolution': size of DTI voxels (mm) assumed for calculating volume
%     may be set in MMIL_ProjInfo.csv as DTI_resolution
%     {default = [2 2 2]}
%
% Optional Parameters specific to DTI aseg analysis:
%   'aseg_disp_flag': [0|1] whether to weighted averages were calculated
%     based on MD values and dispvec
%     {default = 0}
%   'aseg_dispfact': multiplicative factor applied to aseg dispersion values
%     {default = 4}
%   'erode_flag': [0|1] whether aseg ROIs were "eroded" by
%     smoothing and thresholding to reduce partial voluming
%     {default = 1}
%   'erode_nvoxels': number of voxels eroded
%     {default = 1}
%   'aseg_aparc_flag': [0|1|2] whether cortical parcellation volume ROIs were used
%     0: aseg only
%     1: aparc only
%     2: aparc+aseg
%     {default = 0}
%
% Optional Parameters specific to DTI wmparc analysis:
%   'wmparc_aparc_flag': [0|1] whether to use cortical parcellation ROIs
%     0: wmparc only
%     1: wmparc and aparc
%     {default = 0}
%
% Optional Parameters specific to DTI cortical surface analysis:
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast
%     {default = [-1,1]}
%   'gwnorm_flag': [0|1] when calculating gray-white contrast, normalize by mean
%     {default = 1}
%
% Optional Parameters specific to MRI analysis:
%   'MRI_projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast
%     {default = [-0.2,0.2]}
%
% Other Optional Parameters:
%   'union_flag': [0|1] output results for the union of subjects
%     with MRI or DTI, otherwise only for the intersection
%     {default = 1}
%   'continfo_flag': [0|1] include container-related information
%     e.g. 'Manufacturer','ManufacturersModelName','DeviceSerialNumber',
%       'MagneticFieldStrength','MMPS_version','ProcDate'
%     {default = 0}
%   'subjinfo_flag': [0|1] include subject information
%     e.g. Age, Sex, Site, Group
%     this file must exist:
%       {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%     {default = 0}
%   'outdir': output directory for summary csv files
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'ROI_Summaries'}
%   'outstem': output file stem
%     {default = 'MRI_DTI'}
%   'concat_outfix': string attached to summary file
%     {default = 'all'}
%   'save_mat_flag': [0|1] save compiled ROI data in a mat file
%     {default = 0}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  11/07/12 by Don Hagler
% Last Mod: 01/22/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
results = [];

% check input parameters
parms = check_input(ProjID,varargin);

% check whether anything to do
if ~parms.thick_flag && ~parms.area_flag && ~parms.cortT1w_flag && ...
   ~parms.subvol_flag && ~parms.subT1w_flag && ...
   ~parms.fiber_flag && ~parms.aseg_flag && ~parms.cortsurf_flag
  fprintf('%s: nothing to do\n',mfilename);
  return;
end;

% check whether output files exist
parms.fname_out = sprintf('%s/%s_%s.csv',...
  parms.outdir,parms.outstem,parms.concat_outfix);
if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

% compile MRI analysis results for each subject
results = compile_results(parms);

% create summary csv files for all analysis types
summarize_analysis(parms,results);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'QC_raw',true,[false true],...
    'QC_recon',true,[false true],...
    'QC_DTI',true,[false true],...
    'qcflag',true,[false true],...
... % specify which MRI results to compile
    'thick_flag',true,[false true],...
    'area_flag',true,[false true],...
    'cortT1w_flag',true,[false true],...
    'subvol_flag',true,[false true],...
    'subT1w_flag',true,[false true],...
    'fuzzy_flag',true,[false true],...
    'fuzzy_order',18,[2,4,12,18],...
... % specify which DTI results to compile
    'measlist',{'FA','MD','LD','TD','T2w','T1w'},[],...
    'scalefacts',[],[],...
    'inputlist',[],[],...
    'motion_flag',true,[false true],...
    'fiber_flag',true,[false true],...
    'fiber_vol_flag',false,[false true],...
    'aseg_flag',true,[false true],...
    'wmparc_flag',false,[false true],...
    'cortsurf_flag',true,[false true],...
... % specify DTI data used for tensor calculations and analysis
    'snums_flag',1,[0:3],...
    'snum_index',[],[],...
    'infix',[],[],...
    'auto_infix_flag',true,[false true],...
    'revflag',0,[0:2],...
    'nob0_flag',false,[false true],...
    'min_bval',1,[],...
    'flex_flag',false,[false true],...
    'min_nb0',1,[],...
    'min_ndirs',6,[],...
    'full_fstem_flag',false,[false true],...    
    'DT_analdir','DTanalysis',[],...
    'DT_outfix',[],[],...
... % DTI fiber analysis
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'atlas_flag',2,[0:4],...
    'weighted_avg_flag',true,[false true],...
    'thresh_FA',0,[],...
    'thresh_prob',0,[],...
    'xcg_flag',false,[false true],...
    'masksf_flag',false,[false true],...
    'fiber_atlasname',[],[],...
    'fname_fiber_legend',[],[],...
... % DTI fiber volume
    'fiber_vol_norm_flag',false,[false true],...
    'fiber_vol_norm_code',20003,[1 Inf],...
    'resolution',[2 2 2],[],...
... % DTI aseg analysis
    'aseg_disp_flag',false,[false true],...
    'aseg_dispfact',4,[1e-6,1e6],...
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
    'aseg_aparc_flag',0,[0,1,2],...
... % DTI wmparc parameters
    'wmparc_aparc_flag',false,[false true],...    
... % DTI cortsurf analysis
    'projdist_list',[-1,1],[-5,5],...
    'gwnorm_flag',true,[false true],...
... % MRI analysis
    'MRI_projdist_list',[-0.2,0.2],[-5,5],...
    'MRI_erode_flag',false,[false true],...
    'MRI_erode_nvoxels',1,[1:100],...
... % other
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'union_flag',true,[false true],...
    'outdir','ROI_Summaries',[],...
    'outstem','MRI_DTI',[],...
    'concat_outfix','all',[],...
    'save_mat_flag',false,[false true],...
    'verbose',1,[0:2],...
    'forceflag',false,[false true],...
...
    'analysis_outdir','analysis',[],...
    'subjinfo_tags',{'Age','Sex','Site','Group'},[],...
    'info_tags',{'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
... % hidden from MMIL_Summarize_MRI_Analysis
    'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79],[1,Inf],...
    'extra_roilist',[10001:10003,20001:20003,20009,20010],[1,Inf],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'fname_fscolorlut',[],[],...
... % hidden from MMIL_Summarize_DTI_Analysis
    'DT_measlist',{'FA','MD','LD','TD','b0N','b0'},[],...
    'resample_flag',true,[false true],...
    'regT1flag',1,[0:2],...
    'xcg_suffix','xcg',[],...
    'masksf_suffix','masksf',[],...
    'fname_fscolorlut',[],[],...
...    'aseg_roilist',[1:28,40:60],[],...
    'aparc_roilist',[1001:1034,2001:2034],[],...
    'exclude_roilist',[1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59,77:79],[],...
...
    'motion_tags',{'mean_motion','mean_trans','mean_rot'},[],...
    'T2w_tags',{'T2w_sf','T2w_r'},[],...
...
    'ProjInfo_tags',{'min_bval','flex_flag','xcg_flag','masksf_flag','fseg_flag',...
                     'cortsurf_flag','atlas_flag','revflag','snums_flag',...
                     'fiber_atlasname','resample_flag','regT1flag'},[],...
    'compile_MRI_tags',{'StudyInfo','RootDirs','qcflag','analysis_outdir',...
                        'projdist_list','thick_flag','area_flag',...
                        'cortT1w_flag','subvol_flag','subT1w_flag',...
                        'fuzzy_flag','fuzzy_order','concat_flag',...
                        'continfo_flag','subjinfo_flag','verbose',...
                        'baseflag','aseg_roilist','extra_roilist',...
                        'aparc_roilist','fname_fscolorlut',...
                        'required_containers','QC_raw','QC_recon',...
                        'modality','erode_nvoxels','erode_flag',...
                        'hemilist','fuzzy_fstem'},[],...
    'compile_DTI_tags',{'StudyInfo','RootDirs','measlist','scalefacts',...
                        'inputlist','motion_flag','fiber_flag','fiber_vol_flag',...
                        'aseg_flag','wmparc_flag','cortsurf_flag',...
                        'concat_flag','snums_flag','snum_index','infix',...
                        'auto_infix_flag','revflag','nob0_flag',...
                        'min_bval','flex_flag',...
                        'min_nb0','min_ndirs','full_fstem_flag',...
                        'DT_analdir','DT_outfix',...
                        'fibers','atlas_flag','weighted_avg_flag','thresh_FA',...
                        'thresh_prob','xcg_flag','masksf_flag','fiber_atlasname',...
                        'fname_fiber_legend','fiber_vol_norm_flag',...
                        'fiber_vol_norm_code','resolution',...
                        'aseg_disp_flag','aseg_dispfact','erode_flag',...
                        'erode_nvoxels','aseg_aparc_flag','wmparc_aparc_flag',...
                        'projdist_list','gwnorm_flag','continfo_flag',...
                        'subjinfo_flag','verbose','DT_measlist','qcflag',...
                        'required_containers','QC_raw','QC_DTI','QC_recon',...
                        'resample_flag','regT1flag','xcg_suffix',...
                        'masksf_suffix','fname_fscolorlut',...
                        'aseg_roilist','aparc_roilist','exclude_roilist'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);

  ProjInfo = MMIL_Get_ProjInfo(ProjID);
  
  % for arg names present in both varargin and ProjInfo,
  %   the options values will appear in merged_args
  %   i.e. command-line takes precedence, then ProjInfo, then defaults
  if ~isempty(ProjInfo)
    % get all 'DTI_' parms, strip 'DTI_'
    ProjInfo_args = MMIL_Args(ProjInfo,'DTI');
    % convert from args cell array to parms struct
    ProjInfo_parms = mmil_args2parms(ProjInfo_args,[],0);
    % keep only those in ProjInfo_tags
    ProjInfo_args = mmil_parms2args(ProjInfo_parms,parms.ProjInfo_tags);
    % merge arguments, giving varargin precedence over ProjInfo_args
    merged_args = mmil_merge_args(options,ProjInfo_args);
    % check that parameters fit allowed range, use defaults if not supplied
    parms = mmil_args2parms(merged_args,parms_filter);
  end;

  % add subjinfo_tags to info_tags only if subjinfo_flag = 1
  if parms.subjinfo_flag
    parms.info_tags = cat(2,parms.subjinfo_tags,parms.info_tags);
  end;

  parms.concat_flag = 1;

  % check outdir
  if mmil_isrelative(parms.outdir)
    parms.outdir = [getenv('HOME') '/MetaData/' parms.ProjID '/' parms.outdir];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_results(parms)
  results = [];
  fname_out = sprintf('%s/%s_roi_data.mat',parms.outdir,parms.outstem);
  if ~exist(fname_out,'file') || parms.forceflag || ~parms.save_mat_flag
    mmil_mkdir(parms.outdir);
    MRI_results = compile_MRI_results(parms);
    DTI_results = compile_DTI_results(parms);
    results = merge_results(parms,MRI_results,DTI_results);
    if parms.save_mat_flag, save(fname_out,'results'); end;
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_MRI_results(parms)
  all_args = mmil_parms2args(parms);
  MRI_args = MMIL_Args(parms,'MRI'); % get MRI parms, strip 'MRI_'
  merged_args = mmil_merge_args(MRI_args,all_args); % replace parms
  tmp_parms = mmil_args2parms(merged_args,[],0);
  args = mmil_parms2args(tmp_parms,parms.compile_MRI_tags);
  if parms.verbose==2
    fprintf('%s: compiling MRI results...\n',mfilename);
  end;
  results = MMIL_Compile_MRI_Analysis(parms.ProjID,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_DTI_results(parms)
  args = mmil_parms2args(parms,parms.compile_DTI_tags);
  if parms.verbose==2
    fprintf('%s: compiling DTI results...\n',mfilename);
  end;
  results = MMIL_Compile_DTI_Analysis(parms.ProjID,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = merge_results(parms,MRI_results,DTI_results);
  results = [];
  if parms.verbose==2
    fprintf('%s: merging MRI and DTI results...\n',mfilename);
  end;
  MRI_VisitIDs = MRI_results.VisitIDs;
  DTI_VisitIDs = DTI_results.VisitIDs;
  if parms.union_flag
    VisitIDs = union(MRI_VisitIDs,DTI_VisitIDs);
  else
    VisitIDs = intersect(MRI_VisitIDs,DTI_VisitIDs);
  end;
  % merge StudyInfo
  StudyInfo = [];
  for s=1:length(VisitIDs)
    VisitID = VisitIDs{s};
    ind = find(strcmp(VisitID,DTI_VisitIDs));
    if isempty(ind)
      ind = find(strcmp(VisitID,MRI_VisitIDs));
      tmp_info = MRI_results.StudyInfo(ind);
    else
      tmp_info = DTI_results.StudyInfo(ind);
    end;
    tags = cat(2,{'SubjID','VisitID','STRUCT_VisitID'},parms.info_tags);
    for t=1:length(tags);
      tag = tags{t};
      StudyInfo(s).(tag) = mmil_getfield(tmp_info,tag);
    end;
  end;
  % add fields to results
  results.StudyInfo = StudyInfo;
  results.RootDirs = DTI_results.RootDirs;
  results.SubjIDs = {StudyInfo.SubjID};
  results.VisitIDs = {StudyInfo.VisitID};
  results.STRUCT_VisitIDs = {StudyInfo.STRUCT_VisitID};
  results.nsubs = length(VisitIDs);

  % copy motion measures from DTI_results to results
  if parms.motion_flag
    for i=1:length(parms.motion_tags)
      tag = parms.motion_tags{i};
      [tmp,ind_DTI] = intersect(DTI_VisitIDs,results.VisitIDs);
      results.DTI.(tag) = nan(1,results.nsubs);
      results.DTI.(tag)(ind_DTI) = DTI_results.(tag);
    end;
  end;
  
  % copy T2w_sf and T2w_r from DTI_results to results
  if ismember('T2w',parms.measlist)
    for i=1:length(parms.T2w_tags)
      tag = parms.T2w_tags{i};
      [tmp,ind_DTI] = intersect(DTI_VisitIDs,results.VisitIDs);
      results.DTI.(tag) = nan(1,results.nsubs);
      results.DTI.(tag)(ind_DTI) = DTI_results.(tag);
    end;
  end;

  % init All struct
  results.nrois = MRI_results.all.nrois + DTI_results.all.nrois;
  results.data = nan(results.nsubs,results.nrois);
  results.roicodes = ...
    cat(1,MRI_results.all.roicodes,DTI_results.all.roicodes);
  results.roinames = cell(results.nrois,1);
  r = 1;
  % prepend roinames with 'MRI_' or 'DTI_'
  for k=1:MRI_results.all.nrois
    results.roinames{r} = ['MRI_' MRI_results.all.roinames{k}];
    r = r + 1;
  end;
  for k=1:DTI_results.all.nrois
    results.roinames{r} = ['DTI_' DTI_results.all.roinames{k}];
    r = r + 1;
  end;
  % merge data
  ind_MRI_ROIs = [1:MRI_results.all.nrois];
  ind_DTI_ROIs = [MRI_results.all.nrois+1:results.nrois];
  for s=1:results.nsubs
    VisitID = VisitIDs{s};
    ind_MRI_visit = find(strcmp(VisitID,MRI_VisitIDs));
    if ~isempty(ind_MRI_visit)
      results.data(s,ind_MRI_ROIs) = ...
        MRI_results.all.data(ind_MRI_visit,:);
    end;
    ind_DTI_visit = find(strcmp(VisitID,DTI_VisitIDs));
    if ~isempty(ind_DTI_visit)
      results.data(s,ind_DTI_ROIs) = ...
        DTI_results.all.data(ind_DTI_visit,:);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_analysis(parms,results)
  fname = parms.fname_out;
  if ~exist(fname,'file') || parms.forceflag
    mmil_mkdir(parms.outdir);
    if parms.verbose==2
      fprintf('%s: writing %s...\n',mfilename,fname);
    end;
    % initialize data cell array
    data = cell(results.nsubs,0);
    col_labels = cell(1,0);
    % add SubjIDs as first column    
    data = cat(2,data,results.SubjIDs');
    col_labels = cat(2,col_labels,'SubjID');
    % add VisitIDs as second column
    data = cat(2,data,results.VisitIDs');
    col_labels = cat(2,col_labels,'VisitID');
    % add info from StudyInfo
    if parms.continfo_flag || parms.subjinfo_flag
      for t=1:length(parms.info_tags)
        tag = parms.info_tags{t};
        if isfield(results.StudyInfo,tag)
          info = {results.StudyInfo.(tag)}';
          data = cat(2,data,info);
          col_labels = cat(2,col_labels,tag);
        end;
      end;
    end;
    % add ROI data
    tmp_data = results.data;
    data = cat(2,data,num2cell(tmp_data));
    col_labels = cat(2,col_labels,mmil_rowvec(results.roinames));
    % add motion measures
    if parms.motion_flag
      for i=1:length(parms.motion_tags)
        tag = parms.motion_tags{i};
        data = cat(2,data,num2cell(results.DTI.(tag)'));
        col_labels = cat(2,col_labels,['DTI_' tag]);
      end;
    end;
    % add T2w measures
    if ismember('T2w',parms.measlist)
      for i=1:length(parms.T2w_tags)
        tag = parms.T2w_tags{i};
        data = cat(2,data,num2cell(results.DTI.(tag)'));
        col_labels = cat(2,col_labels,['DTI_' tag]);
      end;
    end;  
    % add column labels
    data = cat(1,col_labels,data);
    % write file
    mmil_write_csv(fname,data);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
