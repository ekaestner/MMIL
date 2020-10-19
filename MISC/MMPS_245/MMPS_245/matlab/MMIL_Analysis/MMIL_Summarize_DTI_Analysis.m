function results = MMIL_Summarize_DTI_Analysis(ProjID,varargin)
% function results = MMIL_Summarize_DTI_Analysis(ProjID,[options])
%
% Purpose: create summary spreadsheets with DTI ROI measures
%
% Usage:
%  results = MMIL_Summarize_DTI_Analysis(ProjID,'key1', value1,...);
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
% Optional Parameters that specify which analysis results to compile:
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
%   'motion_flag': [0|1] compile head motion measures
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
%   'gwcsurf_flag': [0|1] compile cortical gray/white ROI results
%     {default = 0}
%   'concat_flag': [0|1|2] summarize concatenated results for all analysis types
%     if 0, create separate summaries for each analysis type and measure
%     if 1, create a single summary for all analysis types and measures
%     if 2, create both individual and concatenated summaries
%     {default = 1}
%   'concat_outfix': string attached to concatenated summary file
%     {default = 'all'}
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
%  'max_bval': maximum b-value used in tensor fit
%     {default = Inf}
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%  'min_ndirs': minimum number of diffusion directions to be included
%     {default = 6}
%  'full_fstem_flag': [0|1] full stem included in DT measures analysis output
%     otherwise, use only meas name (e.g. 'FA' instead of 'DTI_scans_1_2_...FA')
%     if 0, options above (snums_flag, snum_index, infix, etc.) are ignored
%     {default = 0}
%
% Optional Parameters specific to fiber analysis:
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
%   'fiber_disp_flag': [0|1] whether to calculate weighted averages
%     based on MD values and dispvec
%     {default = 0}
%   'fiber_dispfact': multiplicative factor applied to fiber dispersion values
%     {default = 4}
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
% Optional Parameters specific to aseg analysis:
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
% Optional Parameters specific to wmparc analysis:
%   'wmparc_aparc_flag': [0|1] whether to use cortical parcellation ROIs
%     0: wmparc only
%     1: wmparc and aparc
%     {default = 0}
%
% Optional Parameters specific to cortical surface analysis:
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast
%     {default = [-1,1]}
%
% Optional Parameters specific to gwcsurf analysis:
%   'gwcsurf_outfix': suffix added to output files for gwcsurf analysis
%     {default = 'gwcsurf'}
%
% Optional Parameters for cortsurf and gwcsurf analysis:
%   'gwnorm_flag': [0|1] when calculating gray-white contrast, normalize by mean
%     {default = 1}
%   'area_flag': [0|1] when calculating mean values, weight by cortical area
%     requires that MRI analysis have been run
%     {default = 1}
%
% Other Optional Parameters:
%   'DT_analdir': name of DTI analysis folder relative to ContainerPath
%     {default = 'DTanalysis'}
%   'DT_outfix': string attached to DT calculation and analysis file names
%     {default = []}
%   'continfo_flag': [0|1] include container-related information
%     e.g. 'Manufacturer','ManufacturersModelName','DeviceSerialNumber',
%       'MagneticFieldStrength','MMPS_version','ProcDate'
%     {default = 0}
%   'subjinfo_flag': [0|1] include subject information
%     e.g. Age, Sex, Site, Group
%     this file must exist:
%       {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%     {default = 0}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%   'outdir': output directory for summary csv files
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'ROI_Summaries'}
%   'outstem': output file stem
%     {default = 'DTI'}
%   'full_DTcalc_outfix_flag': [0|1] whether to include strings summarizing
%     the details of DTI calculation in the output summary files
%     {default = 0}
%   'short_DTcalc_outfix': string to summarize
%     the details of DTI calculation in the output summary files
%     used if full_DTcalc_outfix_flag = 0
%     {default = []}
%   'save_mat_flag': [0|1] save compiled ROI data in a mat file
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  02/03/11 by C Roddey
% Prev Mod: 04/14/16 by Don Hagler
% Last Mod: 08/03/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
results = [];

% check input parameters
parms = check_input(ProjID,varargin);

% check whether anything to do
if ~parms.fiber_flag && ~parms.aseg_flag &&...
   ~parms.wmparc_flag && ~parms.cortsurf_flag &&...
   ~parms.gwcsurf_flag
  fprintf('%s: nothing to do\n',mfilename);
  return;
end;

% check whether output files exist
if ~parms.forceflag
  runflag = check_output(parms);
  if ~runflag, return; end;
end;

% compile DTI analysis results for each subject
results = compile_results(parms);

% create summary csv files for individual roitypes and measuress
if parms.concat_flag~=1
  if parms.fiber_flag
    summarize_roitype_analysis(parms,results,'fiber');
  end;
  if parms.aseg_flag
    summarize_roitype_analysis(parms,results,'aseg');
  end;
  if parms.wmparc_flag
    summarize_roitype_analysis(parms,results,'wmparc');
  end;
  if parms.cortsurf_flag
    for p=1:parms.nprojdist
      summarize_roitype_analysis(parms,results,'cortsurf',p);
    end;
    if parms.contrast_flag
      summarize_roitype_analysis(parms,results,'contrast');
    end;
  end;
  if parms.gwcsurf_flag
    for p=1:parms.ngwlayers
      summarize_roitype_analysis(parms,results,'gwcsurf',p);
    end;
  end;
end;

% create summary csv files for all roitypes and measuress
if parms.concat_flag
  summarize_concat_analysis(parms,results);
end;

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
... % specify which results to compile
    'measlist',{'FA','MD','LD','TD','T2w','T1w'},[],...
    'scalefacts',[],[],...
    'inputlist',[],[],...
    'motion_flag',true,[false true],...
    'fiber_flag',true,[false true],...
    'fiber_vol_flag',false,[false true],...
    'aseg_flag',true,[false true],...
    'wmparc_flag',false,[false true],...
    'cortsurf_flag',true,[false true],...
    'gwcsurf_flag',false,[false true],...
    'concat_flag',1,[0:2],...
    'concat_outfix','all',[],...
... % specify DTI data used for tensor calculations and analysis
    'snums_flag',1,[0:3],...
    'snum_index',[],[],...
    'infix',[],[],...
    'auto_infix_flag',true,[false true],...
    'revflag',0,[0:2],...
    'nob0_flag',false,[false true],...
    'min_bval',1,[],...
    'max_bval',Inf,[100,Inf],...
    'flex_flag',false,[false true],...
    'min_nb0',1,[],...
    'min_ndirs',6,[],...
    'full_fstem_flag',false,[false true],...    
... % fiber analysis
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'atlas_flag',2,[0:4],...
    'weighted_avg_flag',true,[false true],...
    'thresh_FA',0,[],...
    'thresh_prob',0,[],...
    'fiber_disp_flag',false,[false true],...
    'fiber_dispfact',4,[1e-6,1e6],...
    'xcg_flag',false,[false true],...
    'masksf_flag',false,[false true],...
    'fiber_atlasname',[],[],...
    'fname_fiber_legend',[],[],...
... % fiber volume
    'fiber_vol_norm_flag',false,[false true],...
    'fiber_vol_norm_code',20003,[1 Inf],...
    'resolution',[2 2 2],[],...
... % aseg analysis
    'aseg_disp_flag',false,[false true],...
    'aseg_dispfact',4,[1e-6,1e6],...
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
    'aseg_aparc_flag',0,[0,1,2],...
... % wmparc parameters
    'wmparc_aparc_flag',false,[false true],...    
... % cortsurf analysis
    'projdist_list',[-1,1],[-5,5],...
... % gwcsurf analysis
    'gwcsurf_outfix','gwcsurf',[],...
... % cortsurf and gwcsurf analysis
    'gwnorm_flag',true,[false true],...
    'area_flag',true,[false true],...
... % other
    'DT_analdir','DTanalysis',[],...
    'DT_outfix',[],[],...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'verbose',1,[0:2],...
    'outdir','ROI_Summaries',[],...
    'outstem','DTI',[],...
    'full_DTcalc_outfix_flag',false,[false true],...
    'short_DTcalc_outfix',[],[],...
    'save_mat_flag',false,[false true],...
    'forceflag',false,[false true],...
... % hidden
    'fnum',1,[1,Inf],...
    'DT_outdir','DTcalc',[],...
    'RSI_analdir','RSIanalysis',[],...
    'RSI_outfix',[],[],...
    'RSI_outdir','RSIcalc',[],...
    'DT_measlist',{'FA','MD','LD','TD','b0','b0N','T2w'},[],...
    'RSI_measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
                    'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'},[],...
    'multiframe_measlist',{'F0','N0','F2','N2','F4','N4',...
                           'FD','ND','FT','NT'},[],...
    'required_containers',{'proc_dti','fsurf'},[],...
    'resample_flag',true,[false true],...
    'regT1flag',1,[0:2],...
    'aseg_disp_suffix','dwtd',[],...
    'fiber_disp_suffix','dwtd',[],...
    'xcg_suffix','xcg',[],...
    'masksf_suffix','masksf',[],...
    'fname_fscolorlut',[],[],...
    'aseg_roilist',[1:28,40:60],[],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'exclude_roilist',[1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59],[],...
    'wmparc_roilist',[3001:3003,3005:3034,4001:4003,4005:4034],[1,Inf],...
...
    'motion_tags',{'mean_motion','mean_trans','mean_rot'},[],...
    'T2w_tags',{'T2w_sf','T2w_r'},[],...
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'ProjInfo_tags',{'min_bval','max_bval','flex_flag',...
                     'xcg_flag','masksf_flag','fseg_flag',...
                     'wmparc_flag','cortsurf_flag','gwcsurf_flag',...
                     'atlas_flag','revflag','snums_flag',...
                     'fiber_atlasname','resample_flag','regT1flag',...
                     'RSI_outdir','RSI_outfix','DT_outdir','DT_outfix'},[],...
    'excl_tags',{'excl_tags','ProjInfo_tags','ProjID',...
                 'fnum','multiframe_measlist',...
                 'concat_outfix','outdir','outstem',...
                 'full_DTcalc_outfix_flag','short_DTcalc_outfix',...
                 'save_mat_flag',...
                 'forceflag','combined_measlist','DTI_flags',...
                 'motion_tags','T2w_tags','info_tags'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  try
    ProjInfo = MMIL_Get_ProjInfo(ProjID);
  catch me
    fprintf('%s: WARNING: %s\n',mfilename,me.message);
    ProjInfo = [];
  end;
  
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
 
  parms.compile_tags = setdiff(fieldnames(parms),parms.excl_tags);

  if mmil_isrelative(parms.outdir)
    parms.outdir = [getenv('HOME') '/MetaData/' parms.ProjID '/' parms.outdir];
  end;

  if parms.cortsurf_flag
    % whether to calculate gray-white contrast depends on if 2 projdist
    parms.nprojdist = length(parms.projdist_list);
    if parms.nprojdist==2
      parms.contrast_flag = 1;
      parms.projdist_list = sort(parms.projdist_list); % white matter first
    else
      parms.contrast_flag = 0;
    end;
  else
    parms.contrast_flag = 0;
  end;

  if parms.gwcsurf_flag
    parms.ngwlayers = 3;
  end;

  % check for multi-frame measures
  parms = check_multiframe(parms);

  parms.nmeas = length(parms.measlist);
  parms = combine_measlist(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_multiframe(parms)
  if (length(parms.fnum)==1 && parms.fnum==1) ||...
     isempty(intersect(parms.measlist,parms.multiframe_measlist))
    return;
  end;
  measlist = parms.measlist;
  j = 1;
  parms.measlist = [];
  for m=1:length(measlist)
    meas = measlist{m};
    if ~ismember(meas,parms.multiframe_measlist)
      fnum = 1;
    else
      fnum = parms.fnum;
    end;
    for f=fnum
      if f==1
        parms.measlist{j} = meas;
      else
        multi_meas = sprintf('%ss%d',meas,f);
        if ismember(meas,parms.RSI_measlist)
          parms.RSI_measlist = union(parms.RSI_measlist,multi_meas);
        elseif ismember(meas,parms.DT_measlist)
          parms.DT_measlist = union(parms.DT_measlist,multi_meas);
        end;
        parms.measlist{j} = multi_meas;
      end;
      j = j + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = combine_measlist(parms)
  parms.combined_measlist = parms.measlist;
  parms.DTI_flags = zeros(1,parms.nmeas);
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    if ismember(meas,parms.DT_measlist) ||...
       ismember(meas,parms.RSI_measlist)
      parms.DTI_flags(m) = 1;
    end;
  end;
  m = parms.nmeas + 1;
  for i=1:length(parms.inputlist)
    [tmp,meas] = fileparts(parms.inputlist{i});
    parms.combined_measlist{m} = meas;
    m = m + 1;
  end;
  parms.nmeas = length(parms.combined_measlist);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runflag = check_output(parms)
  runflag = 0;
  if parms.concat_flag~=1
    if parms.fiber_flag
      all_exist = check_output_files(parms,'fiber');
      if ~all_exist
        runflag = 1;
        return;
      end;
    end;
    if parms.aseg_flag
      all_exist = check_output_files(parms,'aseg');
      if ~all_exist
        runflag = 1;
        return;
      end;
    end;
    if parms.wmparc_flag
      all_exist = check_output_files(parms,'wmparc');
      if ~all_exist
        runflag = 1;
        return;
      end;
    end;
    if parms.cortsurf_flag
      for p=1:parms.nprojdist
        all_exist = check_output_files(parms,'cortsurf',p);
        if ~all_exist
          runflag = 1;
          return;
        end;
      end;
      if parms.contrast_flag
        all_exist = check_output_files(parms,'contrast');
        if ~all_exist
          runflag = 1;
          return;
        end;
      end;
    end;
    if parms.gwcsurf_flag
      for p=1:parms.ngwlayers
        all_exist = check_output_files(parms,'gwcsurf',p);
        if ~all_exist
          runflag = 1;
          return;
        end;
      end;
    end;
  end;
  if parms.concat_flag
    outfix = set_concat_outfix(parms);
    fname = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
    if ~exist(fname,'file')
      runflag = 1;
      return;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_exist = check_output_files(parms,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  all_exist = 1;
  [measlist,nmeas,DTI_flags] = set_measlist(parms,roitype);
  for m=1:nmeas
    meas = measlist{m};
    DTI_flag = DTI_flags(m);
    switch roitype
      case 'fiber'
        outfix = set_fiber_outfix(parms,DTI_flag);
      case 'aseg'
        outfix = set_aseg_outfix(parms,DTI_flag);
      case 'wmparc'
        outfix = set_wmparc_outfix(parmx,DTI_flag);
      case 'cortsurf'  
        outfix = set_cortsurf_outfix(parms,DTI_flag,p);
      case 'contrast'  
        outfix = set_cortsurf_outfix(parms,DTI_flag,0);
      case 'gwcsurf'  
        outfix = set_gwcsurf_outfix(parms,DTI_flag,p);
    end;
    fname = sprintf('%s/%s_%s_%s.csv',parms.outdir,parms.outstem,meas,outfix);
    if ~exist(fname,'file')
      all_exist = 0;
      return;
    end;
  end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [measlist,nmeas,DTI_flags] = set_measlist(parms,roitype)
  measlist = parms.combined_measlist;
  DTI_flags = parms.DTI_flags;
  if strcmp(roitype,'fiber') && parms.fiber_vol_flag
    measlist = cat(2,measlist,'vol');
    DTI_flags = cat(2,DTI_flags,1);
  end;
  nmeas = length(measlist);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_results(parms)
  fname_out = sprintf('%s/%s_roi_data.mat',parms.outdir,parms.outstem);
  if ~exist(fname_out,'file') || parms.forceflag || ~parms.save_mat_flag
    mmil_mkdir(parms.outdir);
    parms.concat_flag = (parms.concat_flag>0);
    args = mmil_parms2args(parms,parms.compile_tags);
    results = MMIL_Compile_DTI_Analysis(parms.ProjID,args{:});
    if parms.save_mat_flag, save(fname_out,'results'); end;
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_csv(fname,results,parms,roitype,meas,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  if ~exist('meas','var'), meas = []; end;
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
  % add motion measures
  if parms.motion_flag
    for i=1:length(parms.motion_tags)
      tag = parms.motion_tags{i};
      data = cat(2,data,num2cell(results.(tag)'));
      col_labels = cat(2,col_labels,tag);
    end;
  end;
  % add T2w measures
  if ismember('T2w',parms.measlist)
    for i=1:length(parms.T2w_tags)
      tag = parms.T2w_tags{i};
      data = cat(2,data,num2cell(results.(tag)'));
      col_labels = cat(2,col_labels,tag);
    end;
  end;
  % add ROI data
  if strcmp(roitype,'all')
    tmp_data = results.(roitype).data;
  else
    tmp_data = results.(roitype).(meas).data(:,:,p);
  end;
  data = cat(2,data,num2cell(tmp_data));
  col_labels = cat(2,col_labels,mmil_rowvec(results.(roitype).roinames));
  % add column labels
  data = cat(1,col_labels,data);
  % write file
  mmil_write_csv(fname,data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_roitype_analysis(parms,results,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  [measlist,nmeas,DTI_flags] = set_measlist(parms,roitype);
  for m=1:nmeas
    meas = measlist{m};
    DTI_flag = DTI_flags(m);
    switch roitype
      case 'fiber'
        outfix = set_fiber_outfix(parms,DTI_flag);
      case 'aseg'
        outfix = set_aseg_outfix(parms,DTI_flag);
      case 'wmparc'
        outfix = set_wmparc_outfix(parms,DTI_flag);
      case 'cortsurf'  
        outfix = set_cortsurf_outfix(parms,DTI_flag,p);
      case 'contrast'  
        outfix = set_cortsurf_outfix(parms,DTI_flag,0);
      case 'gwcsurf'
        outfix = set_gwcsurf_outfix(parms,DTI_flag,p);
    end;
    fname = sprintf('%s/%s_%s_%s.csv',parms.outdir,parms.outstem,meas,outfix);
    if ~exist(fname,'file') || parms.forceflag
      mmil_mkdir(parms.outdir);
      if parms.verbose==2
        fprintf('%s: writing %s...\n',mfilename,fname);
      end;
      write_csv(fname,results,parms,roitype,meas,p);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_concat_analysis(parms,results)
  outfix = set_concat_outfix(parms);
  fname = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,outfix);
  if ~exist(fname,'file') || parms.forceflag
    mmil_mkdir(parms.outdir);
    if parms.verbose==2
      fprintf('%s: writing %s...\n',mfilename,fname);
    end;
    write_csv(fname,results,parms,'all');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_fiber_outfix(parms,DTI_flag)
  outfix = 'fibers';
  switch parms.atlas_flag
    case 0 % manual
    case 1 % loc only, count atlas
      outfix = [outfix '_loc_countatlas'];
    case 2 % loc+dir, count atlas
      outfix = [outfix '_countatlas'];
    case 3 % loc only, mask atlas
      outfix = [outfix '_loc_atlas'];
    case 4 % loc+dir, mask atlas
      outfix = [outfix '_atlas'];
  end;
  if parms.weighted_avg_flag
    outfix = [outfix '_wtd'];
  end;
  if parms.fiber_disp_flag
    outfix = sprintf('%s_%s%0.1f',...
      outfix,parms.fiber_disp_suffix,parms.fiber_dispfact);
  end;
  if parms.xcg_flag
    outfix = [outfix '_' parms.xcg_suffix];
  end;
  if parms.masksf_flag
    outfix = [outfix '_' parms.masksf_suffix];
  end;
  if parms.thresh_prob>0 && parms.atlas_flag>0
    outfix = sprintf('%s_pthresh%0.2f',outfix,parms.thresh_prob);
  end;
  if parms.thresh_FA>0
    outfix = sprintf('%s_FAthresh%0.2f',outfix,parms.thresh_FA);
  end;
  if ~isempty(parms.fiber_atlasname)
    outfix = [outfix '_' parms.fiber_atlasname];
  end
  if DTI_flag
    outfix = add_DTI_outfix(outfix,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_aseg_outfix(parms,DTI_flag)
  switch parms.aseg_aparc_flag
    case 0
      outfix = 'aseg';
    case 1
      outfix = 'aparc';
    case 2
      outfix = 'aparc+aseg';
  end;
  if parms.erode_flag
    outfix = [outfix '_erode'];
  end
  if parms.aseg_disp_flag
    outfix = sprintf('%s_%s%0.1f',...
      outfix,parms.aseg_disp_suffix,parms.aseg_dispfact);
  end;
  if DTI_flag
    outfix = add_DTI_outfix(outfix,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_wmparc_outfix(parms,DTI_flag)
  outfix = 'wmparc';
  if DTI_flag
    outfix = add_DTI_outfix(outfix,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_cortsurf_outfix(parms,DTI_flag,p)
  outfix = 'cortsurf';
  if DTI_flag
    outfix = add_DTI_outfix(outfix,parms);
  end;
  if p==0
    outfix = [outfix '_contrast'];
  else
    projdist = parms.projdist_list(p);
    if projdist<0
      outfix = sprintf('%s_white%0.1f',outfix,projdist);
    else
      outfix = sprintf('%s_gray+%0.1f',outfix,projdist);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_gwcsurf_outfix(parms,DTI_flag,p)
  outfix = 'gwcsurf';
  if DTI_flag
    outfix = add_DTI_outfix(outfix,parms);
  end;
  switch p
    case 1
      outfix = [outfix '_wm'];
    case 2
      outfix = [outfix '_gm'];
    case 3
      outfix = [outfix '_gwc'];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = add_DTI_outfix(outfix,parms)
  if parms.full_DTcalc_outfix_flag
    if parms.snums_flag
      outfix = sprintf('%s_snums%d',outfix,parms.snums_flag);
    end;
    if ~isempty(parms.snum_index)
      outfix = sprintf('%s_snumidx%d',outfix,parms.snum_index);
    end;
    if parms.nob0_flag
      outfix = [outfix '_nob0'];
    end;
    if parms.min_bval ~= 1
      outfix = sprintf('%s_minb%d',outfix,parms.min_bval);
    end
    if parms.flex_flag
      outfix = [outfix '_flex'];
    end;
    if parms.min_ndirs ~= 6
      outfix = sprintf('%s_mind%d',outfix,parms.min_ndirs);
    end
  elseif ~isempty(parms.short_DTcalc_outfix)
    outfix = [outfix '_' parms.short_DTcalc_outfix];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfix = set_concat_outfix(parms)
  outfix = parms.concat_outfix;
return;

