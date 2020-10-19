function results = MMIL_Compile_DTI_Analysis(ProjID,varargin)
%function results = MMIL_Compile_DTI_Analysis(ProjID,[options])
%
% Usage:
%  results = MMIL_Compile_DTI_Analysis(ProjID,'key1', value1,...);
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
%   'QC_raw' : [0|1] use good raw QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_RawQC.mat
%     {default = 1}
%   'QC_recon' : [0|1] use recon QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_FSReconQC.csv
%     {default = 1}
%   'QC_DTI' : [0|1] use DTIreg QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_DTIQC.csv
%     {default = 1}
%   'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
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
%     {default = 0}
%   'wmparc_flag': [0|1] extract wmparc ROI results
%     {default = 0}
%   'cortsurf_flag': [0|1] compile cortical surface ROI results
%     {default = 0}
%   'gwcsurf_flag': [0|1] compile cortical gray/white ROI results
%     {default = 0}
%   'area_flag': [0|1] when calculating mean values, weight by cortical area
%     requires that MRI analysis have been run
%     {default = 1}
%   'concat_flag': [0|1] concatenate results for all analysis types
%     {default = 0}
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
%  'DT_analdir': name of DTI analysis folder relative to ContainerPath
%     {default = 'DTanalysis'}
%  'DT_outfix': string attached to DT calculation and analysis file names
%     {default = []}
%
% Optional Parameters specific to how fiber analysis was done:
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
% Optional Parameters specific to how aseg analysis was done:
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
% Optional Parameters specific to how cortical surface analysis was done:
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
%
% Other Optional Parameters:
%   'continfo_flag': [0|1] include container-related information in StudyInfo
%     e.g. 'Manufacturer','ManufacturersModelName','DeviceSerialNumber',
%       'MagneticFieldStrength','MMPS_version','ProcDate'
%     {default = 0}
%   'subjinfo_flag': [0|1] include subject information in StudyInfo
%     e.g. Sex, Site, DOB, Group
%     this file must exist:
%       {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%     {default = 0}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%
% Output:
%   results: struct containing compiled DTI analysis results with fields:
%     StudyInfo: struct array containing info for each subject
%     RootDirs: struct containing paths of root directories
%     SubjIDs: cell array of subject IDs
%     VisitIDs: cell array of visit IDs
%     STRUCT_VisitIDs: cell array of structural (i.e. FreeSurfer) visit IDs
%     VisitNumbers: vector of visit numbers
%     nsubs: number of subjects/visits in StudyInfo
%     nmeas: number of measures in measlist
%     fiber: struct containing compiled data for fiber ROIs (if fiber_flag = 1)
%       nrois: number of fiber ROIs
%       roicodes: ROI codes
%       roinames: ROI names
%       FA: struct containing data for FA (other measures have similar fields)
%         data: matrix with size [nsubs,nrois] containing ROI averages
%         std: matrix of standard deviations
%         nvals: matrix of number of voxels for each subject and ROI
%         nvals_valid: number of values included in ROI averages
%         nvals_invalid: number of values excluded from ROI averages
%     aseg       (if aseg_flag = 1)
%     wmparc     (if wmparc_flag = 1)
%     cortsurf   (if cortsurf_flag = 1)
%     all   (if concat_flag = 1)
%       nrois: number of ROIs of all types
%       roicodes: ROI codes
%       roinames: ROI names with meas prepended
%       data: matrix of ROI averages for all measures and ROI types
%       std: matrix of standard deviations
%       nvals: matrix of number of voxels for each subject and ROI
%       nvals_valid: number of values included in ROI averages
%       nvals_invalid: number of values excluded from ROI averages
%
% Created:  10/06/12 by Don Hagler
% Prev Mod: 07/31/15 by Don Hagler
% Last Mod: 08/03/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);

results = init_results(parms);

% compile motion measures
if parms.motion_flag
  results = compile_motion(parms,results);
end;

% compile T2w scaling factors
if ismember('T2w',parms.measlist)
  results = compile_T2w_sf(parms,results);
end;

% compile fiber results
if parms.fiber_flag
  results = compile_roi_data(parms,results,'fiber');
  if parms.fiber_vol_flag
    results = calc_fiber_vol(parms,results);
  end;
end;

% compile aseg results
if parms.aseg_flag
  results = compile_roi_data(parms,results,'aseg');
end;

% compile wmparc results
if parms.wmparc_flag
  results = compile_roi_data(parms,results,'wmparc');
end;

% compile cortsurf results
if parms.cortsurf_flag
  results = compile_roi_data(parms,results,'cortsurf');
  % calculate difference between gray and white matter
  if parms.contrast_flag
    results = calc_contrast(parms,results);
  end;
  % compile cortical area for calculating weighted mean values
  if parms.area_flag
    results = compile_cort_area(parms,results);
  end;
  % calculate mean cortical values
  results = calc_mean_cort_vals(parms,results,'cortsurf');
  if parms.contrast_flag
    results = calc_mean_cort_vals(parms,results,'contrast');
  end;
end;

% compile gwcsurf results
if parms.gwcsurf_flag
  results = compile_roi_data(parms,results,'gwcsurf');
  % compile cortical area for calculating weighted mean values
  if parms.area_flag && ~parms.cortsurf_flag
    results = compile_cort_area(parms,results);
  end;
  % calculate mean cortical values
  results = calc_mean_cort_vals(parms,results,'gwcsurf');
end;

% concatenate results
if parms.concat_flag
  results = concat_results(parms,results);
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
    'aseg_flag',false,[false true],...
    'wmparc_flag',false,[false true],...
    'cortsurf_flag',false,[false true],...
    'gwcsurf_flag',false,[false true],...
    'area_flag',true,[false true],...
    'concat_flag',false,[false true],...
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
    'DT_analdir','DTanalysis',[],...
    'DT_outfix',[],[],...
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
... % other
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'verbose',1,[0:2],...
... % hidden
    'DT_outdir','DTcalc',[],...
    'RSI_outdir','RSIcalc',[],...
    'RSI_analdir','RSIanalysis',[],...
    'RSI_outfix',[],[],...
    'DT_measlist',{'FA','MD','LD','TD','b0','b0N','T2w'},[],...
    'RSI_measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
                    'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'},[],...
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
    'ProjInfo_tags',{'min_bval','max_bval','flex_flag',...
                     'xcg_flag','masksf_flag','fseg_flag',...
                     'cortsurf_flag','gwcsurf_flag','atlas_flag',...
                     'revflag','snums_flag',...
                     'fiber_atlasname','resample_flag','regT1flag',...
                     'resolution'},[],...
    'info_tags',{'snums','revflag','min_nb0','min_ndirs',...
                 'min_bval','flex_flag'},[],...
    'fstem_tags',{'snums','infix','revflag','min_bval','max_bval','flex_flag',...
                  'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
    'aseg_code_tags',{'aseg_aparc_flag','aseg_roilist','aparc_roilist',...
                      'exclude_roilist'},[],...
    'wmparc_code_tags',{'wmparc_aparc_flag','wmparc_roilist','aparc_roilist',...
                        'exclude_roilist'},[],...
    'suffix_tags',{'xcg_flag','xcg_suffix',...
                   'masksf_flag','masksf_suffix',...
                   'disp_flag','disp_suffix','dispfact',...
                   'thresh_prob','thresh_FA'},[],...
    'compile_area_tags',{'StudyInfo','RootDirs','QC_raw','QC_recon',...
      'qcflag','thick_flag','area_flag','cortvol_flag','cortT1w_flag',...
      'subvol_flag','subT1w_flag','aparc_flag','fuzzy_flag','fuzzy_dir',...
      'fuzzy_fstem','fuzzy_order','concat_flag','continfo_flag',...
      'subjinfo_flag','analysis_outdir','check_complete_flag',... % ,'projdist_list'
      'verbose','baseflag','aseg_roilist','extra_roilist',...
      'aparc_roilist','fname_fscolorlut','required_containers','modality',...
      'erode_nvoxels','erode_flag','hemilist','fuzzy_name_tags'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(StudyInfo), error('empty StudyInfo'); end;

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

  parms.StudyInfo = StudyInfo;
  parms.RootDirs = RootDirs;
  
  % set infix
  if parms.auto_infix_flag && isempty(parms.infix)
    parms.infix = 'corr';
    if parms.resample_flag
      if parms.regT1flag==2
        parms.infix = [parms.infix '_regT1'];
      else
        parms.infix = [parms.infix '_resDTI'];
      end;
    end;
  end;

  % set DTI-specific fields in StudyInfo
  [parms.StudyInfo,parms.SIflags] = DTI_MMIL_Check_StudyInfo(...
    parms.StudyInfo,parms.RootDirs,...
    'snums_flag',parms.snums_flag,... % DTIScanNums field will get changed
    'qcflag',parms.qcflag,...
    'DTI_flag',1);
  if isempty(parms.StudyInfo), error('empty StudyInfo after check for DTI'); end;
  parms.nsubs = length(parms.StudyInfo);

  % get container-related info
  if parms.continfo_flag
    parms.StudyInfo = MMIL_Get_ContInfo(parms.StudyInfo,parms.RootDirs);
  end;

  % get subject-related info
  if parms.subjinfo_flag
    parms.StudyInfo = ...
      MMIL_Get_SubjInfo(parms.StudyInfo,parms.RootDirs,parms.ProjID);
  end;

  if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;
  parms.nmeas = length(parms.measlist);
  if isempty(parms.scalefacts)
    parms.scalefacts = ones(1,parms.nmeas);
  elseif numel(parms.scalefacts)~=parms.nmeas
    error('number of scalefacts (%d) does not match measlist (%d)',...
      numel(parms.scalefacts),parms.nmeas);
  end;
  parms = combine_measlist(parms);

  % set parameters for fiber analysis
  if parms.fiber_flag
    if parms.fiber_vol_flag
      if length(parms.resolution)==1
        parms.resolution = parms.resolution*ones(1,3);
      elseif length(parms.resolution)~=3
        error('resolution vector must have one or three elements');
      end;
      parms.voxel_volume = prod(double(parms.resolution));
    end;    
    FiberLegend = DTI_MMIL_Read_Fiber_Legend(parms.fname_fiber_legend);
    legend_fibernums = cell2mat({FiberLegend.FiberNumber});
    legend_fibernames = {FiberLegend.FiberName};    
    [tmp,ind_fibers,ind_legend] = intersect(parms.fibers,legend_fibernums);
    parms.fibers = parms.fibers(ind_fibers);
    parms.fiber_roicodes = mmil_colvec(parms.fibers);
    parms.fiber_roinames = mmil_colvec(legend_fibernames(ind_legend));
    parms.fiber_nrois = length(parms.fibers);
    if parms.fiber_nrois==0
      error('no valid fiber numbers');
    end;      
  end;

  if parms.aseg_flag || parms.cortsurf_flag || parms.gwcsurf_flag
    % read roi codes and names from FreeSurfer color lookup table
    if isempty(parms.fname_fscolorlut)
      MMPS_parms = getenv('MMPS_PARMS');
      parms.fname_fscolorlut = [MMPS_parms '/MMIL_FSColorLUT.txt'];
    end;
    if ~exist(parms.fname_fscolorlut,'file')
      error('FreeSurfer color lookup table file %s not found',...
        parms.fname_fscolorlut);
    end;
    [fs_roicodes,fs_roinames] = fs_colorlut(parms.fname_fscolorlut);
  end;

  % set parameters for aseg analysis
  if parms.aseg_flag
    args = mmil_parms2args(parms,parms.aseg_code_tags);
    parms.aseg_roicodes = mmil_aseg_roicodes(args{:});
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aseg_roicodes);
    parms.aseg_roicodes = fs_roicodes(ind_fs);
    parms.aseg_roinames = fs_roinames(ind_fs);
    parms.aseg_nrois = length(parms.aseg_roicodes);
  end;

  % set parameters for wmparc analysis
  if parms.wmparc_flag
    args = mmil_parms2args(parms,parms.wmparc_code_tags);
    parms.wmparc_roicodes = mmil_wmparc_roicodes(args{:});
    [tmp,ind_fs]=intersect(fs_roicodes,parms.wmparc_roicodes);
    parms.wmparc_roicodes = fs_roicodes(ind_fs);
    parms.wmparc_roinames = fs_roinames(ind_fs);
    parms.wmparc_nrois = length(parms.wmparc_roicodes);
  end;

  % set parameters for cortsurf analysis
  if parms.cortsurf_flag
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aparc_roilist);
    parms.cortsurf_roicodes = fs_roicodes(ind_fs);
    parms.cortsurf_roinames = fs_roinames(ind_fs);
    parms.cortsurf_nrois = length(parms.cortsurf_roicodes);
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

  % set parameters for gwcsurf analysis
  if parms.gwcsurf_flag
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aparc_roilist);
    parms.gwcsurf_roicodes = fs_roicodes(ind_fs);
    parms.gwcsurf_roinames = fs_roinames(ind_fs);
    parms.gwcsurf_nrois = length(parms.gwcsurf_roicodes);
    parms.ngwlayers = 3;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = combine_measlist(parms)
  for m=1:parms.nmeas
    meas = parms.measlist{m};
  end;
  m = parms.nmeas + 1;
  for i=1:length(parms.inputlist)
    [tmp,meas] = fileparts(parms.inputlist{i});
    parms.measlist{m} = meas;
    m = m + 1;
  end;
  parms.nmeas = length(parms.measlist);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = set_fstem(parms,s,meas)
  % set input file stem
  fstem = [];
  ContainerPath = sprintf('%s/%s',...
    parms.RootDirs.proc_dti,parms.StudyInfo(s).proc_dti);
  parms.snums = parms.StudyInfo(s).DTIScanNums;
  args = mmil_parms2args(parms,parms.info_tags);
  [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
  if errcode || isempty(ScanInfo)
    fprintf('%s: WARNING: no DTI scans for %s\n',mfilename,ContainerPath);
    return;
  end;
  switch meas
    case parms.DT_measlist
      if parms.full_fstem_flag
        parms.outdir = parms.DT_outdir;
        parms.outfix = parms.DT_outfix;
        args = mmil_parms2args(parms,parms.fstem_tags);
        fstem = DTI_MMIL_Set_DT_fstem(ContainerPath,args{:});
        if isempty(fstem), return; end;
        [tmp_path,fstem,tmp_ext] = fileparts(fstem);
      else
        fstem = parms.DT_outfix;
        if parms.flex_flag
          if isempty(fstem)
            fstem = 'flex';
          else
            fstem = ['flex_' fstem];
          end;
        end;
      end;
      indir = parms.DT_analdir;
    case parms.RSI_measlist
      if parms.full_fstem_flag
        parms.outdir = parms.RSI_outdir;
        parms.outfix = parms.RSI_outfix;
        args = mmil_parms2args(parms,parms.fstem_tags);
        fstem = DTI_MMIL_Set_RSI_fstem(ContainerPath,args{:});
        if isempty(fstem), return; end;
        [tmp_path,fstem,tmp_ext] = fileparts(fstem);
      else
        fstem = parms.RSI_outfix;
        if parms.flex_flag
          if isempty(fstem)
            fstem = 'flex';
          else
            fstem = ['flex_' fstem];
          end;
        end;
      end;
      indir = parms.RSI_analdir;
    otherwise
      fstem = [];
      indir = parms.DT_analdir;
  end;
  if ~isempty(fstem), fstem = [fstem '_']; end;
  fstem = sprintf('%s/%s/%s',ContainerPath,indir,fstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function infix = set_infix(parms,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  switch roitype
    case 'fiber'
      % set file name stems (must match those set by DTI_MMIL_Analyze_Fibers_Exam)
      infix = 'fiber';
      switch parms.atlas_flag
        case 1 % loc only, count atlas
          infix = [infix '_loc_countatlas'];
        case 2 % loc+dir, count atlas
          infix = [infix '_countatlas'];
        case 3 % loc only, mask atlas
          infix = [infix '_loc_atlas'];
        case 4 % loc+dir, mask atlas
          infix = [infix '_atlas'];
      end;
      if parms.weighted_avg_flag
        infix = [infix '_wtd'];
      end;
      parms.disp_flag = parms.fiber_disp_flag;
      parms.dispfact = parms.fiber_dispfact;
      parms.disp_suffix = parms.fiber_disp_suffix;
      if parms.atlas_flag==0
        parms.thresh_prob = 0;
      end;
      args = mmil_parms2args(parms,parms.suffix_tags);
      suffix = dti_set_fiber_suffix(args{:});
      if ~isempty(suffix)
        infix = [infix '_' suffix];
      end;
      if ~isempty(parms.fiber_atlasname)
        infix = [infix '_' parms.fiber_atlasname];
      end
      infix = [infix '_roi_data'];
    case 'aseg'
      switch parms.aseg_aparc_flag
        case 0
          infix = 'aseg';
        case 1
          infix = 'aparc';
        case 2
          infix = 'aparc+aseg';
      end;
      if parms.erode_flag
        infix = sprintf('%s_erode',infix);
        if parms.erode_nvoxels>1,
          infix = sprintf('%s%d',infix,parms.erode_nvoxels);
        end;
      end;
      if parms.aseg_disp_flag
        infix = sprintf('%s_%s%0.1f',...
          infix,parms.aseg_disp_suffix,parms.aseg_dispfact);
      end;
      infix = [infix '_roi_data'];
    case 'wmparc'
      infix = 'wmparc';
      if parms.erode_flag
        infix = sprintf('%s_erode',infix);
        if parms.erode_nvoxels>1,
          infix = sprintf('%s%d',infix,parms.erode_nvoxels);
        end;
      end;
      infix = [infix '_roi_data'];
    case 'cortsurf'
      infix = sprintf('pdist%0.1f_roi_data',...
        parms.projdist_list(p));
    case 'gwcsurf'
      infix = parms.gwcsurf_outfix;
      switch p
        case 1
          infix = [infix '_wm'];
        case 2
          infix = [infix '_gm'];
        case 3
          infix = [infix '_gwc'];
      end;
      infix = [infix '_roi_data'];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results = [];
  results.StudyInfo = parms.StudyInfo;
  results.RootDirs = parms.RootDirs;
  results.SubjIDs = {parms.StudyInfo.SubjID};
  results.VisitIDs = {parms.StudyInfo.VisitID};
  results.STRUCT_VisitIDs = {parms.StudyInfo.STRUCT_VisitID};
  if ~isfield(parms.StudyInfo,'VisitNumber')
    results.VisitNumbers = zeros(parms.nsubs,1);
  else
    results.VisitNumbers = cell2mat({parms.StudyInfo.VisitNumber});
  end;

  % store relevant numbers
  results.nsubs = parms.nsubs;
  results.nmeas = parms.nmeas;

  % initialize motion measures
  if parms.motion_flag
    results.mean_motion = nan(1,results.nsubs);
    results.mean_trans = nan(1,results.nsubs);
    results.mean_rot = nan(1,results.nsubs);
  end;

  % initialize T2w_sf and T2w_r
  if ismember('T2w',parms.measlist)
    results.T2w_sf = nan(1,results.nsubs);
    results.T2w_r = nan(1,results.nsubs);  
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_roi_results(parms,results,roitype)
  results.(roitype) = [];
  nrois = parms.([roitype '_nrois']);
  switch roitype
    case 'cortsurf'
      np = parms.nprojdist;
    case 'gwcsurf'
      np = parms.ngwlayers;
    otherwise
      np = 1;
  end;
  results.(roitype).nrois = nrois;
  if strcmp(roitype,'fiber')
    results.(roitype).roicodes = parms.([roitype '_roicodes']) + 1e4;
    % NOTE: 1e4 added to roicodes to ensure distinction
    %       between fiber and aseg/cortsurf ROIs
  else
    results.(roitype).roicodes = parms.([roitype '_roicodes']);  
  end;
  results.(roitype).roinames = parms.([roitype '_roinames']);
  init_vals = nan(parms.nsubs,nrois,np);
  for m=1:parms.nmeas
    results.(roitype).(parms.measlist{m}).data = init_vals;
    results.(roitype).(parms.measlist{m}).std = init_vals;
    results.(roitype).(parms.measlist{m}).nvals = init_vals;
    results.(roitype).(parms.measlist{m}).nvals_valid = init_vals;
    results.(roitype).(parms.measlist{m}).nvals_invalid = init_vals;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_all_results(parms,results)
  results.all = [];
  results.all.nrois = 0;
  results.all.roicodes = [];
  results.all.roinames = {};
  results.all.data  = [];
  results.all.std   = [];
  results.all.nvals = [];
  results.all.nvals_valid = [];
  results.all.nvals_invalid = [];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_motion(parms,results)
  for s=1:parms.nsubs
    VisitID = parms.StudyInfo(s).VisitID;
    ContainerPath = sprintf('%s/%s',...
      parms.RootDirs.proc_dti,parms.StudyInfo(s).proc_dti);
    parms.snums = parms.StudyInfo(s).DTIScanNums;
    args = mmil_parms2args(parms,parms.info_tags);
    [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
    if errcode || isempty(ScanInfo)
      fprintf('%s: WARNING: no DTI scans for %s\n',mfilename,ContainerPath);
      return;
    end;
    mean_motion = 0; mean_trans = 0; mean_rot = 0;
    nscans = length(SessInfo.snums_DT);
    n = 0;
    for i=1:nscans 
      j = SessInfo.snums_DT(i);
      fname_qmat = sprintf('%s/DTI%d',ContainerPath,j);
      if SessInfo.revflag
        fname_qmat = [fname_qmat '_rev'];
      end;
      if ~isempty(parms.infix)
        fname_qmat = [fname_qmat '_' parms.infix];
      end;
      fname_qmat = [fname_qmat '_qmat.mat'];
      if ~exist(fname_qmat,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname_qmat);
        continue;
      end;
      n = n + 1;
      tmp = load(fname_qmat);
      mean_motion = mean_motion + tmp.mean_motion;
      mean_trans = mean_trans + tmp.mean_trans;
      mean_rot = mean_rot + tmp.mean_rot;
    end;
    if n>=1
      % average across scans
      results.mean_motion(s) = mean_motion / n;
      results.mean_trans(s) = mean_trans / n;
      results.mean_rot(s) = mean_rot / n;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_T2w_sf(parms,results)
  ind_T2w = find(strcmp('T2w',parms.measlist));
  for s=1:parms.nsubs
    fstem = set_fstem(parms,s,'T2w');
    if isempty(fstem)
      continue;
    end;
    fname_sf = [fstem 'T2w_sf.mat'];
    if ~exist(fname_sf,'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_sf);
      continue;
    end;
    tmp = load(fname_sf);        
    results.T2w_sf(s) = tmp.sf;
    results.T2w_r(s) = tmp.r;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_roi_data(parms,results,roitype)
  results = init_roi_results(parms,results,roitype);
  switch roitype
    case 'cortsurf'
      pvec = [1:parms.nprojdist];
    case 'gwcsurf'
      pvec = [1:parms.ngwlayers];
    otherwise
      pvec = 1;
  end;
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    sf = parms.scalefacts(m);
    if parms.verbose==2
      fprintf('%s: compiling %s %s results for %d subjects...\n',...
        mfilename,meas,roitype,parms.nsubs);
    end;
    for s=1:parms.nsubs
      for p=pvec
        roi_data = load_roi_data(parms,s,meas,roitype,p);
        if isempty(roi_data), continue; end;
        [tmp,ind_data,ind_sel] = ...
          intersect([roi_data.roicode],parms.([roitype '_roicodes']));
        if isempty(ind_data)
          if parms.verbose
            fprintf('%s: WARNING: no valid %s ROIs for %s\n',...
              mfilename,roitype,parms.StudyInfo(s).VisitID);
          end;
          continue;
        end;
        results.(roitype).(meas).data(s,ind_sel,p) = ...
          [roi_data(ind_data).avg]*sf;
        results.(roitype).(meas).std(s,ind_sel,p) = ...
          [roi_data(ind_data).stdv]*sf;
        results.(roitype).(meas).nvals(s,ind_sel,p)  = ...
          [roi_data(ind_data).nvals];
        results.(roitype).(meas).nvals_valid(s,ind_sel,p) = ...
          [roi_data(ind_data).nvals_valid];
        results.(roitype).(meas).nvals_invalid(s,ind_sel,p) = ...
          [roi_data(ind_data).nvals_invalid];
      end;
    end
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_fiber_vol(parms,results)
  % calculate fiber ROI volumes
  meas = parms.measlist{1};
  results.fiber.vol = results.fiber.(meas);
  results.fiber.vol.data = results.fiber.(meas).nvals * parms.voxel_volume;
  results.fiber.vol.std = zeros(results.nsubs,results.fiber.nrois);
  % normalize by ICV or total brain volume
  if parms.fiber_vol_norm_flag  
    % load aseg_stats for each subject
    for s=1:parms.nsubs
      VisitID = parms.StudyInfo(s).STRUCT_VisitID;
      subj = parms.StudyInfo(s).fsurf;
      try
        aseg_stats = ...
          fs_read_aseg_stats(subj,parms.RootDirs.fsurf);
      catch
        fprintf('%s: WARNING: failed to read aseg stats for %s\n',...
          mfilename,VisitID);
        normval = nan;
      end;    
      roicodes = cell2mat({aseg_stats.roicode});
      i_norm = find(roicodes==parms.fiber_vol_norm_code);
      if isempty(i_norm)
        fprintf('%s: WARNING: missing fiber_vol_norm_code %d for %s\n',...
          mfilename,VisitID);
        norm_val = nan;
      else
        norm_val = aseg_stats(i_norm).volume;
      end;
      results.fiber.vol.data(s,:) = 100*results.fiber.vol.data(s,:)/norm_val;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_contrast(parms,results)
  % calculate contrast ratio
  results.contrast = [];
  results.contrast.nrois = results.cortsurf.nrois;
  results.contrast.roicodes = results.cortsurf.roicodes;
  results.contrast.roinames = results.cortsurf.roinames;
  for m=1:parms.nmeas
    if parms.gwnorm_flag
      % normalize by mean
      mean_vals = mean(results.cortsurf.(parms.measlist{m}).data,3);
    else
      mean_vals = ones(size(diff_vals));
    end;
    % calculate contrast ratio
    results.contrast.(parms.measlist{m}).data = ...
      -diff(results.cortsurf.(parms.measlist{m}).data,[],3) ./ mean_vals;
    results.contrast.(parms.measlist{m}).std = ...
      sqrt(sum(results.cortsurf.(parms.measlist{m}).std.^2,3)) ./ mean_vals;
    results.contrast.(parms.measlist{m}).nvals = ...
      mean(results.cortsurf.(parms.measlist{m}).nvals,3); % should be identical
    results.contrast.(parms.measlist{m}).nvals_valid = ...
      min(results.cortsurf.(parms.measlist{m}).nvals_valid,[],3);
    results.contrast.(parms.measlist{m}).nvals_invalid = ...
      max(results.cortsurf.(parms.measlist{m}).nvals_invalid,[],3);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function results = compile_cort_area(parms,results)
  parms.thick_flag = 0;
  parms.area_flag = 1;
  parms.cortvol_flag = 0;
  parms.cortT1w_flag = 0;
  parms.subvol_flag = 0;
  parms.subT1w_flag = 0;
  parms.aparc_flag = 1;
  parms.fuzzy_flag = 0;
  parms.concat_flag = 0;
  parms.subjinfo_flag = 0;
  args = mmil_parms2args(parms,parms.compile_area_tags);
  tmp_results = MMIL_Compile_MRI_Analysis(parms.ProjID,args{:});
  % handle potential for subjects wtih missing MRI data (e.g. incomplete recon)
  DTI_VisitIDs = results.STRUCT_VisitIDs;
  MRI_VisitIDs = tmp_results.VisitIDs;
  [tmp,ind_MRI,ind_DTI] = intersect(MRI_VisitIDs,DTI_VisitIDs);
  results.cort_area = tmp_results.cort_area;
  results.cort_area.data = nan(results.nsubs,results.cort_area.nrois);
  results.cort_area.data(ind_DTI,:) = tmp_results.cort_area.data(ind_MRI,:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_mean_cort_vals(parms,results,roitype)
  if ~exist('roitype','var') || isempty(roitype)
    roitype = 'cortsurf';
  end;
  % get roicodes from results
  roicodes = results.(roitype).roicodes;
  roinames = results.(roitype).roinames;
  codes_all = parms.aparc_roilist;
  codes_lh = codes_all(codes_all<2000);
  codes_rh = codes_all(codes_all>=2000);
  % identify ROIs from aparc
  [tmp,ind_lh]=intersect(roicodes,codes_lh);
  [tmp,ind_rh]=intersect(roicodes,codes_rh);
  [tmp,ind_all]=intersect(roicodes,codes_all);
  % calculate mean values for each measure
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    % prepare cortical area values for calculating weighted averages
    if parms.area_flag
      np = size(results.(roitype).(meas).data,3);
      % get area data
      area_lh = repmat(results.cort_area.data(:,ind_lh),[1,1,np]);
      area_rh = repmat(results.cort_area.data(:,ind_rh),[1,1,np]);
      area_all = repmat(results.cort_area.data(:,ind_all),[1,1,np]);
    else
      area_lh = ones(size(results.(roitype).(meas).data(:,ind_lh,:)));
      area_rh = ones(size(results.(roitype).(meas).data(:,ind_rh,:)));
      area_all = ones(size(results.(roitype).(meas).data(:,ind_all,:)));
    end;
    % calculate mean of ROI measures
    statlist = {'data','std'};
    for s=1:length(statlist)
      stat = statlist{s};
      % get data
      tmp_lh = results.(roitype).(meas).(stat)(:,ind_lh,:);
      tmp_rh = results.(roitype).(meas).(stat)(:,ind_rh,:);
      tmp_all = results.(roitype).(meas).(stat)(:,ind_all,:);
      % calculate mean value weighted by area
      tmp_lh = mmil_wtd_mean(tmp_lh,area_lh,2);
      tmp_rh = mmil_wtd_mean(tmp_rh,area_rh,2);
      tmp_all = mmil_wtd_mean(tmp_all,area_all,2);
      % concatenate to matrix
      results.(roitype).(meas).(stat) = ...
        cat(2,results.(roitype).(meas).(stat),tmp_lh,tmp_rh,tmp_all);
    end;
    % calculate sum of nvals
    statlist = {'nvals','nvals_valid','nvals_invalid'};
    for s=1:length(statlist)
      stat = statlist{s};
      % get data
      tmp_lh = results.(roitype).(meas).(stat)(:,ind_lh,:);
      tmp_rh = results.(roitype).(meas).(stat)(:,ind_rh,:);
      tmp_all = results.(roitype).(meas).(stat)(:,ind_all,:);
      % calculate sum
      tmp_lh = sum(tmp_lh,2);
      tmp_rh = sum(tmp_rh,2);
      tmp_all = sum(tmp_all,2);
      % concatenate to matrix
      results.(roitype).(meas).(stat) = ...
        cat(2,results.(roitype).(meas).(stat),tmp_lh,tmp_rh,tmp_all);
    end;
  end;
  % set roicodes
  new_roicodes = [1040;2040;2050];
  new_roinames = {'ctx-lh-mean';'ctx-rh-mean';'ctx-mean'};
  % concatenate to roicodes, and roinames
  results.(roitype).roicodes = cat(1,results.(roitype).roicodes,new_roicodes);
  results.(roitype).roinames = cat(1,results.(roitype).roinames,new_roinames);
  results.(roitype).nrois = length(results.(roitype).roicodes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = concat_results(parms,results)
  results = init_all_results(parms,results);
  roitypes = {};
  if parms.fiber_flag,    roitypes{end+1} = 'fiber';    end;
  if parms.aseg_flag,     roitypes{end+1} = 'aseg';     end;
  if parms.wmparc_flag,   roitypes{end+1} = 'wmparc';   end;
  if parms.cortsurf_flag, roitypes{end+1} = 'cortsurf'; end;
  if parms.contrast_flag,   roitypes{end+1} = 'contrast'; end;
  if parms.gwcsurf_flag, roitypes{end+1} = 'gwcsurf'; end;
  for t=1:length(roitypes)
    roitype = roitypes{t};
    roicodes = results.(roitype).roicodes;
    switch roitype
      case 'cortsurf'
        pvec = [1:parms.nprojdist];
      case 'gwcsurf'
        pvec = [1:parms.ngwlayers];
      otherwise
        pvec = 1;
    end;
    if strcmp(roitype,'fiber') && parms.fiber_vol_flag
      measlist = cat(2,parms.measlist,'vol');
    else
      measlist = parms.measlist;
    end;    
    nmeas = length(measlist);
    for p=pvec
      for m=1:nmeas
        meas = measlist{m};
        results.all.roicodes = cat(1,results.all.roicodes,roicodes);
        roinames = results.(roitype).roinames;
        for i=1:length(roinames)
          stem = [roitype '_' meas];
          switch roitype
            case 'cortsurf'
              projdist = parms.projdist_list(p);
              if projdist < 0
                pstr = sprintf('white%0.1f',projdist);
              else
                pstr = sprintf('gray+%0.1f',projdist);
              end;
              stem = [stem '_' pstr];
            case 'gwcsurf'
              switch p
                case 1
                  stem = [stem '_wm'];
                case 2
                  stem = [stem '_gm'];
                case 3
                  stem = [stem '_gwc'];
              end;
          end;
          roinames{i} = [stem '-' roinames{i}];
        end;
        results.all.roinames = cat(1,results.all.roinames,roinames);
        results.all.data = ...
          cat(2,results.all.data,results.(roitype).(meas).data(:,:,p));
        results.all.std = ...
          cat(2,results.all.std,results.(roitype).(meas).std(:,:,p));
        results.all.nvals = ...
          cat(2,results.all.nvals,results.(roitype).(meas).nvals(:,:,p));
        results.all.nvals_valid = ...
          cat(2,results.all.nvals_valid,results.(roitype).(meas).nvals_valid(:,:,p));
        results.all.nvals_invalid = ...
          cat(2,results.all.nvals_invalid,results.(roitype).(meas).nvals_invalid(:,:,p));
      end;
    end;
  end;
  results.all.nrois = length(results.all.roicodes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = load_roi_data(parms,s,meas,roitype,p)
  if ~exist('p','var') || isempty(p), p = 1; end;
  roi_data = [];
  % set input file name
  fstem = set_fstem(parms,s,meas);
  if isempty(fstem), return; end;
  infix = set_infix(parms,roitype,p);
  fname_in = sprintf('%s%s_%s.mat',fstem,meas,infix);
  if ~exist(fname_in,'file')
    if parms.verbose
      fprintf('%s: WARNING: missing %s for %s\n',...
        mfilename,fname_in,parms.StudyInfo(s).VisitID);
    end;
    return;
  end;
  % load ROI data
  if parms.verbose==2
    fprintf('%s: loading %s results for %s...\n',...
      mfilename,meas,parms.StudyInfo(s).VisitID);
  end;
  load(fname_in);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

