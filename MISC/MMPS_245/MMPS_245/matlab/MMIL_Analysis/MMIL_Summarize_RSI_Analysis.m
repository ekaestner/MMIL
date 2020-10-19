function results = MMIL_Summarize_RSI_Analysis(ProjID,varargin)
% function results = MMIL_Summarize_RSI_Analysis(ProjID,[options])
%
% Purpose: create summary spreadsheets with RSI ROI measures
%
% Usage:
%  results = MMIL_Summarize_RSI_Analysis(ProjID,'key1', value1,...);
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
%   'measlist': list of RSI "measures" to extract for each ROI
%       T  : norm of all parameters
%       F0 : 0th order FOD
%       N0 : F0 normalized by T
%       F2 : norm of 2nd order FOD
%       N2 : F2 normalized by T
%       F4 : norm of 4th order FOD
%       N4 : F4 normalized by T
%       FD : norm of directional FOD components (F2 and F4 combined)
%       ND : FD normalized by T (a.k.a. "neurite density")
%       FT : norm of total FOD (F0, F2, and F4 combined)
%       NT : FT normalized by T (a.k.a. "volume fraction")
%       Ir : restricted isotropic
%       NIr: Ir normalized by T (restricted isotropic volume fraction)
%       Ih : hindered isotropic
%       NIh: Ih normalized by T (hindered isotropic volume fraction)
%       If : free isotropic
%       NIf: If normalized by T (free isotropic volume fraction)
%       AU : angular uncertainty
%     {default = {'T','F0','N0','F2','N2','F4','N4','FD','ND',...
%                 'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'}}
%   'scalefacts': scaling factors applied to each measure in measlist
%     if empty, will use 1 for each
%     if not empty, length must match length of measlist
%     {default = []}
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
%   'wmparc_flag': [0|1] extract wmparc ROI results
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
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%  'min_ndirs': minimum number of diffusion directions to be included
%     {default = 6}
%  'full_fstem_flag': [0|1] full stem included in RSI measures analysis output
%     otherwise, use only meas name (e.g. 'F2' instead of 'DTI_scans_1_2_...F2')
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
%   'RSI_analdir': name of RSI analysis folder relative to ContainerPath
%     {default = 'RSIanalysis'}
%   'RSI_outfix': string attached to RSI calculation and analysis file names
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
%     {default = 'RSI'}
%   'save_mat_flag': [0|1] save compiled ROI data in a mat file
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  02/06/13 by Don Hagler
% Last Mod: 04/14/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
results = [];

% check input parameters
parms = check_input(ProjID,varargin);

% summarize data
args = mmil_parms2args(parms,parms.tags);
MMIL_Summarize_DTI_Analysis(ProjID,args{:});

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
    'measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
                'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'},[],...
    'scalefacts',[],[],...
    'fiber_flag',true,[false true],...
    'fiber_vol_flag',false,[false true],...
    'aseg_flag',true,[false true],...
    'wmparc_flag',false,[false true],...
    'cortsurf_flag',true,[false true],...
    'gwcsurf_flag',false,[false true],...
    'concat_flag',1,[0:2],...
    'concat_outfix','all',[],...
... % specify DTI data used for RSI calculations and analysis
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
    'RSI_analdir','RSIanalysis',[],...
    'RSI_outfix',[],[],...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'verbose',1,[0:2],...
    'outdir','ROI_Summaries',[],...
    'outstem','RSI',[],...
    'save_mat_flag',false,[false true],...
    'forceflag',false,[false true],...
... % hidden
    'fnum',1,[1,Inf],...
    'DT_analdir','DTanalysis',[],...
    'DT_outfix',[],[],...
    'DT_outdir','DTcalc',[],...
    'RSI_outdir','RSIcalc',[],...
    'aseg_disp_suffix','dwtd',[],...
    'fiber_disp_suffix','dwtd',[],...
    'xcg_suffix','xcg',[],...
    'masksf_suffix','masksf',[],...
    'resample_flag',true,[false true],...
    'regT1flag',1,[0:2],...
...
    'motion_tags',{'mean_motion','mean_trans','mean_rot'},[],...
    'T2w_tags',{'T2w_sf','T2w_r'},[],...
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'ProjInfo_tags',{'min_bval','flex_flag',...
                     'xcg_flag','masksf_flag','fseg_flag',...
                     'wmparc_flag','cortsurf_flag','gwcsurf_flag',...
                     'atlas_flag','revflag','snums_flag',...
                     'fiber_atlasname','resample_flag','regT1flag',...
                     'RSI_outdir','RSI_outfix','DT_outdir','DT_outfix'},[],...
    'excl_tags',{'excl_tags','ProjInfo_tags','ProjID'},[],...
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

  parms.tags = setdiff(fieldnames(parms),parms.excl_tags);
return;

