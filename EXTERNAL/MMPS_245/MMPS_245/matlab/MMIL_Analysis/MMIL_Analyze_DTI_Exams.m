function MMIL_Analyze_DTI_Exams(ProjID,varargin)
%function MMIL_Analyze_DTI_Exams(ProjID,[options])
%
% Usage:
%  MMIL_Analyze_DTI_Exams(ProjID,'key1', value1,...);
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
%  'batchrootdir': top level directory containing output batch job directories
%     {default = /home/$USER/batchdirs}
%  'batchname': name of output batchdir
%     {default = 'MMIL_Analyze_DTI_Exams'}
%
% Optional Parameters:
%   'outdir': output directory relative to ContainerPath
%     {default = 'DTanalysis'}
%   'fiber_flag': [0|1] extract white matter fiber ROI results
%     {default = 1}
%   'aseg_flag': [0|1] extract aseg ROI results
%     {default = 1}
%   'wmparc_flag': [0|1] extract wmparc ROI results
%     {default = 0}
%   'cortsurf_flag': [0|1] paint to cortical surface
%     and extract aparc ROI results
%     {default = 1}
%   'gwcsurf_flag': [0|1] paint to cortical surface for gray and white matter
%     and extract aparc ROI results
%     {default = 0}
%   'brainmask_flag': [0|1] create DTI resolution brain mask from aseg
%     {default = 1}
%   'measlist': list of DT "measures" to extract for each ROI
%      e.g. 'FA', 'MD', 'LD', 'TD', 'b0', 'b0N', 'T2w', 'T1w'
%      'b0' is the b=0 volume collected with diffusion scan
%      'b0N' is b=0 volume normalized by mean signal in ventricles
%      'T2w' is b=0 volume normalized by linear fit to MD
%      'T1w' is the T1-weighted FreeSurfer nu.mgz volume
%       use 'none' to extract none of these measures
%     {default = {'FA','MD','LD','TD','T2w','T1w'}}
%   'scalefacts': scaling factors applied to each measure in measlist
%     {default = [1 1000 1000 1000 1 1/256]}
%   'inputlist': list of input files in addition to or in place of measlist
%     relative to ContainerPath e.g. 'T1ratio.mgz'
%     if dimensions and vox2ras matrices do not match DTI,
%       they will be assumed to be in register with T1 volume
%       and resampled to DTI resolution for aseg and fiber analysis
%     {default = []}
%   'resT1flag': [0|1] conduct all analyses at T1 resolution
%     otherwise aseg and fiber analysis is done in DTI resolution
%     {default = 0}
%   'resT1_outfix': output string added to volumes resampled to T1 space
%     {default = 'resT1'}
%   'resDTI_outfix': output string added to volumes resampled to DTI space
%     {default = 'resDTI'}
%   'b0N_roilist': vector of aseg ROI numbers
%       to use to normalize b=0 image for b0N measure
%     {default = [4,43]} (lateral ventricles)
%   'b0N_erode_flag': [0|1] erode aseg ROIs used to normalize b=0 image
%     {default = 1}
%   'T2w_scalefact': additional scaling factor applied when creating T2w
%     {default = 1000}
%   'b0_norm_DTcalc_flag': [0|1] save T2w or b0N in DTcalc dir
%     otherwise in outdir
%     {default = 0}
%   'disp_scalefact': scaling factor applied to MD when
%     comparing values to dispersion values
%     relevant if aseg_disp_flag=1 or fiber_disp_flag=1
%     {default = 1000}
%   'minval': minimum value; if NaN, include all voxels
%     {default = 1e-6}
%   'full_fstem_flag': [0|1] include full stem for DT measures in output names
%     otherwise, use only meas name (e.g. 'FA' instead of 'DTI_scans_1_2_...FA')
%     {default = 0}
%   'DT_outdir': output directory for DT calculations
%     {default = 'DTcalc'}
%   'DT_outfix': string attached to DT calculation file names
%     {default = []}
%   'csv_flag': [0|1] output results in csv files
%     {default = 1}
%   'verbose': [0|1] display status meassages
%     {default = 1}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Optional Parameters specific to fiber tract ROIs:
%   'fibers': fiber numbers to include in analysis
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'atlas_flag': whether to use atlas fibers and if so,
%                 what type of atlas fibers
%     0 - manually assisted fiber tracts generated with DTIStudio
%       (DTIStudio_fiber_masks directory must exist
%        -- imported with Import_DTIStudio_FiberMasks)
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%       (fiber_masks_from_atlas directory must exist
%        -- created with AtlasTrack_Fibers)
%     If input supplied is a vector, will loop over each element and
%       run analysis with each atlas type (e.g. 'atlas_flag',[0,1,2])
%     {default = 2}
%   'weighted_avg_flag': [0|1] whether to calculated weighted averages for each ROI
%     using fiber counts or fiber probabilities to weight contribution from each voxel
%     {default = 1}
%   'thresh_prob': mask atlas-derived fiber ROIs by applying probability threshold
%     If input supplied is a vector, will loop over each element and
%       run analysis with each thresh_prob (e.g. 'thresh_prob',[0,0.1,0.2,0.3])
%     {default = 0}
%   'thresh_FA': mask fiber ROIs by applying threshold to FA image
%     If input supplied is a vector, will loop over each element and
%       run analysis with each atlas type (e.g. 'thresh_FA',[0,0.15])
%     {default = 0}
%   'fiber_disp_flag': [0|1] whether to calculate weighted averages
%     based on MD values and dispvec
%     {default = 0}
%   'fiber_dispvec': vector of MD ROI dispersion values (MAD estimates)
%     may be a single value used for all fibers
%       or length must match 'fibers'
%     if empty, will use fiber_disp_fname
%     {default = []}
%   'fiber_disp_fname': file name containing ROI dispersion results
%     i.e. output of DTI_MMIL_Calc_ROI_Dispersion
%     must exist if fiber_disp_flag = 1 and fiber_dispvec is empty
%     full path or relative to atlases or MetaData dir
%     {default = 'DTI_ROI_Dispersion/DTI_fiber_dispersion.mat'}
%   'fiber_disp_atlas_flag': [0|1] specify relative location of fiber_disp_fname
%     0: relative to MMPS atlases dir
%     1: relative to /home/{user}/MetaData/{ProjID}
%     {default = 1}
%   'fiber_disp_med_flag': [0|1] use median value from fiber_disp_fname
%     {default = 0}
%   'fiber_dispfact': multiplicative factor applied to fiber dispersion values
%     {default = 4}
%   'fiber_disp_xcg_flag': [0|1] interaction between disp weighting and xcg
%     0: calculate dispersion weighting everywhere, but mask out xcg voxels
%     1: calculate dispersion weighting only in xcg voxels
%     {default = 0}
%   'xcg_flag': [0|1] exclude CSF and gray-mattter from fiber ROIs
%     {default = 0}
%   'masksf_flag': [0|1] exclude voxels with multiple fibers
%     {default = 0}
%   'fseg_flag': [0|1] save fibers as fseg segmentation volume
%     {default = 0}
%   'fseg_resT1flag': [0|1] gemerate fiber segmentation in T1 resolution
%     but do not save individual resampled fibers (if resT1flag = 0)
%     ignored if resT1flag = 1
%     {default = 0}
%   'create_paths_flag': [0|1] to generate DTI Studio format fiber paths
%     {default = 0}
%   'atlastrack_dir': name of AtlasTrack output directory
%     relative to ContainerPath
%     {default = 'AtlasTrack'} 
%   'mapsdir': directory containing fiber probability map files
%     relative to atlastrack_dir
%     {default = 'fiber_maps'}
%   'pathsdir': directory containing fiber paths relative to atlastrack_dir
%     {default = 'fiber_paths'}
%   'fiber_atlasname': name of fiber atlas used, added to fiber_dir
%     if empty, default atlas used, no extra string included
%     {default=[])
%   'fiber_ATdir_flag': [0|1] save transformed fibers in atlastrack_dir
%     otherwise in outdir
%     {default = 0}
%
% Optional Parameters specific to aseg ROIs:
%   'aseg_disp_flag': [0|1] whether to calculate weighted averages
%     based on MD values and dispvec
%     {default = 0}
%   'aseg_disp_fname': file name containing ROI dispersion results
%     i.e. output of DTI_MMIL_Calc_ROI_Dispersion
%     must exist if aseg_disp_flag = 1
%     full path or relative to atlases or MetaData dir
%     {default = 'DTI_ROI_Dispersion/DTI_aseg_dispersion.mat'}
%   'aseg_disp_atlas_flag': [0|1] specify relative location of aseg_disp_fname
%     0: relative to MMPS atlases dir
%     1: relative to /home/{user}/MetaData/{ProjID}
%     {default = 1}
%   'aseg_dispfact': multiplicative factor applied to aseg dispersion values
%     {default = 4}
%   'erode_flag': [0|1] whether to "erode" ROIs
%     by smoothing and thresholding (to reduce edge effects)
%     {default = 1}
%   'erode_nvoxels': number of voxels to erode (integer)
%     {default = 1}
%   'aseg_aparc_flag': [0|1|2] whether to use cortical parcellation ROIs
%     0: aseg only
%     1: aparc only   NOTE: for best results with aparc, use erode_flag=0
%     2: aparc+aseg
%     {default = 0}
%
% Optional Parameters specific to wmparc ROIs:
%   'wmparc_aparc_flag': [0|1] whether to use cortical parcellation ROIs
%     0: wm only
%     1: wm and ctx
%     {default = 0}
%
% Optional Parameters specific to cortical surface:
%   'fnames_aparc': cell array of annotation files (one for each hemisphere)
%     relative to fspath/label
%     if empty, will use ?h.aparc.annot files in fspath/label
%     {default = []}
%   'cortsurf_resT1flag': [0|1] resample volume to T1 resolution before painting
%     implicitly 1 if resT1flag = 1
%     {default = 1}
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast maps
%     {default = [-1,1]}
%   'gwnorm_flag': [0|1] for gray/white contrast, normalize difference by mean
%     {default = 1}
%   'smoothsteps': number of smoothing iterations on individual subject surface
%     {default = 0}
%   'sphere_flag': [0|1] whether to resample to spherical atlas
%     {default = 0}
%   'sphsmoothsteps': number of smoothing iterations on sphere
%     {default = 0}
%
% Optional Parameters specific to gwcsurf analysis:
%   'gwcsurf_outfix': suffix added to output files for gwcsurf analysis
%     {default = 'gwcsurf'}
%   'pmin': minimum distance (mm) to project along surface vertex normal
%     {default = 0}
%   'pmax': maximum distance (mm) to project along surface vertex normal
%     {default = 5}
%   'pstep': step distance (mm) along surface vertex normal
%     {default = 1}
%   'interpmethod': interpolation method ('nearest','linear','spline','cubic')
%     {default = 'linear'}
%   'gwmask_flag': [0|1] use cortical ribbon volume for gray/white matter masks
%     to calculate weighted average
%     {default = 0}
%   'tukey_flag': use Tukey's bisquare function to transform weights
%     {default = 1}
%   'tukey_fact': scaling factor used to transform weights with Tukey's bisqaure
%     {default = 0.5}
%    NOTE: options for cortsurf analysis also apply to gwcsurf
%          all except for projdist_list
%
% Optional Parameters for resampling volumes to atlas space
%   'atlas_warp_flag': [0|1] nonlinearly resample input volumes to atlas space
%     {default = 0}
%   'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%   'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%
% Optional Parameters that determine input diffusion data:
%   (must correspond to existing DTcalc output)
%   'snums_flag': [0|1|2|3] which set of "snums" to use
%     0: use all available scans
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%     {default = 1}
%   'snum_index': index specifying which scan number of DTIScanNums
%     (or DTIScanNums2) to use in spread sheet (must have run DT
%     calculations separately for each scan)
%     If empty, use DT measures calculated from all DTIScanNums
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     e.g. 'corr_resDTI' or 'corr_regT1'
%     {default = []}
%   'auto_infix_flag': [0|1] set infix automatically based on typical
%     processing and settings in ProjInfo
%     ignored if infix is not empty
%     {default = 1}
%   'revflag': [0|1|2] specify whether to use sournon-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%   'min_ndirs': minimum number of gradient directions allowed
%     for tensor calculations
%     {default = 6}
%   'min_bval': minimum b-value allowed for tensor calculations
%     {default = 0}
%   'max_bval': maximum b-value used in tensor fit
%     {default = Inf}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'min_nb0': minimum number of b=0 images required for tensor calculations
%     {default = 1}
%   'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%
% NOTE: 'atlas_flag', 'thresh_FA', and 'thresh_prob' may optionally be supplied
%   as vectors.  If so, this function will loop over each vector, performing
%   analyses for all combinations.  For example, if each of them has 2 elements,
%   the fiber ROI analysis will be run each of 8 different ways (2*2*2).
%
% Created:  02/12/07 by Don Hagler
% Prev Mod: 07/31/16 by Don Hagler
% Last Mod: 08/03/17 by Don Hagler
%

%% todo: wmparc_disp_flag, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchrootdir',[],[],...
  'batchname','MMIL_Analyze_DTI_Exams',[],...
...
  'outdir','DTanalysis',[],...
  'fiber_flag',true,[false true],...
  'aseg_flag',true,[false true],...
  'wmparc_flag',false,[false true],...
  'cortsurf_flag',true,[false true],...
  'gwcsurf_flag',false,[false true],...
  'brainmask_flag',true,[false true],...
  'measlist',{'FA','MD','LD','TD','T2w','T1w'},[],...
  'scalefacts',[1,1000,1000,1000,1,1.0/256],[],...
  'inputlist',[],[],...
  'resT1flag',false,[false true],...
  'resT1_outfix','resT1',[],...
  'resDTI_outfix','resDTI',[],...
  'b0N_roilist',[4,43],[1,1000],...
  'b0N_erode_flag',true,[false true],...
  'T2w_scalefact',1000,[1e-10,1e10],...
  'b0_norm_DTcalc_flag',false,[false true],...
  'disp_scalefact',1000,[1e-10,1e10],...
  'full_fstem_flag',false,[false true],...
  'DT_outdir','DTcalc',[],...
  'DT_outfix',[],[],...
  'csv_flag',true,[false true],...
  'minval',1e-6,[],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
... % fiber parameters
  'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
  'atlas_flag',2,[],...
  'weighted_avg_flag',true,[false true],...
  'thresh_prob',0,[],...
  'thresh_FA',0,[],...
  'fiber_disp_flag',false,[false true],...
  'fiber_dispvec',[],[],...
  'fiber_disp_fname','DTI_ROI_Dispersion/DTI_fiber_dispersion.mat',[],...
  'fiber_disp_atlas_flag',true,[false true],...
  'fiber_disp_med_flag',false,[false true],...
  'fiber_dispfact',4,[1e-6,1e6],...
  'fiber_disp_xcg_flag',false,[false true],...
  'xcg_flag',false,[false true],...
  'masksf_flag',false,[false true],...
  'fseg_flag',false,[false true],...
  'fseg_resT1flag',false,[false true],...
  'create_paths_flag',false,[false true],...
  'atlastrack_dir','AtlasTrack',[],...
  'mapsdir','fiber_maps',[],...
  'pathsdir','fiber_paths',[],...
  'fiber_atlasname',[],[],...
  'fiber_ATdir_flag',false,[false true],...
... % aseg parameters
  'aseg_disp_flag',false,[false true],...
  'aseg_disp_fname','DTI_ROI_Dispersion/DTI_aseg_dispersion.mat',[],...
  'aseg_disp_atlas_flag',true,[false true],...
  'aseg_dispfact',4,[1e-6,1e6],...
  'erode_flag',true,[false true],...
  'erode_nvoxels',1,[1:100],...
  'aseg_aparc_flag',0,[0,1,2],...
... % wmparc parameters
  'wmparc_aparc_flag',false,[false true],...    
... % cortsurf parameters
  'fnames_aparc',[],[],...
  'cortsurf_resT1flag',true,[false true],...
  'projdist_list',[-1,1],[-5,5],...
  'gwnorm_flag',true,[false true],...
  'smoothsteps',0,[0,Inf],...
  'sphere_flag',false,[false true],...
  'sphsmoothsteps',0,[0,Inf],...
  ... % gwcsurf parameters
  'gwcsurf_outfix','gwcsurf',[],...
  'pmin',0,[0,10],...
  'pmax',5,[0,10],...
  'pstep',1,[0.001,10],...
  'interpmethod','linear',{'nearest','linear','spline','cubic'},...
  'gwmask_flag',false,[false true],...
  'tukey_flag',true,[false true],...
  'tukey_fact',0.5,[1e-10,1e10],...
... % resample atlas parameters
  'atlas_warp_flag',false,[false true],...
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],...    
... % data parameters
  'snums_flag',1,[0:3],...
  'snum_index',[],[],...
  'infix',[],[],...
  'auto_infix_flag',true,[false true],...
  'revflag',0,[0,1,2],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'max_bval',Inf,[100,Inf],...
  'flex_flag',false,[false true],...  
  'min_nb0',1,[],...
  'nob0_flag',false,[false true],...
... % hidden parameters
  'fnum',1,[1,Inf],...
  'RSI_outdir','RSIcalc',[],...
  'RSI_outfix',[],[],...
  'aseg_disp_suffix','dwtd',[],...
  'fiber_disp_suffix','dwtd',[],...
  'xcg_suffix','xcg',[],...
  'xcg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
  'masksf_suffix','masksf',[],...
  'manual_mapsdir','DTIStudio_fiber_masks',[],...
  'manual_pathsdir','DTIStudio_fiber_paths',[],...
  'qcflag',true,[false true],...
  'required_rootdirs',{'proc_dti','fsurf'},[],...
  'QC_raw',true,[false true],...
  'QC_DTI',true,[false true],...
  'QC_recon',true,[false true],...
  'resample_flag',true,[false true],...
  'regT1flag',1,[0:2],...
...
  'excl_tags',{'excl_tags','ProjInfo_tags',...
               'StudyInfo','RootDirs','auto_infix_flag','batchrootdir','batchname',...
               'snums_flag','snum_index','qcflag','required_rootdirs',...
               'QC_raw','QC_DTI','QC_recon',...
               'resample_flag','regT1flag',...
               'aseg_disp_fname','aseg_disp_atlas_flag'...
               'fiber_disp_fname','fiber_disp_atlas_flag','fiber_disp_med_flag'},[],...
  'ProjInfo_tags',{'min_bval','max_bval','flex_flag','xcg_flag','masksf_flag',...
                   'fseg_flag','fseg_resT1flag','wmparc_flag','cortsurf_flag',...
                   'gwcsurf_flag','atlas_flag','revflag','snums_flag','fiber_atlasname',...
                   'resample_flag','regT1flag','RSI_outdir','RSI_outfix',...
                   'DT_outdir','DT_outfix'},[],...
};
parms = mmil_args2parms(varargin,parms_filter);



args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args

  % get all 'DTI_' parms, strip 'DTI_'
  ProjInfo_args = MMIL_Args(ProjInfo,'DTI');
  % convert from args cell array to parms struct
  ProjInfo_parms = mmil_args2parms(ProjInfo_args,[],0);
  % keep only those in ProjInfo_tags
  ProjInfo_args = mmil_parms2args(ProjInfo_parms,parms.ProjInfo_tags);
  % merge arguments, giving varargin precedence over ProjInfo_args
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

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

% get aseg dispersion values
if parms.aseg_flag && parms.aseg_disp_flag
  if mmil_isrelative(parms.aseg_disp_fname)
    if parms.aseg_disp_atlas_flag
      parms.aseg_disp_fname = ...
        [getenv('MMPS_DIR') '/atlases/' parms.aseg_disp_fname];
    else
      parms.aseg_disp_fname = ...
        [getenv('HOME') '/MetaData/' ProjID '/' parms.aseg_disp_fname];
    end;
  end;
  if ~exist(parms.aseg_disp_fname)
    error('file %s not found (required because aseg_disp_flag = 1)');
  end;
  results = [];
  load(parms.aseg_disp_fname);
  if isempty(results)
    error('aseg dispersion results not found in %s',parms.aseg_disp_fname);
  end;
  parms.aseg_dispvec = results.dispersion;
  parms.aseg_disp_roicodes = results.roicodes;
end;

% get aseg dispersion values
if parms.fiber_flag && parms.fiber_disp_flag && isempty(parms.fiber_dispvec)
  if mmil_isrelative(parms.fiber_disp_fname)
    if parms.fiber_disp_atlas_flag
      parms.fiber_disp_fname = ...
        [getenv('MMPS_DIR') '/atlases/' parms.fiber_disp_fname];
    else
      parms.fiber_disp_fname = ...
        [getenv('HOME') '/MetaData/' ProjID '/' parms.fiber_disp_fname];
    end;
  end;
  if ~exist(parms.fiber_disp_fname)
    error('file %s not found (required because fiber_disp_flag = 1)');
  end;
  results = [];
  load(parms.fiber_disp_fname);
  if isempty(results)
    error('fiber dispersion results not found in %s',parms.fiber_disp_fname);
  end;
  dispvec = results.dispersion;
  med_disp = median(dispvec);
  if parms.fiber_disp_med_flag
    parms.fiber_dispvec = med_disp;
  else
    % match disp values for selected fibers
    nfibers = length(parms.fibers);
    parms.fiber_dispvec = zeros(1,nfibers);
    roicodes = results.roicodes - 1e4;
    [tmp,ind_fibers,ind_disp] = intersect(parms.fibers,roicodes);
    ind_missing = setdiff(1:nfibers,ind_fibers);
    if isempty(ind_fibers)
      fprintf('%s: WARNING: fibers do not match roicodes in fiber_disp_fname %s\n',...
        mfilename,parms.fiber_disp_fname);
    else
      parms.fiber_dispvec(ind_fibers) = dispvec(ind_disp);  
    end;
    if ~isempty(ind_missing)
      fprintf('%s: WARNING: fibers %s not found in fiber_disp_fname %s, using median value %0.3f\n',...
        mfilename,sprintf('%d ',parms.fibers(ind_missing)),...
        parms.fiber_disp_fname,med_disp);
      parms.fiber_dispvec(ind_missing) = med_disp;
    end;
  end;
end;

if ~isempty(parms.inputlist)
  DTI_flag = 0; % do not require DTI data
else
  DTI_flag = 1; % require DTI data
end;

% set DTI-specific fields in StudyInfo
StudyInfo = DTI_MMIL_Check_StudyInfo(...
  StudyInfo,RootDirs,...
  'snums_flag',parms.snums_flag,... % DTIScanNums field will get changed
  'qcflag',parms.qcflag,...
  'DTI_flag',DTI_flag);
if isempty(StudyInfo), error('empty StudyInfo after check for DTI'); end;

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;

if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

% create jobs for each subject in StudyInfo
j = 1;
for i=1:length(StudyInfo)
  VisitID = StudyInfo(i).VisitID;
  tmp_parms = parms;
  tmp_parms.snums = StudyInfo(i).DTIScanNums;
  ContainerPath = sprintf('%s/%s',RootDirs.proc_dti,StudyInfo(i).proc_dti);
  FSContainerPath = sprintf('%s/%s',RootDirs.fsurf,StudyInfo(i).fsurf);
  if isempty(RootDirs.fsurf) || isempty(StudyInfo(i).fsurf) || ~exist(FSContainerPath,'dir')
    FSContainerPath = [];
  end;
  if isempty(StudyInfo(i).proc_dti) || ~exist(ContainerPath,'file')
    continue;
  end;
  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = sprintf('%s/%s.m',batchdir,jobID);
  tags = setdiff(fieldnames(tmp_parms),parms.excl_tags);
  mmil_write_script(jobfname,'MMIL_Analyze_DTI_Exam',...
    {ContainerPath,FSContainerPath},tags,tmp_parms);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

% check available disk space
MMIL_Check_Usage(RootDirs.proc_dti);

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);


