function errcode = MMIL_Analyze_DTI_Exam(ContainerPath,FSContainerPath,varargin)
%function errcode = MMIL_Analyze_DTI_Exam(ContainerPath,FSContainerPath,[options])
%
% Purpose:
%  Extract ROI averages of diffusion tensor (DT) measures
%    like FA, MD, etc.
%  ROIs include fiber tracts (from DTI Studio or probabilistic atlas)
%    and those derived from freesurfer's subcortical segmentation (aseg)
%    and cortical surface parcellation (aparc)
%
% Usage:
%  MMIL_Analyze_DTI_Exam(ContainerPath,FSContainerPath,'key1', value1,...)
%
% Required Parameters:
%   ContainerPath: full path of directory containing processed
%     diffusion data (mgh format) and fiber tract ROI files
%   FSContainerPath: full path of directory containing freesurfer recon
%
% Optional Parameters:
%   'outdir': output directory
%     full path or relative to ContainerPath
%     {default = 'DTanalysis'}
%   'fiber_flag': [0|1] extract white matter fiber ROI results
%     {default = 1}
%   'aseg_flag': [0|1] extract aseg ROI results
%     {default = 1}
%   'wmparc_flag': [0|1] extract wmparc ROI results
%     {default = 0}
%   'cortsurf_flag': [0|1] paint to cortical surface
%     and extract aparc ROI results
%     {default = 0}
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
%     can be full path or relative to ContainerPath e.g. 'T1ratio.mgz'
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
%   'minval': minimum value; if NaN, include all voxels
%     {default = 1e-6}
%   'full_fstem_flag': [0|1] include full stem for DT measures in output names
%     otherwise, use only meas name (e.g. 'FA' instead of 'DTI_scans_1_2_...FA')
%     use this to compare results for different snums, min_bval, etc.
%     {default = 0}
%   'csv_flag': [0|1] output results in csv files
%     {default = 1}
%   'verbose': [0|1] display status meassages
%     {default = 0}
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
%     may be a single value used for all fibers or length must match 'fibers'
%     {default = 0.1}
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
%     absolute or else relative to ContainerPath
%     {default = 'AtlasTrack'} 
%   'mapsdir': directory containing fiber probability map files
%     relative to atlastrack_dir
%     {default = 'fiber_maps'}
%   'pathsdir': directory containing fiber paths relative to atlastrack_dir
%     {default = 'fiber_paths'}
%   'fiber_atlasname': name of fiber atlas used, added to mapsdir and pathsdir
%     if empty, default atlas used, no extra string included
%     {default = [])
%   'fiber_ATdir_flag': [0|1] save transformed fibers in atlastrack_dir
%     otherwise in outdir
%     {default = 0}
%
% Optional Parameters specific to aseg ROIs:
%   'aseg_disp_flag': [0|1] whether to calculate weighted averages
%     based on MD values and dispvec
%     {default = 0}
%   'aseg_dispvec': vector of MD ROI dispersion values (MAD estimates)
%     required if aseg_disp_flag = 1
%     {default = []}
%   'aseg_disp_roicodes': vector of ROI codes corresponding to dispvec elements
%     required if aseg_disp_flag = 1
%     {default = []}
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
%     if empty, will use ?h.aparc.annot files in fspath/label
%     if not full path, assumed to be relative to fspath/label
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
%   'sphere_flag': [0|1[ whether to resample to spherical atlas
%     {default = 0}
%   'sphsmoothsteps': number of smoothing iterations on sphere
%     {default = 0}
%
% Optional Parameters specific to gwcsurf analysis:
%   'gwcsurf_outfix': suffix added to output files for gwcsurf analysis
%     {default = []}
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
%   'snums': list of scan numbers used for tensor calculations
%     if empty (or unspecified), used all DTI scans in container
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     e.g. 'corr_resDTI' or 'corr_regT1'
%     {default = []}
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
% Prev Mod: 08/03/17 by Don Hagler
% Last Mod: 10/19/17 by Don Hagler
%

%% todo: wmparc_disp_flag, etc.
%% todo: option to use FA for warp to atlas?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
errcode = 0;

[parms,errcode] = check_input(ContainerPath,FSContainerPath,varargin);
if errcode, return; end;

if ~parms.fiber_flag && ~parms.aseg_flag &&...
   ~parms.wmparc_flag &&...
   ~parms.cortsurf_flag &&...
   ~parms.gwcsurf_flag &&...
   ~parms.brainmask_flag &&...
   ~parms.atlas_warp_flag
  fprintf('%s: WARNING: nothing to do\n',mfilename);
  return;
end;
mmil_mkdir(parms.outdir);

if parms.brainmask_flag, create_brainmask(parms); end;

if parms.fiber_flag, parms = transform_fibers(parms); end;

[parms,errcode] = prep_data(parms);
if errcode, return; end;

errcode = analyze_data(parms);

if parms.atlas_warp_flag
  errcode = atlas_warp_data(parms);
end;

if errcode
  fprintf('\n%s: WARNING: one or more analysis steps failed\n\n',mfilename);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_input(ContainerPath,FSContainerPath,options)
  errcode = 0;
  parms = mmil_args2parms(options,{...
    'cpath',ContainerPath,[],...
    'fspath',FSContainerPath,[],...
  ...
    'outdir','DTanalysis',[],...
    'fiber_flag',true,[false true],...
    'aseg_flag',true,[false true],...
    'wmparc_flag',false,[false true],...
    'cortsurf_flag',false,[false true],...
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
    'csv_flag',true,[false true],...
    'minval',1e-6,[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % fiber parameters
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'atlas_flag',2,[],...
    'weighted_avg_flag',true,[false true],...
    'thresh_prob',0,[],...
    'thresh_FA',0,[],...
    'fiber_disp_flag',false,[false true],...
    'fiber_dispvec',0.1,[],...
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
    'aseg_dispvec',[],[],...
    'aseg_disp_roicodes',[],[],...
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
    'gwcsurf_outfix',[],[],...
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
    'snums',[],[],...
    'infix',[],[],...
    'revflag',0,[0,1,2],...
    'min_ndirs',6,[],...
    'min_bval',1,[],...
    'max_bval',Inf,[100,Inf],...
    'flex_flag',false,[false true],...
    'min_nb0',1,[],...
    'nob0_flag',false,[false true],...
  ... % hidden parameters
    'fnamestem','DTI',[],...
    'fnum',1,[1,Inf],...
    'DT_outdir','DTcalc',[],...
    'DT_outfix',[],[],...
    'DT_measlist',{'FA','MD','LD','TD','b0','b0N','T2w'},[],...
    'RSI_outdir','RSIcalc',[],...
    'RSI_outfix',[],[],...
    'RSI_measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
                    'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'},[],...
    'multiframe_measlist',{'F0','N0','F2','N2','F4','N4',...
                           'FD','ND','FT','NT'},[],...
    'manual_mapsdir','DTIStudio_fiber_masks',[],...
    'manual_pathsdir','DTIStudio_fiber_paths',[],...
    'fseg_fstem','fseg',[],...
    'aseg_disp_suffix','dwtd',[],...
    'fiber_disp_suffix','dwtd',[],...
    'xcg_suffix','xcg',[],...
    'xcg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
    'masksf_suffix','masksf',[],...
    'count_flag',true,[false true],... % applies to manual fibers only
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'fname_colorlut',[],[],...
  ... % sets of parameter names to pass to other functions
    'fstem_tags',{'snums','infix','revflag','min_bval','max_bval','flex_flag',...
                  'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
    'rsi_fstem_tags',{'snums','infix','revflag','min_bval','flex_flag',...
                      'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
    'reg_tags',{'infix','revflag'},[],...
    'aseg_tags',{'outdir','outstem','fname_out','csv_flag','fname_aseg',...
                 'aseg_aparc_flag','fname_vals','dispvec',...
                 'disp_roicodes','disp_scalefact','dispfact','disp_suffix',...
                 'erode_flag','erode_nvoxels',...
                 'scalefact','frames','minval','M_reg','res_outfix',...
                 'fname_colorlut','verbose','forceflag',...
                 'aseg_roilist','aparc_roilist','exclude_roilist'},[],...
    'mask_tags',{'outdir','M_T1_to_DTI','fname_DTI','volsz_DTI','M_DTI',...
                 'smooth1','thresh1','smooth2','thresh2','smooth3',...
                 'nii_flag','nii_out_orient','xcg_flag','xcg_suffix',...
                 'xcg_codes',...
                 'res_outfix','interpm','forceflag'},[],...
    'fiber_tags',{'outdir','outstem','csv_flag','fibers','atlas_flag',...
                  'resT1flag','disp_flag','disp_suffix','dispfact',...
                  'xcg_flag','xcg_suffix','masksf_flag','masksf_suffix',...
                  'fiber_atlasname','fiber_infix','fiber_ext','M_reg',...
                  'res_outfix','fname_FA','thresh_FA','thresh_prob',...
                  'weighted_avg_flag','scalefact','minval','verbose',...
                  'forceflag','outfix','frames',},[],...
    'cortsurf_tags',{'outdir','outstem','fnames_aparc','fnames_weights',...
                     'weights_thresh','csv_flag','M_reg','resT1flag',...
                     'res_outfix','projdist_list',...
                     'gwnorm_flag','smoothsteps','sphere_flag',...
                     'sphsmoothsteps',...
                     'mask_midbrain_flag',...
                     'outtype','scalefact',...
                     'minval','fname_colorlut',...
                     'verbose','forceflag','frames'},[],...
    'gwcsurf_tags',{'outdir','outstem','fnames_aparc','fnames_weights',...
                    'weights_thresh','csv_flag','M_reg','resT1flag',...
                    'res_outfix','pmin','pmax','pstep','interpmethod',...
                    'avg_flag','gwmask_flag','tukey_flag','tukey_fact',...
                    'gwnorm_flag','smoothsteps','sphere_flag',...
                    'sphsmoothsteps','outtype','scalefact','minval',...
                    'fname_colorlut','verbose','forceflag','frames',...
                    'wm_outstem','gm_outstem'},[],...
    'transform_tags',{'outdir','fiber_outdir','atlas_flag','fibers',...
                      'resT1flag','M_T1_to_DTI','volsz_T1','M_T1',...
                      'save_mgz_flag','verbose','forceflag',...
                      'xcg_flag','fname_aseg','xcg_suffix',...
                      'masksf_flag','masksf_suffix',...
                      'disp_flag','fname_vals','disp_suffix',...
                      'disp_scalefact','dispvec','dispfact','disp_xcg_flag',...
                      'thresh_FA','fname_FA','thresh_prob',...
                      'fseg_flag','fseg_fibers','fseg_resT1flag',...
                      'fseg_xcg_flag','fseg_thresh_prob','fname_fseg',...
                      'create_paths_flag','fname_V0','min_fiberlen',...
                      'thresh_angle','path_suffix','paths_outdir',... 
                      'xcg_codes','roicode_base','count_flag',...
                      'prob_exponent','orient_ref'},[],...
    'warp_tags',{'atlasdir','atlasname'...
                 'smoothflag','sampling','nK','tstep','astep',...
                 'scales','ns','sf','thresh','stdflag','maskflag',...
                 'stdbgval','fname_reg','fname_T1',...
                 'outdir','outstem','outstem_T1',...
                 'interpm','bclamp','padding','vxlmap_flag','forceflag'},[],...
    'wmparc_code_tags',{'wmparc_aparc_flag','wmparc_roilist','aparc_roilist',...
                        'exclude_roilist'},[],...
  });

  % check that number of measures matches number of scalefacts
  if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;
  if strcmp(parms.measlist{1},'none')
    parms.measlist = [];
  elseif length(parms.scalefacts)==1
    parms.scalefacts = parms.scalefacts*ones(length(parms.measlist),1);
  elseif length(parms.measlist)~=length(parms.scalefacts)
    error('number of measures in measlist does not match number of scalefacts');
  end;

  % check that inputlist is a cell array
  if ~isempty(parms.inputlist) && ~iscell(parms.inputlist)
    parms.inputlist = {parms.inputlist};
  end;

  % check that processed data container exists
  if ~exist(parms.cpath,'dir')
    fprintf('%s: ERROR: %s not found\n',mfilename,parms.cpath);
    errcode = 1; return;
  end;

  % check that Freesurfer recon exists
  if ~exist(parms.fspath,'dir')
    fprintf('%s: ERROR: %s not found\n',mfilename,parms.fspath);
    errcode = 1; return;
  end

  % check DTcalc input files
  if ~isempty(intersect(parms.measlist,parms.DT_measlist)) ||...
     any(parms.thresh_FA~=0) ||...
     parms.aseg_disp_flag || parms.fiber_disp_flag
    parms = check_DTcalc(parms);
  end;

  % check RSIcalc input files
  if ~isempty(intersect(parms.measlist,parms.RSI_measlist))
    parms = check_RSIcalc(parms);
  end;

  % check input files for measlist and inputlist
  [parms,errcode] = check_input_files(parms);
  if errcode, return; end;
  
  % load registration information
  [parms,errcode] = load_RegInfo(parms);
  if errcode, return; end;

  % check fibers
  if parms.fiber_flag
    parms = check_orig_fibers(parms);
  end;

  % check that aseg exists
  if parms.aseg_flag || parms.brainmask_flag ||...
     (parms.fiber_flag && parms.xcg_flag)
    parms = check_aseg(parms);
  end;

  % check that wmparc exists
  if parms.wmparc_flag
    parms = check_wmparc(parms);
  end;

  % check cortsurf
  if parms.cortsurf_flag
    parms = check_cortsurf(parms);
  end;

  % check gwcsurf
  if parms.gwcsurf_flag
    parms = check_gwcsurf(parms);
  end;

  % set outdir
  if mmil_isrelative(parms.outdir)
    parms.outdir = [parms.cpath '/' parms.outdir];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_DTcalc(parms)
  invalid_flag = 0;
  if parms.verbose
    fprintf('%s: checking DT calculations...\n',mfilename);
  end;
  % check that DT calculations exist
  tmp_parms = parms;
  tmp_parms.outdir = parms.DT_outdir;
  tmp_parms.outfix = parms.DT_outfix;
  args = mmil_parms2args(tmp_parms,parms.fstem_tags);
  [parms.full_DT_fstem,parms.snums] = ...
    DTI_MMIL_Set_DT_fstem(parms.cpath,args{:});
  if isempty(parms.snums)
    fprintf('%s: WARNING: no valid DTI scans\n',mfilename);
    parms.measlist = setdiff(parms.measlist,parms.DT_measlist);
    invalid_flag = 1;
  end;
  if ~invalid_flag
    fname_in = [parms.full_DT_fstem '_meas.mat'];
    if ~exist(fname_in,'file')
      fprintf('%s: WARNING: %s not found\n',mfilename,fname_in);
      parms.measlist = setdiff(parms.measlist,parms.DT_measlist);
      invalid_flag = 1;
    end;
  end;
  if invalid_flag
    if parms.aseg_disp_flag
      parms.aseg_flag = 0;
    end;
    if parms.fiber_disp_flag
      parms.fiber_flag = 0;
    end;
    return;
  end;
  [tmp,parms.DT_fstem] = fileparts(parms.full_DT_fstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_RSIcalc(parms)
  if parms.verbose
    fprintf('%s: checking RSI calculations...\n',mfilename);
  end;
  % check that RSI calculations exist
  tmp_parms = parms;
  tmp_parms.outdir = parms.RSI_outdir;
  tmp_parms.outfix = parms.RSI_outfix;
  args = mmil_parms2args(tmp_parms,parms.rsi_fstem_tags);
  [parms.full_RSI_fstem,parms.snums] = ...
    DTI_MMIL_Set_RSI_fstem(parms.cpath,args{:});
  if isempty(parms.snums)
    fprintf('%s: WARNING: no valid DTI scans\n',mfilename);
    parms.measlist = setdiff(parms.measlist,parms.RSI_measlist);
    return;
  end;
  fname_in = [parms.full_RSI_fstem '_meas.mat'];
  if ~exist(fname_in,'file')
    fprintf('%s: WARNING: %s not found\n',mfilename,fname_in);
    parms.measlist = setdiff(parms.measlist,parms.RSI_measlist);
    return;
  end;
  [tmp,parms.RSI_fstem] = fileparts(parms.full_RSI_fstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_input_files(parms)
  errcode = 0;
  measlist = parms.measlist;
  scalefacts = parms.scalefacts;
  parms.fnamelist = [];
  parms.measlist = [];
  parms.full_measlist = [];
  parms.scalefacts = [];
  parms.fnumlist = [];
  % check files for each DT measure
  j = 1;
  for i=1:length(measlist)
    meas = measlist{i};
    scalefact = scalefacts(i);
    if ismember(meas,parms.multiframe_measlist)
      fnum = parms.fnum;
    else
      fnum = 1;
    end;
    for f=fnum
      full_meas = meas;
      switch meas
        case {'T1','T1w'}
          fname = sprintf('%s/mri/nu.mgz',parms.fspath);
        case {'T2w','b0N'}
          fname = sprintf('%s_b0.mgz',parms.full_DT_fstem);
          if parms.full_fstem_flag
            full_meas = [parms.DT_fstem '_' meas];
          else
            if ~isempty(parms.DT_outfix)
              full_meas = [parms.DT_outfix '_' meas];
            end;
            if parms.flex_flag
              full_meas = ['flex_' full_meas];
            end;
          end;
        case parms.DT_measlist
          fname = sprintf('%s_%s.mgz',parms.full_DT_fstem,meas);
          if parms.full_fstem_flag
            full_meas = [parms.DT_fstem '_' meas];
          else
            if ~isempty(parms.DT_outfix)
              full_meas = [parms.DT_outfix '_' meas];
            end;
            if parms.flex_flag
              full_meas = ['flex_' full_meas];
            end;
          end;
        case parms.RSI_measlist
          fname = sprintf('%s_%s.mgz',parms.full_RSI_fstem,meas);
          if parms.full_fstem_flag
            full_meas = [parms.RSI_fstem '_' meas];
          else
            if ~isempty(parms.RSI_outfix)
              full_meas = [parms.RSI_outfix '_' meas];
            end;
            if parms.flex_flag
              full_meas = ['flex_' full_meas];
            end;
          end;
        otherwise
          fprintf('%s: WARNING: invalid measure %s\n',mfilename,meas);
          continue;
      end;
      if f>1
        full_meas = sprintf('%ss%d',full_meas,f);
      end;
      if exist(fname,'file')
        parms.fnamelist{j} = fname;
        parms.measlist{j} = meas;
        parms.full_measlist{j} = full_meas;
        parms.scalefacts(j) = scalefact;
        parms.fnumlist(j) = f;
        j = j + 1;
      else
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
        continue;
      end;
    end;
  end;
  % check files for each input file
  for i=1:length(parms.inputlist)
    fname = parms.inputlist{i};
    if isempty(fname)
      fprintf('%s: WARNING: inputlist{%d} is empty\n',mfilename,i);
      continue;
    end;
    if mmil_isrelative(fname), fname = [parms.cpath '/' fname]; end;
    if exist(fname,'file')
      [tmp,fstem] = fileparts(fname);
      parms.fnamelist{j} = fname;
      parms.measlist{j} = fstem;
      parms.full_measlist{j} = fstem;
      parms.scalefacts(j) = 1;
      parms.fnumlist(j) = 1;
      j = j + 1;
    else
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
    end;
  end;
  parms.ninputs = length(parms.fnamelist);
  if parms.ninputs==0
    fprintf('%s: ERROR: no valid input files found\n',mfilename);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_orig_fibers(parms)
  if parms.verbose
    fprintf('%s: checking original fibers...\n',mfilename);
  end;
  if mmil_isrelative(parms.atlastrack_dir)
    parms.atlastrack_dir = [parms.cpath '/' parms.atlastrack_dir];
  end;
  for atlas_flag = parms.atlas_flag
    % set fiber_dir, check that it exists
    if atlas_flag
      fiber_dir = [parms.atlastrack_dir '/' parms.mapsdir];
      if ~isempty(parms.fiber_atlasname)
        fiber_dir = [fiber_dir '_' parms.fiber_atlasname];
      end;
    else
      fiber_dir = [parms.cpath '/' parms.manual_mapsdir];
    end;
    if ~exist(fiber_dir,'dir')
      fprintf('\n%s: WARNING: fiber_dir %s not found\n',...
        mfilename,fiber_dir);
      parms.fiber_flag = 0;
    else
      [all_exist_flag,some_exist_flag,fiber_infix,fiber_ext] = ...
        check_fibers(fiber_dir,parms.fibers,atlas_flag,0,0,...
          0,[],...
          0,[],1,...
          parms.count_flag);
      if ~some_exist_flag
        fprintf('%s: WARNING: no %s fiber files with fiber infix %s found in %s\n',...
          fiber_ext,fiber_infix,fiber_dir);
        parms.fiber_flag = 0;
      end;
    end;
  end;

  % check fiber dispersion weighting parameters
  if parms.fiber_disp_flag
    if length(parms.fiber_dispvec) ~= 1 &&...
       length(parms.fiber_dispvec) ~= length(parms.fibers)
      error('mismatch between fiber_dispvec and fibers');
    end;
    parms.fname_MD = sprintf('%s_MD.mgz',parms.full_DT_fstem);
    if ~exist(parms.fname_MD,'file')
      fprintf('%s: WARNING: MD file %s not found\n',...
        mfilename,parms.fname_MD);
      parms.fiber_disp_flag = 0;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [all_exist_flag,some_exist_flag,fiber_infix,fiber_ext] = ...
 check_fibers(fiber_dir,fibers,atlas_flag,thresh_FA,thresh_prob,xcg_flag,...
    xcg_suffix,disp_flag,disp_suffix,dispfact,count_flag)

  all_exist_flag = 1;
  some_exist_flag = 1;

  [fiber_infix,fiber_ext] = dti_set_fiber_infix('atlas_flag',atlas_flag,...
    'thresh_FA',thresh_FA,'thresh_prob',thresh_prob,...
    'xcg_flag',xcg_flag,'xcg_suffix',xcg_suffix,...
    'disp_flag',disp_flag,'disp_suffix',disp_suffix,'dispfact',dispfact,...
    'count_flag',count_flag);
  flist = dir(sprintf('%s/fiber_*_%s%s',fiber_dir,fiber_infix,fiber_ext));
  if isempty(flist)
    all_exist_flag = 0;
    some_exist_flag = 0;
    return;
  end;
  for f=1:length(fibers)
    fname = sprintf('%s/fiber_%03d_%s%s',...
      fiber_dir,fibers(f),fiber_infix,fiber_ext);
    if ~exist(fname,'file')
      all_exist_flag = 0;
      return;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_aseg(parms)
  if parms.verbose
    fprintf('%s: checking aseg...\n',mfilename);
  end;
  % set name of aseg
  if parms.aseg_aparc_flag
    parms.aseg_name = 'aparc+aseg';
  else
    parms.aseg_name = 'aseg';
  end;
  parms.fname_aseg = sprintf('%s/mri/%s.mgz',parms.fspath,parms.aseg_name);
  if ~exist(parms.fname_aseg,'file')
    fprintf('%s: WARNING: aseg file %s not found\n',...
      mfilename,parms.fname_aseg);
    parms.aseg_flag = 0;
    parms.brainmask_flag = 0;
    if parms.xcg_flag, parms.fiber_flag = 0; end;
  end

  % check aseg dispersion weighting parameters
  if parms.aseg_flag && parms.aseg_disp_flag
    if isempty(parms.aseg_dispvec)
      error('aseg_dispvec required if aseg_disp_flag = 1');
    end;
    if isempty(parms.aseg_disp_roicodes)
      error('aseg_disp_roicodes required if aseg_disp_flag = 1');
    end;
    if length(parms.aseg_dispvec) ~= length(parms.aseg_disp_roicodes)
      error('mismatch between aseg_dispvec and aseg_disp_roicodes');
    end;
    parms.fname_MD = sprintf('%s_MD.mgz',parms.full_DT_fstem);
    if ~exist(parms.fname_MD,'file')
      fprintf('%s: WARNING: MD file %s not found\n',...
        mfilename,parms.fname_MD);
      parms.aseg_disp_flag = 0;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_wmparc(parms)
  if parms.verbose
    fprintf('%s: checking wmparc...\n',mfilename);
  end;
  % set name of wmparc
  parms.wmparc_name = 'wmparc';
  parms.fname_wmparc = sprintf('%s/mri/%s.mgz',parms.fspath,parms.wmparc_name);
  if ~exist(parms.fname_wmparc,'file')
    fprintf('%s: WARNING: wmparc file %s not found\n',...
      mfilename,parms.fname_wmparc);
    parms.wmparc_flag = 0;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_cortsurf(parms)
  if parms.verbose
    fprintf('%s: checking cortical surface...\n',mfilename);
  end;
  % check aparc
  [parms,exist_flag] = check_aparc(parms);
  parms.cortsurf_flag = exist_flag;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_gwcsurf(parms)
  if parms.verbose
    fprintf('%s: checking cortical surface...\n',mfilename);
  end;
  % check cortical ribbon file
  if parms.gwmask_flag
    fname_ribbon = sprintf('%s/mri/ribbon.mgz',parms.fspath);
    if ~exist(fname_ribbon,'file')
      fprintf('%s: WARNING: ribbon file %s not found\n',fname_ribbon);
      parms.gwcsurf_flag = 0;
    end;
  end;
  % check aparc
  [parms,exist_flag] = check_aparc(parms);
  parms.gwcsurf_flag = exist_flag;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,exist_flag] = check_aparc(parms)
  exist_flag = 1;
  % check fnames_aparc
  parms.nhemi = length(parms.hemilist);
  if isempty(parms.fnames_aparc)
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      parms.fnames_aparc{h} = sprintf('%s/label/%s.aparc.annot',...
        parms.fspath,hemi);
    end;
  else
    if ~iscell(parms.fnames_aparc)
      parms.fnames_aparc = {parms.fnames_aparc};
    end;
    if length(parms.fnames_aparc) ~= parms.nhemi
      error('must have %d elements in fnames_aparc (have %d)',...
        parms.nhemi,length(parms.fnames_aparc));
    end;
  end;
  for h=1:parms.nhemi
    if mmil_isrelative(parms.fnames_aparc{h})
      parms.fnames_aparc{h} = [parms.fspath '/label/' parms.fnames_aparc{h}];
    end;
    if ~exist(parms.fnames_aparc{h},'file')
      fprintf('%s: WARNING: aparc file %s not found\n',parms.fnames_aparc{h});
      exist_flag = 0;
      break;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = load_RegInfo(parms)
  errcode = 0;
  % load registration info
  if parms.verbose
    fprintf('%s: loading registration info from %s...\n',mfilename,parms.cpath);
  end;
  args = mmil_parms2args(parms,parms.reg_tags);
  [RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(parms.cpath,args{:});
  if errcode
    fprintf('%s: ERROR: *_regT1.mat file not found\n',mfilename);
    errcode = 1; return;
  end;
  parms.M_T1_to_DTI = RegInfo.M_T1_to_T2;
  parms.volsz_DTI = RegInfo.volsz_T2(1:3);
  parms.M_DTI = RegInfo.M_T2;
  parms.volsz_T1 = RegInfo.volsz_T1(1:3);
  parms.M_T1 = RegInfo.M_T1;
  parms.fname_T1 = RegInfo.fname_T1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_brainmask(parms)
  % create brain mask from aseg, resample to DTI resolution
  parms.res_outfix = parms.resDTI_outfix;
  args = mmil_parms2args(parms,parms.mask_tags);
  dti_create_asegmask(parms.fspath,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = transform_fibers(parms)
  orig_parms = parms;
  if parms.resT1flag
    fiber_outdir = [fiber_outdir '_resT1'];
    paths_outdir = [paths_outdir '_resT1'];
    parms.fname_fseg = [parms.fseg_fstem '_resT1'];
    parms.create_paths_flag = 0; % not currently supported
  elseif parms.fseg_resT1flag
    parms.fname_fseg = [parms.fseg_fstem '_resT1'];
  else
    parms.fname_fseg = [parms.fseg_fstem '_resDTI'];
  end;
  if parms.xcg_flag
    parms.fname_fseg = [parms.fname_fseg '_' parms.xcg_suffix];
  end;
  if parms.fiber_disp_flag
    parms.disp_flag = 1;
    parms.fname_vals = parms.fname_MD;
    parms.dispvec = parms.fiber_dispvec;
    parms.dispfact = parms.fiber_dispfact;
    parms.disp_xcg_flag = parms.fiber_disp_xcg_flag;
    parms.disp_suffix = parms.fiber_disp_suffix;
    parms.fname_fseg = sprintf('%s_%s%0.1f',...
      parms.fname_fseg,parms.disp_suffix,parms.dispfact);
  end;
  parms.fname_fseg = [parms.fname_fseg '.mgz'];
  % loop over multiple atlas flag, thresh_prob, and thresh_FA
  atlas_flag_list = parms.atlas_flag;
  thresh_prob_list = parms.thresh_prob;
  thresh_FA_list = parms.thresh_FA;
  for a=1:length(atlas_flag_list)
    parms.atlas_flag = atlas_flag_list(a);
    if parms.atlas_flag
      if parms.fiber_ATdir_flag, parms.outdir = parms.atlastrack_dir; end;
      parms.fiber_dir = [parms.atlastrack_dir '/' parms.mapsdir];
      parms.fiber_outdir = [parms.outdir '/' parms.mapsdir];
      parms.paths_outdir = [parms.outdir '/' parms.pathsdir];
      if ~isempty(parms.fiber_atlasname)
        parms.fiber_dir = [parms.fiber_dir '_' parms.fiber_atlasname];
        parms.fiber_outdir = [parms.fiber_outdir '/' parms.fiber_atlasname];
        parms.pathsdir = [parms.pathsdir '/' parms.fiber_atlasname];
      end;
    else
      if parms.fiber_ATdir_flag, parms.outdir = parms.cpath; end;
      parms.fiber_dir = [parms.cpath '/' parms.manual_mapsdir];
      parms.fiber_outdir = [parms.outdir '/' parms.manual_mapsdir];
      parms.paths_outdir = [parms.outdir '/' parms.manual_pathsdir];
    end;
    for p=1:length(thresh_prob_list)
      parms.thresh_prob = thresh_prob_list(p);
      for f=1:length(thresh_FA_list)
        parms.thresh_FA = thresh_FA_list(f);
        all_exist_flag = ...
          check_fibers(parms.fiber_dir,parms.fibers,parms.atlas_flag,...
            parms.thresh_FA,parms.thresh_prob,...
            parms.xcg_flag,parms.xcg_suffix,...
            parms.fiber_disp_flag,parms.fiber_disp_suffix,parms.fiber_dispfact,...
            parms.count_flag);
        if ~all_exist_flag
          args = mmil_parms2args(parms,parms.transform_tags);
          try
            dti_transform_fibers(parms.fiber_dir,args{:})
          catch
            fprintf('\n%s: WARNING: dti_transform_fibers failed for %s:\n%s\n\n',...
              mfilename,parms.cpath,lasterr);
            parms.fiber_flag = 0;
          end;
        end;
      end;
    end;
  end;
  orig_parms.fiber_flag = parms.fiber_flag;
  parms = orig_parms;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = prep_data(parms)
  errcode = 0;
  valid_flag = ones(1,parms.ninputs);
  parms.fnamelist_resT1 = cell(1,parms.ninputs);
  parms.fnamelist_resDTI = cell(1,parms.ninputs);
  parms.resT1_flags = zeros(1,parms.ninputs);
  for i=1:parms.ninputs
    fname = parms.fnamelist{i};
    meas = parms.measlist{i};
    full_meas = parms.full_measlist{i};
    fnum = parms.fnumlist(i);
    % determine whether input has DTI or T1 resolution
    [volsz,M] = read_volsz(fname,parms);
    if all(volsz==parms.volsz_DTI) && all(M(:)==parms.M_DTI(:))
      resDTI_flag = 1;
    else
      resDTI_flag = 0;
    end;
    if all(volsz==parms.volsz_T1) && all(M(:)==parms.M_T1(:))
      resT1_flag = 1;
    else
      resT1_flag = 0;
    end;
    if ~resDTI_flag && ~resT1_flag
      fprintf('%s: WARNING: input file %s has neither DTI nor T1 resolution',...
        mfilename,fname);
      valid_flag(i) = 0;
      continue;
    end;
    parms.resT1_flags(i) = resT1_flag;
    % normalize b0 image
    switch meas
      case 'T2w'
        fname = norm_b0_MD(fname,parms);
        parms.fnamelist{i} = fname;
      case 'b0N'
        fname = norm_b0_aseg(fname,parms);
        parms.fnamelist{i} = fname;
    end;
    % resample to T1 or DTI resolution as necessary
    if resT1_flag
      parms.fnamelist_resT1{i} = fname;
      if ~parms.resT1flag
        fname = resample_data(fname,fnum,full_meas,parms,0);
        parms.fnamelist_resDTI{i} = fname;
      end;
    else
      parms.fnamelist_resDTI{i} = fname;
      if parms.resT1flag ||...
         ((parms.cortsurf_flag || parms.gwcsurf_flag) && parms.cortsurf_resT1flag) ||...
         parms.atlas_warp_flag
        fname = resample_data(fname,fnum,full_meas,parms,1);
        parms.fnamelist_resT1{i} = fname;
      end;
    end;
  end;
  ind_keep = find(valid_flag);
  parms.fnamelist = parms.fnamelist(ind_keep);
  parms.measlist = parms.measlist(ind_keep);
  parms.full_measlist = parms.full_measlist(ind_keep);
  parms.fnamelist_resT1 = parms.fnamelist_resT1(ind_keep);
  parms.fnamelist_resDTI = parms.fnamelist_resDTI(ind_keep);
  parms.ninputs = length(parms.fnamelist);
  if parms.ninputs==0
    fprintf('%s: ERROR: no valid input files found\n',mfilename);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = analyze_data(parms)
  errcode = 0;
  for i=1:parms.ninputs
    full_meas = parms.full_measlist{i};
    scalefact = parms.scalefacts(i);
    resT1_flag = parms.resT1_flags(i);
    if parms.fiber_flag
      [fname,parms.frames] = set_input(parms,i,0);
      errcode = fiber_analysis(fname,full_meas,scalefact,parms);
      if errcode, parms.fiber_flag = 0; end;
    end;
    if parms.aseg_flag
      [fname,parms.frames] = set_input(parms,i,0);
      errcode = aseg_analysis(fname,full_meas,scalefact,parms);
      if errcode, parms.aseg_flag = 0; end;
    end;
    if parms.wmparc_flag
      [fname,parms.frames] = set_input(parms,i,0);
      errcode = wmparc_analysis(fname,full_meas,scalefact,parms);
      if errcode, parms.wmparc_flag = 0; end;
    end;
    if parms.cortsurf_flag
      [fname,parms.frames] = set_input(parms,i,1);
      errcode = cortsurf_analysis(fname,full_meas,scalefact,resT1_flag,parms);
      if errcode, parms.cortsurf_flag = 0; end;
    end;
    if parms.gwcsurf_flag
      [fname,parms.frames] = set_input(parms,i,1);
      errcode = gwcsurf_analysis(fname,full_meas,scalefact,resT1_flag,parms);
      if errcode, parms.gwcsurf_flag = 0; end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname,fnum] = set_input(parms,i,cortsurf_flag)
  if ~exist('cortsurf_flag','var'), cortsurf_flag = 0; end;
  resT1_flag = parms.resT1_flags(i);
  if parms.resT1flag ||...
     (cortsurf_flag && (resT1_flag || parms.cortsurf_resT1flag))
    fname = parms.fnamelist_resT1{i};
    if resT1_flag
      fnum = parms.fnumlist(i);
    else
      % already extracted correct frame before resampling
      fnum = 1;
    end;
  else
    fname = parms.fnamelist_resDTI{i};
    if resT1_flag
      % already extracted correct frame before resampling
      fnum = 1;
    else
      fnum = parms.fnumlist(i);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = atlas_warp_data(parms)
  errcode = 0;

  % locate existing reg file in AtlasTrack dir
  if mmil_isrelative(parms.atlastrack_dir)
    parms.atlastrack_dir = [parms.cpath '/' parms.atlastrack_dir];
  end;
  fname_reg = sprintf('%s/T1_dctReg2Atlas.mat',parms.atlastrack_dir);

  if exist(fname_reg,'file')
    % use existing reg file in AtlasTrack dir
    parms.fname_reg = fname_reg;
  else
    % register T1-weighted MRI to atlas and warp T1 to atlas
    if ~exist(parms.fname_T1,'file')
      fprintf('%s: WARNING: T1 file %s not found, skipping atlas warp\n',...
        mfilename,parms.fname_T1);
      return;
      errcode = 1;
    end;
    parms.outstem_T1 = 'T1';
    args = mmil_parms2args(parms,parms.warp_tags);
    [fname_atl,parms.fname_reg] = mmil_warp_to_atlas(parms.fname_T1,args{:});
  end;

  % warp each input volume to atlas
  for i=1:parms.ninputs
    full_meas = parms.full_measlist{i};
    fname = parms.fnamelist_resT1{i};
    parms.outstem = full_meas;
    args = mmil_parms2args(parms,parms.warp_tags);
    mmil_warp_to_atlas(fname,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = resample_data(fname,fnum,full_meas,parms,resT1_flag)
  % resample data to either DTI or T1 resolution using M_T1_to_DTI
  if ~exist('resT1_flag','var'), resT1_flag = 0; end;
  if resT1_flag
    resstr = 'T1';
    nvox_ref = parms.volsz_T1;
    M_ref = parms.M_T1;
    M_reg = parms.M_T1_to_DTI;
    res_outfix = parms.resT1_outfix;
  else
    resstr = 'DTI';
    nvox_ref = parms.volsz_DTI;
    M_ref = parms.M_DTI;
    M_reg = inv(parms.M_T1_to_DTI);
    res_outfix = parms.resDTI_outfix;
  end;
  fname_out = [parms.outdir '/' full_meas '_' res_outfix '.mgz'];
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: resampling %s to %s resolution...\n',...
        mfilename,fname,resstr);
    end;
    [vol,M] = fs_load_mgh(fname,[],fnum);
    [vol,M] = mmil_resample_vol(vol,M,...
      'nvox_ref',nvox_ref,'M_ref',M_ref,'M_reg',M_reg);
    fs_save_mgh(vol,fname_out,M);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = fiber_analysis(fname,full_meas,scalefact,parms)
  errcode = 0;
  % calculate averages for fiber ROIs
  if parms.verbose
    fprintf('%s: fiber analysis for %s...\n',mfilename,fname);
  end;
  parms.outstem = [parms.outdir '/' full_meas];
  if parms.fiber_disp_flag
    parms.disp_flag = 1;
    parms.dispvec = parms.fiber_dispvec;
    parms.dispfact = parms.fiber_dispfact;
    parms.disp_suffix = parms.fiber_disp_suffix;
  end;
  parms.scalefact = scalefact;
  atlas_flag_list = parms.atlas_flag;
  thresh_prob_list = parms.thresh_prob;
  thresh_FA_list = parms.thresh_FA;
  for a=1:length(atlas_flag_list)
    parms.atlas_flag = atlas_flag_list(a);
    if parms.atlas_flag
      orig_fiber_dir = [parms.atlastrack_dir '/' parms.mapsdir];
      if parms.fiber_ATdir_flag
        fiber_dir = orig_fiber_dir;
      else
        fiber_dir = [parms.outdir '/' parms.mapsdir];
      end;
      if ~isempty(parms.fiber_atlasname)
        orig_fiber_dir = [parms.orig_fiber_dir '_' parms.fiber_atlasname];
        fiber_dir = [parms.fiber_dir '_' parms.fiber_atlasname];
      end;
    else
      orig_fiber_dir = [parms.cpath '/' parms.manual_mapsdir];
      if parms.fiber_ATdir_flag
        fiber_dir = orig_fiber_dir;
      else
        fiber_dir = [parms.outdir '/' parms.manual_mapsdir];
      end;
    end;
    for p=1:length(thresh_prob_list)
      parms.thresh_prob = thresh_prob_list(p);
      for f=1:length(thresh_FA_list)
        parms.thresh_FA = thresh_FA_list(f);
        if ~parms.fiber_ATdir_flag
          all_exist_flag = ...
            check_fibers(orig_fiber_dir,parms.fibers,parms.atlas_flag,...
              parms.thresh_FA,parms.thresh_prob,...
              parms.xcg_flag,parms.xcg_suffix,...
              parms.fiber_disp_flag,parms.fiber_disp_suffix,parms.fiber_dispfact,...
              parms.count_flag);
          if all_exist_flag, fiber_dir = orig_fiber_dir; end;
        end;
        args = mmil_parms2args(parms,parms.fiber_tags);
        try
          mmil_fiber_analysis(fname,fiber_dir,args{:});
        catch
          fprintf('\n%s: WARNING: fiber analysis failed:\n%s\n\n',...
            mfilename,lasterr);
          errcode = 1;
          return;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = aseg_analysis(fname,full_meas,scalefact,parms)
  errcode = 0;
  % calculate averages for aseg ROIs
  if parms.verbose
    fprintf('%s: aseg analysis for %s...\n',mfilename,fname);
  end;
  parms.outstem = [parms.outdir '/' full_meas];
  if parms.aseg_disp_flag
    parms.fname_vals = parms.fname_MD;
    parms.dispvec = parms.aseg_dispvec;
    parms.disp_roicodes = parms.aseg_disp_roicodes;
    parms.disp_suffix = parms.aseg_disp_suffix;
    parms.dispfact = parms.aseg_dispfact;
  end;
  parms.scalefact = scalefact;
  parms.res_outfix = parms.resDTI_outfix;
  parms.M_reg = parms.M_T1_to_DTI;
  args = mmil_parms2args(parms,parms.aseg_tags);
  try
    mmil_aseg_analysis(fname,parms.fspath,args{:});
  catch
    fprintf('\n%s: WARNING: aseg analysis failed:\n%s\n\n',...
      mfilename,lasterr);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = wmparc_analysis(fname,full_meas,scalefact,parms)
  errcode = 0;
  % calculate averages for wmparc ROIs
  if parms.verbose
    fprintf('%s: wmparc analysis for %s...\n',mfilename,fname);
  end;
  parms.outstem = [parms.outdir '/' full_meas];
  parms.scalefact = scalefact;
  parms.res_outfix = parms.resDTI_outfix;
  parms.M_reg = parms.M_T1_to_DTI;
  args = mmil_parms2args(parms,parms.wmparc_code_tags);
  parms.aseg_roilist = mmil_wmparc_roicodes(args{:});
  parms.fname_aseg = parms.fname_wmparc;
  parms.aseg_aparc_flag = 0;
  args = mmil_parms2args(parms,parms.aseg_tags);
  try
    mmil_aseg_analysis(fname,parms.fspath,args{:});
  catch
    fprintf('\n%s: WARNING: aseg analysis failed:\n%s\n\n',...
      mfilename,lasterr);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = cortsurf_analysis(fname,full_meas,scalefact,resT1_flag,parms)
  errcode = 0;
  % calculate averages for cortsurf ROIs
  if parms.verbose
    fprintf('%s: cortsurf analysis for %s...\n',mfilename,fname);
  end;
  parms.outstem = [parms.outdir '/' full_meas];
  parms.outfix = 'cortsurf';
  parms.scalefact = scalefact;
  if resT1_flag || parms.resT1flag || parms.cortsurf_resT1flag
    parms.M_reg = [];
    parms.resT1flag = 0;
  else
    parms.M_reg = parms.M_T1_to_DTI;
  end;
  args = mmil_parms2args(parms,parms.cortsurf_tags);
%  try
    mmil_cortsurf_analysis(fname,parms.fspath,args{:});
%  catch
%    fprintf('\n%s: WARNING: cortsurf analysis failed:\n%s\n\n',...
%      mfilename,lasterr);
%    errcode = 1;
%  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = gwcsurf_analysis(fname,full_meas,scalefact,resT1_flag,parms)
  errcode = 0;
  % calculate averages for gwcsurf ROIs
  if parms.verbose
    fprintf('%s: gwcsurf analysis for %s...\n',mfilename,fname);
  end;
  parms.outstem = [parms.outdir '/' full_meas];
  if ~isempty(parms.gwcsurf_outfix)
    parms.outstem = [parms.outstem '_' parms.gwcsurf_outfix];
    parms.wm_outstem = ['wmmask_' parms.gwcsurf_outfix];
    parms.gm_outstem = ['gmmask_' parms.gwcsurf_outfix];
  end;
  parms.scalefact = scalefact;
  if resT1_flag || parms.resT1flag || parms.cortsurf_resT1flag
    parms.M_reg = [];
    parms.resT1flag = 0;
  else
    parms.M_reg = parms.M_T1_to_DTI;
  end;
  args = mmil_parms2args(parms,parms.gwcsurf_tags);
%  try
    mmil_gwcsurf_analysis(fname,parms.fspath,args{:});
%  catch
%    fprintf('\n%s: WARNING: gwcsurf analysis failed:\n%s\n\n',...
%      mfilename,lasterr);
%    errcode = 1;
%  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = norm_b0_MD(fname,parms)
  % normalize b=0 image to MD
  if parms.b0_norm_DTcalc_flag
    fname_out = sprintf('%s_T2w.mgz',parms.full_DT_fstem);
    fname_sf = sprintf('%s_T2w_sf.mat',parms.full_DT_fstem);
  elseif parms.full_fstem_flag
    fname_out = sprintf('%s/%s_T2w.mgz',parms.outdir,parms.DT_fstem);
    fname_sf = sprintf('%s/%s_T2w_sf.mat',parms.outdir,parms.DT_fstem);
  elseif ~isempty(parms.DT_outfix)
    fname_out = sprintf('%s/%s_T2w.mgz',parms.outdir,parms.DT_outfix);
    fname_sf = sprintf('%s/%s_T2w_sf.mat',parms.outdir,parms.DT_outfix);
  else
    fname_out = sprintf('%s/T2w.mgz',parms.outdir);
    fname_sf = sprintf('%s/T2w_sf.mat',parms.outdir);
  end;
  fname_MD = sprintf('%s_MD.mgz',parms.full_DT_fstem);
  mmil_normvol_linfit(fname,fname_out,fname_MD,...
    'fname_sf',fname_sf,...
    'scalefact',parms.T2w_scalefact,...
    'forceflag',parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = norm_b0_aseg(fname,parms)
  % resample aseg to DTI space
  aseg_name = [parms.aseg_name '_' parms.resDTI_outfix];
  fname_aseg_res = [parms.outdir '/' aseg_name '.mgz'];
  if ~exist(fname_aseg_res,'file') || parms.forceflag
    [vol,M] = fs_load_mgh(parms.fname_aseg);
    [vol,M] = mmil_resample_vol(vol,M,...
      'nvox_ref',parms.volsz_DTI,'M_ref',parms.M_DTI,...
      'interpm',0,'M_reg',inv(parms.M_T1_to_DTI));
    fs_save_mgh(vol,fname_aseg_res,M);
  end;
  % erode aseg
  if parms.b0N_erode_flag
    aseg_name = [aseg_name '_erode'];
    if parms.erode_nvoxels>1
      aseg_name = sprintf('%s%d',aseg_name,parms.erode_nvoxels);
    end;
    fname_aseg_erode = [parms.outdir '/' aseg_name '.mgz'];  
    fs_erode_aseg(fname_aseg_res,fname_aseg_erode,...
      'nvoxels',parms.erode_nvoxels);
    fname_aseg_res = fname_aseg_erode;
  end;
  % normalize b=0 image to ventricles
  if parms.b0_norm_DTcalc_flag
    fname_out = sprintf('%s_b0N.mgz',parms.full_DT_fstem);
  elseif parms.full_fstem_flag
    fname_out = sprintf('%s/%s_b0N.mgz',parms.outdir,parms.DT_fstem);
  else
    fname_out = sprintf('%s/b0N.mgz',parms.outdir);
  end;
  mmil_normvol_asegroi(fname,fname_out,fname_aseg_res,...
    parms.b0N_roilist,parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [volsz,M] = read_volsz(fname,parms)
  volsz = []; M = [];
  [M,volsz] = mmil_load_mgh_info(fname,parms.forceflag,parms.outdir);
  volsz = volsz(1:3);
return;

