function MMIL_Process_Exams(ProjID,varargin)
% function MMIL_Process_Exams(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%       (e.g. read from VisitInfo.csv file with MMIL_Read_StudyInfo)
%    Use this option to limit processing to specific subjects
%    Should have a field 'VisitID' that indicates the name of the
%       input data directory (i.e. in orig data directory)
%    If VisitID field is missing, will use SubjID field instead
%    May also specify subject-specific parameter values
%     but note that these will override ProjInfo and command line input
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: raw, proc
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'procstep': [0|1|2] at which stage to start processing
%    0 = start from orig, 1 = start from raw, 2 = start from proc
%    {default = 0}
%  'newflag': [0|1] whether to create jobs only for new data
%    ignored if procstep=2 or forceflag=1
%    {default = 0}
%  'unpackflag': [0|1] unpack dicoms from orig to raw only
%    {default = 0}
%  'preprocflag': [0|1] convert raw dicoms to mgz and stop
%    {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%  'batchname': name of output batchdir
%     {default = 'MMIL_Process_Exams'}
%
% Optional Dicom Parameters
%  'DCM_linkflag': [0|1] whether to create links instead of copying from
%     orig to raw
%     {default = 1}
%  'DCM_classify_set': string specifying default classify file
%     in MMPS/parms/classify e.g. MMIL_Series_Classify_{set}.csv
%     'All', 'GE', 'Philips', 'Siemens', or 'Strict'
%     Note: 'Strict' is same as 'All', excluding those that
%       rely on SeriesDescription
%     May also use 'none', which will not use classification rules
%       instead rely on function mmil_classify_by_code which
%       does nothing unless overridden with custom-modified version of it
%     {default = 'All'}
%  'DCM_classify_file': full path name of csv file with rules
%     for classifying series of dicoms
%     if empty, and file called {ProjID}_Series_Classify.csv is
%       found in ProjInfo/{ProjID}, will use that
%     if specified, or file exists in ProjInfo/{ProjID}
%       will override DCM_classify_set
%     {default = []}
%
% Optional Structural Parameters
%  'STRUCTflag': [0|1] whether to process structural data
%    {default = 1}
%  'STRUCT_T1type': which type of T1 series ('MPR' or 'hiFA') to use
%     as reference (i.e. which one gets registered to atlas and 
%     to which other scans are registered)
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%  'STRUCT_gradunwarp_flag': [0|1|2] whether to correct for gradient warping
%     if 2, processing is aborted if gradwarp info is missing
%     {default = 1}
%  'STRUCT_wmbc_flag': [0|1] whether to correct for intensity bias-field
%     using white matter segmentation and sparse smoothing
%    {default = 0}
%  'STRUCT_nu_flag': [0|1] whether to correct for intensity bias-field 
%     using nu (N3) correction
%     {default = 0}
%  'STRUCT_tal_flag': [0|1] register to Talairach space
%    for standardization of white matter intensity values
%    during nu correction
%    {default = 0}
%  'STRUCT_atlasflag': [0|1] whether to resample output in atlas space
%              (rigid-body registration only)
%    {default = 1}
%  'STRUCT_nativeflag': [0|1] if STRUCT_atlasflag=0
%     whether to keep native resolution
%     otherwise, resample to 1mm, 256^3, LIA
%     {default = 0}
%  'STRUCT_rawQCflag': [0|1] whether to use manual raw QC info
%    {default = 0}
%  'STRUCT_minmax': vector of minimum and maximum number of scans required
%    {default = [1 Inf]}
%  'STRUCT_scantypes': cell array of scan types to register
%     {default = {'MPR','XetaT2','FLASHhi','FLASHlo','MEDIChi','MEDIClo'}}
%  'STRUCT_BEMflag': [0|1] whether to generate BEM surfaces
%    {default = 0}
%  'STRUCT_first_flag': [0|1] whether to run FSL's FIRST
%     to generate hippocampal volumes
%    {default = 0}
%  'STRUCT_first_nu_flag': [0|1] whether to use FreeSurfer nu.mgz
%     when running FSL's FIRST (requires FreeSurfer recon already exists)
%    {default = 0}
%  'STRUCT_oldreg_flag': [0|1] use old registration mriRegister
%     otherwise use mmil_reg in mmil_average_volumes
%     {default = 0}
%  'STRUCT_prereg_flag': [0|1] whether to use reg for rigid registration
%     before running dct registration to atlas
%     {default = 1}  
%  'STRUCT_atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%  'STRUCT_atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%
% Optional DTI Parameters
%  'DTIflag': [0|1] whether to process DTI data
%    {default = 0}
%  'DTI_snums': list of DTI scan numbers to use
%     if empty, use all DTI scans in container
%     {default = []}
%  'DTI_snums_flag': [0|1|2|3] which set of "snums" to use
%     0: use all available DTI scans in container
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%     {default = 1}
%  'DTI_B0unwarp_flag': [0|1|2] whether to correct for B0 distortions
%    if 2, processing aborted if B0 unwarping fails (e.g. no rev scans)
%    {default = 1}
%  'DTI_motion_B0uw_flag': [0|1] whether to estimate head motion
%     and apply to B0dx field
%     {default = 0}
%  'DTI_ecc_flag': [0|1] whether to perform eddy current correction
%     {default=1}
%  'DTI_censor_flag': [0|1] reject bad slices based on tensor fit
%    {default=1}
%  'DTI_mc_flag': [0|1] whether to perform motion correction
%     {default=1}
%  'DTI_gradunwarp_flag': [0|1|2] whether to correct for gradient warping
%     if 2, processing is aborted if gradwarp info is missing
%    {default = 1}
%  'DTI_resample_flag': [0|1] whether to resample processed diffusion data
%     e.g. to make isotropic, to correct inter-scan motion, to register to T1
%    {default = 1}
%  'DTI_regT1flag': [0,1,2] whether to register DTI data to T1
%    0=do not register, 1=register only, 2=apply registration to diffusion data
%    {default = 1}
%  'DTI_export_flag': [0|1|2] export processed diffusion data 
%      0: do not export
%      1: in export_type format
%      2: with customized steps for FSL
%      {default = 0}
%  'DTI_CSD_tracto_flag': [0|1] run CSD tractography on processed diffusion data
%    {default = 0}
%  'DTI_calcDT_flag': [0|1] whether to perform diffusion tensor calculations
%    {default = 1}
%  'DTI_calcRSI_flag': [0|1] whether to perform Restriction Spectrum Imaging
%    calculations on processed diffusion data
%    {default = 0}
%  'DTI_ATLflag': [0|1] whether to do DTI fiber atlas tracking
%    {default = 0}
%  'DTI_ATL_regFA_flag': [0|1] use FA image for dct morph for AtlasTrack
%    {default = 0}
%  'DTI_revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use either non-rev and/or rev data
%     {default = 0}
%  'DTI_min_ndirs': minimum number of gradient directions allowed
%    for tensor calculations
%    {default = 6}
%  'DTI_min_bval': minimum b-value allowed for tensor calculations
%    {default = 0}
%  'DTI_max_bval': maximum b-value used for tensor calculations
%    {default = Inf}
%  'DTI_flex_flag': [0|1] DTI_flex scans included in tensor fit
%    {default = 0}
%  'DTI_min_nb0': minimum number of b=0 images required for tensor calculations
%    Note: some scans (e.g. 1-dir diffusion) do have a b=0 image but they
%      report 0 b=0 images in the dicom header
%    {default = 1}
%  'DTI_censor_min_ndirs': minimum number of diffusion directions (not including
%    b=0 images) required for tensor fit after censoring
%    will not do censoring if remaining ndirs < min_ndirs
%    {default = 12}
%  'DTI_censor_thresh': error threshold for censoring bad frames
%    normalized to median error for each slice; higher values mean less censoring
%    {default=3.2}
%  'DTI_optimize_B0uw_flag': [0|1] search for optimal B0 unwarp parameters
%     kernelWidthMax and lambda2
%    {default = 0}
%  'DTI_driftcorr': [0|1] estimate drift correction with eddy current correction
%     {default = 0}
%  'DTI_motion_B0uw_iters': number of iterations to estimate motion
%     and B0 displacement
%    {default = 2}
%  'DTI_min_trans': minimum translation (mm) for estimating motion
%    {default = 0.05}
%  'DTI_min_rot': minimum rotation (degrees) for estimating motion
%    {default = 0.05}
%  'DTI_nvoxels': vector of number of resampled voxels [nx,ny,nz]
%    {default = [120 120 70]}
%  'DTI_resolution': desired voxel sizes of resampled data (x y z)
%    {default = [2 2 2]}
%  'DTI_std_orient' : specify the resampled output orientation 
%    {default = []}  
%  'DTI_smooth': 3D isotropic full width half max blurring kernel (mm) applied
%     to processed diffusion data before tensor calculations
%     {default = 0}
%  'DTI_rot': rotation applied to processed diffusion data (x,y,z deg)
%     {default = [0,0,0]}
%  'DTI_trans': translation applied to processed diffusion data (x,y,z mm)
%     {default = [0,0,0]}
%  'DTI_bbregflag': [0|1] whether to register DTI to T1 using
%     FreeSurfer's bbregister
%    {default = 0}
%  'DTI_nii_out_orient': output orientation for nii format
%    if empty, keep native orientation
%    {default = []}
%  'DTI_nob0_flag': [0|1] exclude b=0 images from tensor fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'DTI_mask_DTmeas_flag': whether or not to mask DTmeas output volumes
%     {default = 1}
%  'DTI_DT_regT1flag': [0|1|2] whether to resample DT output
%    0=not resampled, 1=registered T1, 2=resampled T1
%    {default = 0}
%  'DTI_DT_nonlin_flag': [0|1] use nonlinear optimization for tensor fit
%     with initial parameters from linear fit
%     {default = 0}
%  'DTI_mask_RSImeas_flag': whether or not to mask RSImeas output volumes
%     {default = 1}
%  'DTI_RSI_regT1flag': [0|1|2] whether to resample RSI output
%    0=not resampled, 1=registered T1, 2=resampled T1
%    {default = 0}
%  'DTI_datflag': [0|1] whether to convert DTI FA and V0 to dat format
%    {default = 0}
%  'DTI_niiflag': [0|1] whether to convert DTI measures to nii format
%    {default = 0}
%  'DTI_outfix': string added to DTI file name after correction
%    {default = 'corr'}
%  'DTI_deoblique_flag': [0|1] whether to resample oblique slices to on-axis
%     ignored if native_flag = 1
%     {default = 1}
%  'DTI_fibers': fiber numbers to generate
%    ignored if DTI_ATLflag = 0
%     {default = [101:110,115:123,133:138,141:150]}
%  'DTI_subdiv_fibers': vector of fiber subdivision numbers
%    ignored if DTI_ATLflag = 0
%     {default = [1014,1024,1231,1232]}
%  'DTI_divide_fibers_flag': [0|1] divide fibers into subdivisions
%    ignored if DTI_ATLflag = 0
%     {default = 1}
%  'DTI_combine_fibers_flag': [0|1] combine fibers (e.g. left and right hemi)
%    ignored if DTI_ATLflag = 0
%     {default = 1}
%  'DTI_locflag': [0|1] use location information only for AtlasTrack
%    ignored if DTI_ATLflag = 0
%    {default = 0}
%  'DTI_xcg_flag': [0|1] exclude CSF and gray-mattter from fiber ROIs
%    ignored if DTI_ATLflag = 0
%    {default = 0}
%  'DTI_xcg_suffix': file name string included for CSF/gray excluded
%     ignored if DTI_ATLflag = 0 or DTI_xcg_flag = 0
%     {default = 'xcg'}
%  'DTI_xcg_codes': FreeSurfer aseg codes used to define xcg mask
%     ignored if DTI_ATLflag = 0 or DTI_xcg_flag = 0
%     {default = [0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63]}
%  'DTI_masksf_flag': [0|1] voxels with multiple fibers excluded
%    ignored if DTI_ATLflag = 0
%    {default = 0}
%  'DTI_masksf_suffix': file name string included for multiple fibers excluded
%     ignored if DTI_ATLflag = 0 or DTI_masksf_flag = 0
%     {default = 'masksf'}
%  'DTI_fseg_flag': [0|1] create fiber segmentation volume
%    ignored if DTI_ATLflag = 0
%    {default = 1}
%  'DTI_fiber_resT1flag': [0|1] resample AtlasTrack fibers to T1 resolution
%    {default=0}
%  'DTI_fiber_atlasdir': full path containing DTI fiber atlas files
%    ignored if DTI_ATLflag = 0
%    {default = [getenv('MMPS_DIR') '/atlases/DTI_Atlas/AllSubjects']}
%  'DTI_fiber_atlasname': name of the atlas fibers to be used
%    is not empty, output directory names will include this
%    ignored if DTI_ATLflag = 0
%    {default=[]} 
%
% Optional BOLD Parameters
%  'BOLDflag': [0|1] whether to process BOLD data
%    {default = 0}
%  'BOLD_snums': list of BOLD scan numbers to use
%     if empty, use all BOLD scans in container
%     {default = []}
%  'BOLD_outfix': string attached to BOLD file names after correction
%    {default = 'corr'}
%  'BOLD_B0unwarp_flag': [0|1] whether to correct for B0 distortions
%    if 2, processing aborted if B0 unwarping fails (e.g. no rev scans)
%    {default = 1}
%  'BOLD_optimize_B0uw_flag': [0|1] search for optimal B0 unwarp parameters
%     kernelWidthMax and lambda2
%    {default = 0}
%  'BOLD_regref_B0uw_flag': [0|1] whether to register reference image (fname_ref)
%     to fname_for (if revflag = 0) or fname_rev (if revflag = 1)
%     and apply transformation to B0dx
%    {default = 1}
%  'BOLD_motion_B0uw_flag': [0|1] whether to estimate head motion
%    and apply to B0dx field
%    {default = 0}
%  'BOLD_motion_B0uw_iters': number of iterations to estimate motion
%     and B0 displacement
%     {default = 2}
%  'BOLD_tshift_flag': [0|1] whether to correct for differences in
%    slice timing using AFNI's 3dTshift
%    {default = 1}
%  'BOLD_mc_flag': [0|1] whether to correct for motion using
%     AFNI's 3dvolreg
%     {default = 1}
%  'BOLD_gradunwarp_flag': [0|1|2] whether to correct for gradient distortions
%     if 2, processing is aborted if gradwarp info is missing
%    {default = 1}
%  'BOLD_regT1flag': [0,1,2] whether to register BOLD data to T1
%     0=do not register, 1=register only, 2=apply registration to BOLD data
%     {default = 1}
%  'BOLD_bbregflag': [0|1] whether to register BOLD to T1 using
%     FreeSurfer's bbregister
%    {default = 0}
%  'BOLD_tpattern': slice time pattern
%    Allowed values: {'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z'}
%    {default = 'alt+z'}
%  'BOLD_skipTRs': number of TRs at beginning of scan to be ignored in time shifting
%    {default = 0}
%  'BOLD_resample_flag': [0|1] whether to resample processed BOLD data
%     e.g. to make isotropic, to correct inter-scan motion, to register to T1
%     {default = 1}
%  'BOLD_native_flag': resample but keep original nvoxels and resolution
%     {default = 0}
%  'BOLD_nvoxels': vector of number of resampled voxels [nx,ny,nz]
%    ignored if BOLD_native_flag = 1
%    {default = [120 120 70]}
%  'BOLD_resolution': desired voxel sizes of resampled data (x y z)
%    ignored if BOLD_native_flag = 1
%    {default = [2 2 2]}
%  'BOLD_deoblique_flag': [0|1] whether to resample oblique slices to on-axis
%    ignored if BOLD_native_flag = 1
%    {default = 1}
%  'BOLD_rot': rotation applied to processed BOLD data (x,y,z deg)
%     {default = [0,0,0]}
%  'BOLD_trans': translation applied to processed BOLD data (x,y,z mm)
%     {default = [0,0,0]}
%  'BOLD_smooth': 3D isotropic full width half max blurring kernel (mm) applied
%     to processed BOLD data
%     {default = 0}
%  'BOLD_export_flag': [0|1|2] export processed BOLD data (e.g. as nii format)
%     with value of 2, use customized steps to prepare for use with FSL
%    {default = 0}
%  'BOLD_export_smooth': blurring kernel FWHM applied to exported data
%    Only applies if export_flag = 2
%    {default = 5}
%  'BOLD_export_type': 'BRIK', 'nii', 'mgh', 'mgz'
%    {default = 'nii'}
%  'BOLD_nii_out_orient': output orientation for nii format
%     used if BOLD_export_flag = 1
%    {default = []}
%  'BOLD_extract_tseries_flag': [0|1] crate csv file containing averages
%     time series for aseg and aparc ROIs for each BOLD scan
%    {default = 0}
%
% Created:  11/17/09 by Don Hagler
% Prev Mod: 11/06/17 by Feng Xue
% Last Mod: 11/08/17 by Feng Xue
%

%% TODO: reorder options, especially DTI
%%    remove documentation for less important ones?

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'procstep',0,[0:2],...
  'newflag',false,[false true],...
  'unpackflag',false,[false true],...
  'preprocflag',false,[false true],...
  'forceflag',false,[false true],...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchname','MMIL_Process_Exams',[],...
...
  'DCM_linkflag',true,[false true],...
  'DCM_classify_set','All',{'All','GE','Philips','Siemens','Strict','none'},...
  'DCM_classify_file',[],[],...
  'DCM_RejectSeries',[],[],... % to be overridden by StudyInfo
...
  'STRUCTflag',true,[false true],...
  'STRUCT_VisitID',[],[],...
  'STRUCT_T1type',2,[0:3],...
  'STRUCT_gradunwarp_flag',1,[0:2],...
  'STRUCT_wmbc_flag',false,[false true],...
  'STRUCT_nu_flag',false,[false true],...
  'STRUCT_tal_flag',false,[false true],...
  'STRUCT_atlasflag',true,[false true],...
  'STRUCT_nativeflag',false,[false true],...
  'STRUCT_rawQCflag',false,[false true],...
  'STRUCT_minmax',[1 inf],[],...
  'STRUCT_scantypes',{'MPR','XetaT2','FLASHhi','FLASHlo','MEDIChi','MEDIClo'},[],...
  'STRUCT_BEMflag',false,[false true],...
  'STRUCT_first_flag',false,[false true],...
  'STRUCT_first_nu_flag',false,[false true],...
  'STRUCT_oldreg_flag',false,[false true],...
  'STRUCT_prereg_flag',true,[false true],...
  'STRUCT_atlasdir',[],[],...
  'STRUCT_atlasname','T1_Atlas/T1_atlas',[],... 
...
  'DTIflag',false,[false true],...
  'DTI_snums',[],[],...
  'DTI_snums_flag',1,[0:3],...
  'DTI_B0unwarp_flag',1,[0:2],...
  'DTI_motion_B0uw_flag',false,[false true],...
  'DTI_ecc_flag',true,[false true],...
  'DTI_censor_flag',true,[false,true],...
  'DTI_mc_flag',true,[false true],...
  'DTI_gradunwarp_flag',1,[0:2],...
  'DTI_resample_flag',true,[false true],...
  'DTI_regT1flag',1,[0 1 2],...
  'DTI_export_flag',0,[0 1 2],...
  'DTI_CSD_tracto_flag',false,[false true],...
  'DTI_calcDT_flag',true,[false true],...
  'DTI_calcRSI_flag',false,[false true],...
  'DTI_ATLflag',false,[false true],...
...
  'DTI_ATL_regFA_flag',false,[false true],...
  'DTI_revflag',0,[0 1 2],...
  'DTI_min_ndirs',6,[6,Inf],...
  'DTI_min_bval',1,[0,Inf],...
  'DTI_max_bval',Inf,[100,Inf],...
  'DTI_flex_flag',false,[false true],...
  'DTI_min_nb0',1,[0,Inf],...
  'DTI_censor_min_ndirs',12,[],...
  'DTI_censor_thresh',3.2,[],...
  'DTI_optimize_B0uw_flag',false,[false true],...
  'DTI_driftcorr',false,[false true],...
  'DTI_motion_B0uw_iters',2,[1:10],...
  'DTI_min_trans',0.05,[0,1],... % mm
  'DTI_min_rot',0.05,[0,1],... % degrees
  'DTI_nvoxels',[120 120 70],[10,1000],...
  'DTI_resolution',[2 2 2],[0.1,10],...
  'DTI_std_orient',[],[],...  
  'DTI_smooth',0,[],...
  'DTI_rot',[0 0 0],[],...
  'DTI_trans',[0 0 0],[],...
  'DTI_bbregflag',false,[false true],...
  'DTI_nii_out_orient',[],[],...
  'DTI_nob0_flag',false,[false true],...
  'DTI_mask_DTmeas_flag',true,[false true],...   
  'DTI_DT_regT1flag',0,[0 1 2],...
  'DTI_DT_nonlin_flag',false,[false true],...
  'DTI_mask_RSImeas_flag',true,[false true],... 
  'DTI_RSI_regT1flag',0,[0 1 2],...
  'DTI_datflag',false,[false true],...
  'DTI_niiflag',false,[false true],...
  'DTI_outfix','corr',[],...
  'DTI_deoblique_flag',true,[false true],...
  'DTI_fibers',[101:110,115:123,133:138,141:150],[],...
  'DTI_subdiv_fibers',[1014,1024,1231,1232],[],...
  'DTI_divide_fibers_flag',true,[false true],...
  'DTI_combine_fibers_flag',true,[false true],...
  'DTI_locflag',false,[false true],...
  'DTI_xcg_flag',false,[false true],...
  'DTI_xcg_suffix','xcg',[],...
  'DTI_xcg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
  'DTI_masksf_flag',false,[false true],...
  'DTI_masksf_suffix','masksf',[],...
  'DTI_fseg_flag',true,[false true],...
  'DTI_fseg_thresh_prob',[],[],...
  'DTI_fiber_resT1flag',[],[],...
  'DTI_fiber_atlasdir',[],[],...
  'DTI_fiber_atlasname',[],[],...
  'DTI_cleanupflag',true,[false true],...
...
  'DTI_CSD_seed_point_sampling',[2 2 2],[],...
  'DTI_CSD_step_size',1,[],...
  'DTI_CSD_FOD_thresh',0.1,[],...
  'DTI_CSD_angle_thresh',45,[],...
  'DTI_CSD_fiber_length_range',[50 500],[],...
  'DTI_CSD_max_order',4,[2:2:100],...
  'DTI_CSD_FA_thresh',0.8,[0 1],...
...
  'DTI_DT_outdir','DTcalc',[],...
  'DTI_DT_outfix',[],[],...
...
  'DTI_RSI_outdir','RSIcalc',[],...
  'DTI_RSI_outfix',[],[],...
  'DTI_RSI_lambda',0.1,[],...
  'DTI_RSI_iso_free_flag',true,[false true],...
  'DTI_RSI_iso_hindered_flag',true,[false true],...
  'DTI_RSI_iso_restricted_flag',true,[false true],...
  'DTI_RSI_ADC_free',3e-3,[],...
  'DTI_RSI_ADC_hindered',1.5e-3,[],...
  'DTI_RSI_ADC_long',1e-3,[],...
  'DTI_RSI_ADC_trans_min',0,[],...
  'DTI_RSI_ADC_trans_max',0.9e-3,[],...
  'DTI_RSI_num_ADC_trans',5,[],...
  'DTI_RSI_SH_order',4,[2:2:10],...
  'DTI_RSI_norm_flag',false,[false true],...
...
  'BOLDflag',false,[false true],...
  'BOLD_snums',[],[],...
  'BOLD_outfix','corr',[],...
  'BOLD_B0unwarp_flag',1,[0:2],...
  'BOLD_optimize_B0uw_flag',false,[false true],...
  'BOLD_regref_B0uw_flag',true,[false true],...
  'BOLD_motion_B0uw_flag',false,[false true],...
  'BOLD_motion_B0uw_iters',2,[1:10],...
  'BOLD_tshift_flag',true,[false true],...
  'BOLD_mc_flag',true,[false true],...
  'BOLD_gradunwarp_flag',1,[0:2],...
  'BOLD_regT1flag',1,[0:2],...
  'BOLD_bbregflag',false,[false true],...
  'BOLD_tpattern','alt+z',{'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z'},...
  'BOLD_skipTRs',0,[0,Inf],...
  'BOLD_resample_flag',true,[false true],...
  'BOLD_native_flag',false,[false true],...
  'BOLD_nvoxels',[120 120 70],[10,1000],...
  'BOLD_resolution',[2 2 2],[0.1,10],...
  'BOLD_deoblique_flag',true,[false true],...
  'BOLD_rot',[0 0 0],[],...
  'BOLD_trans',[0 0 0],[],...
  'BOLD_smooth',0,[],...
  'BOLD_export_flag',0,[0:2],...
  'BOLD_export_smooth',5,[0,100],...
  'BOLD_export_type','nii',{'nii','BRIK','mgh','mgz'},...
  'BOLD_nii_out_orient',[],[],...
  'BOLD_extract_tseries_flag',false,[false true],...
...
  'BOLD_min_trans',0.05,[0,1],... % mm
  'BOLD_min_rot',0.05,[0,1],... % degrees
  'BOLD_cleanupflag',true,[false true],...
... % B0uw optimization
  'kernelWidthMax',25,[1:100],...
  'lambda2',1100,[1:10000],...
  'kernelWidthMax_vec',[25,31,35],[1:100],...
  'lambda2_vec',[1100,1500,1900],[1:10000],...
  'multi_opt_flag',false,[false true],...
...
  'required_rootdirs',{'orig','raw','proc'},[],...
...
  'info_tags',{'VisitIDs','SubjIDs','StudyInfo','RootDirs','ignore_VisitInfo_flag',...
               'user','numvec_tags'},[],...
};
parms = mmil_args2parms(varargin,parms_filter);

% excl_tags are fields that should not be passed to MMIL_Process_Exam
excl_tags = {'DTI_snums_flag' 'DCM_classify_set'...
  'RootDirs' 'StudyInfo' 'batchname' ...
  'newflag' 'procstep' 'required_rootdirs' 'info_tags'...
};
tags = setdiff(fieldnames(parms),excl_tags);

%args = MMIL_Args(parms,'MMIL_Check_ProjID');
%[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});

args = mmil_parms2args(parms,parms.info_tags);
[ProjInfo,StudyInfo,RootDirs] = MMIL_Quick_Check_ProjID(ProjID,args{:});

if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

% set series classify file
if isempty(parms.DCM_classify_file) && ~isempty(ProjID)
  fname = sprintf('%s/ProjInfo/%s/%s_Series_Classify.csv',...
    RootDirs.home,ProjID,ProjID);
  if exist(fname,'file'), parms.DCM_classify_file = fname; end
end;
if isempty(parms.DCM_classify_file) &&...
   ~strcmp(lower(parms.DCM_classify_set),'none')
  parms.DCM_classify_file = sprintf('%s/classify/MMIL_Series_Classify_%s.csv',...
    getenv('MMPS_PARMS'),parms.DCM_classify_set);
end;
if ~isempty(parms.DCM_classify_file) && ~exist(parms.DCM_classify_file,'file')
  error('classify file %s not found',parms.DCM_classify_file);
end;

if parms.DTIflag
  % set DTIScanNums and DTI field depending on snums_flag
  StudyInfo = DTI_MMIL_Set_StudyInfo_SNums(StudyInfo,parms.DTI_snums_flag);
end;

if parms.forceflag, parms.newflag = 0; end;

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

switch parms.procstep
  case 0
    RootDir = RootDirs.orig;
    regpat = '(?<VisitID>^[^\.].+)';
  case 1
    RootDir = RootDirs.raw;
    regpat = '^MRIRAW_(?<VisitID>[\w-+^]+)_(?<StudyDate>\d{8})\..+';
  case 2
    if parms.DTIflag
      RootDir = RootDirs.proc_dti;
      regpat = '^DTIPROC_(?<VisitID>[\w-+^]+)_(?<StudyDate>\d{8})\..+';
    elseif parms.BOLDflag
      RootDir = RootDirs.proc_bold;
      regpat = '^BOLDPROC_(?<VisitID>[\w-+^]+)_(?<StudyDate>\d{8})\..+';
    else
      RootDir = RootDirs.proc;
      regpat = '^MRIPROC_(?<VisitID>[\w-+^]+)_(?<StudyDate>\d{8})\..+';
    end
end;
dirlist = dir(sprintf('%s/*',RootDir));
j = 1;
for i=3:length(dirlist) % skip . and ..
  ContainerDir = char(dirlist(i).name);
  istarball = 0;
  if isempty(regexp(ContainerDir,'\.tar$'))
    tarlist = dir(sprintf('%s/%s.tar',RootDir,ContainerDir));
    if ~isempty(tarlist), continue; end; %Let's prefer tarball here
  else
    istarball = 1;
  end

  n = regexp(ContainerDir,regpat,'names','once');
  %if isempty(n) || ~dirlist(i).isdir
  if isempty(n) || (~dirlist(i).isdir && ~istarball)
    continue;
  end;
  if parms.procstep==0
    % replace any '.' in VisitID with '_'
    %  (causes problems when inserted in job names and for ADNI)
    VisitID = regexprep(ContainerDir,'\.','_');
    StudyDate = [];
  else
    VisitID = n.VisitID;
    StudyDate = n.StudyDate;
  end;
  ContainerPath = [RootDir '/' ContainerDir];
  
  tmp_parms = parms;
  tmp_parms.DTI_snums = [];
  if ~isempty(StudyInfo)
    % check that VisitID is in StudyInfo  
    ind = find(strcmp(VisitID,{StudyInfo.VisitID}));
    if isempty(ind)
%      fprintf('%s: WARNING: VisitID %s not found in StudyInfo... skipping\n',...
%        mfilename,VisitID);
      continue;
    end;
    if length(ind)>1
      error('%s: VisitID %s is found %d times in StudyInfo - must be unique\n',...
        mfilename,VisitID,length(ind));
      return;
    end;
    tmp_parms.STRUCT_VisitID = mmil_getfield(StudyInfo(ind),'STRUCT_VisitID');
    if parms.DTIflag
      tmp_parms.DTI_snums = mmil_getfield(StudyInfo(ind),'DTIScanNums');
      if ~StudyInfo(ind).DTI, tmp_parms.DTIflag = 0; end;
    end;
    if parms.BOLDflag
      tmp_parms.BOLD_snums = mmil_getfield(StudyInfo(ind),'BOLDScanNums');
    end;
     tmp_parms.DCM_RejectSeries = mmil_getfield(StudyInfo(ind),'DCM_RejectSeries');
  end;
  
  % skip a session if output already exists
  if parms.newflag && parms.procstep<2
    skipflag = 1;
    step = 'processing';
    % check if raw container exists
    [RawContainerPath,RawContainerDir] = MMIL_Get_Container(RootDirs,VisitID,'raw');
    if isempty(RawContainerDir)
      skipflag = 0;
    else
      if parms.unpackflag
        step = 'unpacking';
      elseif parms.preprocflag
        step = 'preprocessing';
      else
        step = 'processing';
      end;
      if ~isempty(regexp(RawContainerDir,'\.tar$'))
        RawContainerDir = regexprep(RawContainerDir,'\.tar$','');
        tmp = regexp(RawContainerDir,'/','split');
        vid = tmp{end};
        cmd = sprintf('tar tf %s --occurrence --wildcards "%s/%s.txt"', RawContainerPath,vid,step);
        [status, ~] = unix(cmd);
        if ~status
          skipflag = 0;
        end;
      else
        fname_progress = sprintf('%s/%s.txt',RawContainerPath,step);
        if ~exist(fname_progress,'file')
          skipflag = 0;
        end;
      end
    end;
    if skipflag
      fprintf('%s: WARNING: %s already completed for %s... skipping\n',...
        mfilename,step,VisitID);
      continue;
    end;
  end;
  
  if parms.procstep == 2
    if istarball
      ContainerPath = regexprep(ContainerPath,'\.tar$','');
      tmp = regexp(ContainerPath,'/','split');
      vid = tmp{end};
      cmd = sprintf('tar tf %s.tar --occurrence --wildcards "%s/DTcalc"', ContainerPath,vid);
      [status, ~] = unix(cmd);
      if ~status, continue; end;
    else
      if exist(strcat([ContainerPath '/DTcalc'])), continue; end;
    end
  end

  % replace values in parms with (non-empty) values from StudyInfo
  for t=1:length(tags)
    tmp_val = mmil_getfield(StudyInfo(ind),tags{t},[]);
    if ~isempty(tmp_val), tmp_parms.(tags{t}) = tmp_val; end;
  end;

  % create script
  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = [batchdir '/' jobID '.m'];
  matfname = [batchdir '/' jobID '.mat'];
  save(matfname,'RootDirs');
  fid = fopen(jobfname,'w');
  fprintf(fid,'load(''%s'')\n',matfname);
  fprintf(fid,'MMIL_Process_Exam(''%s'',...\n',ContainerPath);
  if ~isempty(ProjID), 
     fprintf(fid,'  ''ProjID'',''%s'',...\n',ProjID);
  end
  fprintf(fid,'  ''RootDirs'',RootDirs,...\n');
  mmil_write_tags(fid,tags,tmp_parms);
  fprintf(fid,');\n');
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

% check available disk space
MMIL_Check_Usage(RootDirs.raw);
MMIL_Check_Usage(RootDirs.proc);

fprintf('%%%% Now login to a cluster and run this:\n');
if parms.unpackflag ||...
  (~parms.STRUCTflag && ~parms.DTIflag && ~parms.BOLDflag)
  fprintf('    qmatjobs %s\n',parms.batchname);
elseif parms.DTI_ATLflag
  fprintf('    qmatjobs3 %s\n',parms.batchname);
elseif parms.DTI_CSD_tracto_flag
  if min(parms.DTI_CSD_seed_point_sampling)<2
    fprintf('    qmatjobs3 %s\n',parms.batchname);
  else
    fprintf('    qmatjobs2 %s\n',parms.batchname);
  end;
else
  fprintf('    qmatjobs2 %s\n',parms.batchname);
end;

