function errcode = DTI_MMIL_AtlasTrack_Exam(ContainerPath,varargin)
%function errcode = DTI_MMIL_AtlasTrack_Exam(ContainerPath,[options])
%
% Required Input:
%   ContainerPath: full path of directory containing processed diffusion data
%
% Optional Parameters:
%   'fname_FA': full path of FA volume (fractional anisotropy)
%     if specified, will use this instead of FA in ContainerPath
%     {default = []}
%   'fname_V0': full path of V0 volume (principal Eigenvector)
%     if specified, will use this instead of V0 in ContainerPath
%     {default = []}
%   'DT_outdir': diffusion tensor calculation output directory
%     absolute or relative to ContainerPath
%     {default = 'DTcalc'}
%   'DT_outfix': string attached to diffusion tensor output file names
%     {default = []}
%   'fibers': fiber numbers to generate
%     {default = [101:110,115:123,133:138,141:150]}
%   'subdiv_fibers': vector of fiber subdivision numbers
%     {default = [1014,1024,1231,1232]}
%   'divide_fibers_flag': [0|1] divide fibers into subdivisions
%     {default = 1}
%   'combine_fibers_flag': [0|1] combine fibers (e.g. left and right hemi)
%     {default = 1}
%   'DTIflag': [0|1] require DTI data for atlas tracking
%     otherwise use T1 registration to atlas alone
%     {default = 1}
%   'locflag': [0|1] toggle use of location information alone
%     {default = 0}
%   'resT1flag': [0|1] resample fibers to T1 resolution
%     {default = 0}
%   'xcg_flag': [0|1] generate fiber ROIs excluding CSF and gray matter
%     as defined by FreeSurfer aseg
%     {default = 0}
%   'xcg_suffix': file name string to include for CSF/gray excluded
%     ignored if xcg_flag = 0
%     {default = 'xcg'}
%   'xcg_codes': FreeSurfer aseg codes used to define xcg mask
%     {default = [0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63]}
%   'masksf_flag': [0|1] exclude voxels with multiple fibers
%     Implicitly sets fseg_flag=1
%     {default = 0}
%   'masksf_suffix': file name string to include for multiple fibers excluded
%     ignored if masksf_flag = 0
%     {default = 'masksf'}
%   'fseg_flag': [0|1] generate fiber segmentation volume
%     {default = 1}
%   'save_mgz_flag': [0|1] save fibers in mgz format in addition to sparse
%     {default = 0}
%   'forceflag': overwrite existing output files
%     {default = 0}
%
% Optional Parameters to control selection of data:
%   'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr_resDTI'
%     {default = 'corr_resDTI'}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
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
%   'T1type': which type of T1 series to use as reference
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%   'FSContainerPath': FreeSurfer recon container path
%     used if fnameT1 in regT1 RegInfo is out of date and regT1 was to FreeSurfer
%     and if xcg_flag = 1
%     {default = []}
%
% Optional Parameters to control creation of fiber paths:
%   'create_paths_flag': [0|1] to generate DTI Studio format fiber paths
%     {default = 1}
%   'thresh_FA_flag': [0|1] toggle application of FA threshold of 0.15
%     FA threshold affects probability threshold required to keep
%       average fiber volumes similar to manual fibers that went into atlas
%     Only affects fiber path generation (ignored if create_paths_flag=0)
%     {default = 0}
%   'thresh_FA': FA threshold applied if thresh_FA_flag = 1
%     {default = 0.15}
%   'thresh_prob': fiber probability threshold for fiber path generation
%     if empty, will bet set depending on thresh_FA_flag and other input
%       assuming thresh_FA = 0.15
%       e.g. 0.08 if thresh_FA_flag = 0 and locflag = 0
%            0.10 if thresh_FA_flag = 0 and locflag = 1
%            0.07 if thresh_FA_flag = 1 and locflag = 0
%            0.08 if thresh_FA_flag = 1 and locflag = 1
%     if thresh_FA or thresh_prob are set directly, fiber volumes will vary
%     {default = []}
%   'min_fiberlen': minimum fiber length for fiber path generation
%     (ignored if create_paths_flag=0)
%     {default = 12}
%   'thresh_angle': maximum turning angle for fiber path generation
%     (ignored if create_paths_flag=0)
%     {default = 70}
%
% Optional Parameters to specify atlases:
%   'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%   'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%   'regFA_flag': [0|1] register FA to atlas instead of T1
%     if 1, be sure to specify an FA atlas for atlasname
%     {default = 0}
%   'fiber_atlasdir': full path containing fiber atlas files
%     {default = [getenv('MMPS_DIR') '/atlases/DTI_Atlas/AllSubjects']}
%   'fiber_atlasname': name of fiber atlas appended to mapsdir and pathsdir
%     {default = []}
%
% Optional Parameters to specify output:
%   'outdir': output directory
%     may be full path, otherwise relative to ContainerPath
%     {default = 'AtlasTrack'}
%   'mapsdir': output directory for fiber masks
%     may be full path, otherwise relative to outdir
%     {default = 'fiber_maps'}
%   'pathsdir': output directory fiber paths
%     may be full path, otherwise relative to outdir
%     {default = 'fiber_paths'}
%   'fseg_fstem': output file stem for fiber segmentation volume
%     {default = 'fseg'}
%
% Output:
%   errcode: 0 if success, 1 if error
%
% Created:  04/27/07 by Don Hagler
% Last Mod: 08/03/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters and check for problems
[parms,errcode] = check_input(ContainerPath,varargin);
if errcode, return; end;

% check that DTI data exists
[parms,errcode] = check_data(parms);
if errcode, return; end;

% check that DTI to T1 registration file exists
[parms,errcode] = check_regT1(parms);
if errcode, return; end;

% set thresh_FA and thresh_prob depending on input parameters
[parms,errcode] = set_track_parms(parms);
if errcode, return; end;

% nonlinear registration to atlas and generate fibers from atlas
[parms,errcode] = track_fibers(parms);
if errcode, return; end;

if ~isempty(parms.snums)
  errcode = transform_fibers(parms);
  if errcode, return; end;
end;

if parms.resT1flag
  errcode = transform_fibers(parms,1);
  if errcode, return; end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_input(ContainerPath,args)
  errcode = 0;
  parms_filter = {...
    'fname_FA',[],[],...
    'fname_V0',[],[],...
    'DT_outdir','DTcalc',[],...
    'DT_outfix',[],[],...
    'fibers',[101:110,115:123,133:138,141:150],[],...
    'subdiv_fibers',[1014,1024,1231,1232],[],...
    'divide_fibers_flag',true,[false true],...
    'combine_fibers_flag',true,[false true],...
    'DTIflag',true,[false true],...
    'locflag',false,[false true],...
    'resT1flag',false,[false true],...
    'xcg_flag',false,[false true],...
    'xcg_suffix','xcg',[],...
    'xcg_codes',[0,24,4,5,14,15,43,44,72,75,76,3,8,42,47,31,63],[],...
    'masksf_flag',false,[false true],...
    'masksf_suffix','masksf',[],...
    'fseg_flag',true,[false true],...
    'save_mgz_flag',false,[false true],...
    'forceflag',false,[false true],...
... % control selection of data
    'snums',[],[],...
    'infix','corr_resDTI',[],...
    'revflag',0,[0 1 2],...
    'min_ndirs',6,[],...
    'min_bval',1,[],...
    'max_bval',Inf,[100,Inf],...
    'flex_flag',false,[false true],...
    'min_nb0',1,[],...
    'nob0_flag',false,[false true],...
    'T1type',2,[0:3],...
    'FSContainerPath',[],[],...
... % control creation of fiber paths:
    'create_paths_flag',true,[false true],...
    'thresh_FA_flag',false,[false true],...
    'thresh_FA',0.15,[],...
    'thresh_prob',[],[],...
    'min_fiberlen',12,[1 1000],...
    'thresh_angle',70,[30 90],...
... % specify atlases
    'atlasdir',[],[],...
    'atlasname','T1_Atlas/T1_atlas',[],...    
    'fiber_atlasdir',[],[],...
    'fiber_atlasname',[],[],...
    'regFA_flag',false,[false true],...
... % specify output directories
    'outdir','AtlasTrack',[],...
    'mapsdir','fiber_maps',[],...
    'pathsdir','fiber_paths',[],...
    'fseg_fstem','fseg',[],...
... % hidden parameters
    'ContainerPath',ContainerPath,[],...
    'fname_T1',[],[],...
    'first_only_flag',1,[],... % tensor atlas generated using first eigen vectors only
    'tensor_smooth_sigma',5,[],... % smoothing applied in generation of atlas
    'countflag',true,[false true],... % for atlas_flag > 0
    'count_flag',true,[false true],... % for atlas_flag = 0
    'fnamestem','DTI',[],...
    'orient_ref','LPS',[],...
    'prob_exponent',1,[],...
    'trans_fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'fseg_fibers',[101:110,115:120,123,135:138,143:150],[],... % excluded Fmaj (121), Fmin (122), SLF (133,134), SCS (141,142)
    'roicode_base',10000,[0,Inf],... % value added to fiber codes for fseg
    'fseg_thresh_prob',0.08,[],...
... % hidden parameters for registration to atlas
    'smoothflag',true,[false true],...
    'sampling',[4 4 4],[],...
    'nK',[5 5 5],[],...
    'tstep',0.5,[],...
    'astep',0.25,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
    'ns',64,[],...
    'sf',1,[],...
    'thresh',20,[0,Inf],...
    'stdflag',true,[false true],...
    'maskflag',true,[false true],...
    'stdbgval',75,[],...
... % parameters to pass to various functions
    'track_tags',{'fname_T1','fname_FA','fname_V0','M_T1_to_DTI',...
                  'fibers','subdiv_fibers','divide_fibers_flag',...
                  'combine_fibers_flag',...
                  'create_paths_flag','save_mgz_flag','forceflag',...
                  'thresh_FA','thresh_prob','min_fiberlen','thresh_angle',...
                  'atlasdir','atlasname','fiber_atlasdir','fiber_atlasname',...
                  'outdir','mapsdir','pathsdir','first_only_flag',...
                  'tensor_smooth_sigma','countflag','orient_ref',...
                  'prob_exponent','regFA_flag','locflag',...
                  'DTIflag','ext','smoothflag','sampling','nK',...
                  'tstep','astep','scales','ns','sf','thresh',...
                  'stdflag','maskflag','stdbgval'},[],...
    'fstem_tags',{'snums','infix','revflag','min_bval','max_bval','flex_flag',...
                  'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
    'reg_tags',{'infix','revflag'},[],...
    'T1_tags',{'fname_T1','T1type','FSContainerPath'},[],...
    'transform_tags',{'resT1flag','xcg_flag','masksf_flag','fseg_flag',...
                      'fseg_resT1flag','fseg_xcg_flag','fseg_thresh_prob',...
                      'create_paths_flag','outdir','fiber_outdir',...
                      'paths_outdir','fname_fseg','fibers','fseg_fibers',...
                      'atlas_flag','M_T1_to_DTI','volsz_T1','M_T1',...
                      'fname_aseg','xcg_suffix','masksf_suffix','fname_FA',...
                      'fname_V0','thresh_FA','thresh_prob','min_fiberlen',...
                      'thresh_angle','path_suffix','save_mgz_flag',...
                      'verbose','forceflag',...
                      'xcg_codes','count_flag','prob_exponent','orient_ref',},[],...
  };
  parms = mmil_args2parms(args,parms_filter);
  
  if isempty(parms.fiber_atlasdir)
    parms.fiber_atlasdir = [getenv('MMPS_DIR') '/atlases/DTI_Atlas/AllSubjects'];  
  end;
  if ~exist(parms.fiber_atlasdir,'dir')
    fprintf('%s: ERROR: atlas fiber directory %s not found\n',...
      mfilename,parms.fiber_atlasdir);
    errcode = 1;
    return;
  end;

  if isempty(parms.atlasdir)
    parms.atlasdir = [getenv('MMPS_DIR') '/atlases'];
  end;
  if ~exist(parms.atlasdir,'dir')
    fprintf('%s: ERROR: atlas fiber directory %s not found\n',...
      mfilename,parms.atlasdir);
    errcode = 1;
    return;
  end;
  
  % check that processed data container exists
  if ~exist(parms.ContainerPath,'dir')
    fprintf('%s: ERROR: %s not found\n',mfilename,parms.ContainerPath);
    errcode = 1;
    return;
  end;

  % check that fname_aseg exists if needed for xcg
  if parms.xcg_flag
    if isempty(parms.FSContainerPath)
      fprintf('%s: ERROR: FSContainerPath required if xcg_flag = 1\n',mfilename);
      errcode = 1;
      return;
    end;
    if ~exist(parms.FSContainerPath,'file')
      fprintf('%s: ERROR: %s not found\n',mfilename,parms.FSContainerPath);
      errcode = 1;
      return;
    end;
    parms.fname_aseg = sprintf('%s/mri/aseg.mgz',parms.FSContainerPath);
    if ~exist(parms.fname_aseg,'file')
      fprintf('%s: ERROR: %s not found\n',mfilename,parms.fname_aseg );
      errcode = 1;
      return;
    end;
  end;

  % set atlas_flag for resT1, xcg, fseg, masksf
  if parms.countflag
    if parms.locflag
      parms.atlas_flag = 1;
    else
      parms.atlas_flag = 2;
    end;
  else
    if parms.locflag
      parms.atlas_flag = 3;
    else
      parms.atlas_flag = 4;
    end;
  end;

  % set output directory
  if mmil_isrelative(parms.outdir)
    parms.outdir = [parms.ContainerPath '/' parms.outdir];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_data(parms)
  errcode = 0;
  % check for DTI data
  tparms = parms;
  tparms.outdir = parms.DT_outdir;
  tparms.outfix = parms.DT_outfix;
  args = mmil_parms2args(tparms,parms.fstem_tags);
  [DT_fstem,parms.snums] = DTI_MMIL_Set_DT_fstem(parms.ContainerPath,args{:});
  if parms.DTIflag && isempty(parms.snums)
    fprintf('%s: ERROR: no DTI data found, skipping AtlasTrack\n',...
      mfilename);
    errcode = 1;
    return;
  elseif ~parms.DTIflag
    fprintf('%s: WARNING: running AtlasTrack without DTI data, using T1 only...\n',...
      mfilename);
    parms.locflag = 1;
  end;
  % check that DT calculations exist
  if ~parms.locflag || parms.DTIflag
    if isempty(parms.fname_FA) || isempty(parms.fname_V0)
      fname_in = [DT_fstem '_meas.mat'];
      if ~exist(fname_in,'file')
        fprintf('%s: ERROR: DT calculations output file %s not found\n',...
          mfilename,fname_in);
        errcode = 1;
        return;
      end;
      parms.fname_FA = [DT_fstem '_FA.mgz'];
      [parms.M_DTI,parms.volsz_DTI] = ...
        mmil_load_mgh_info(parms.fname_FA,parms.forceflag);
      parms.fname_V0 = [DT_fstem '_V0.mgz'];
    end;
    if ~exist(parms.fname_FA,'file')
      fprintf('%s: ERROR: DT FA file %s not found\n',...
        mfilename,parms.fname_FA);
      errcode = 1;
      return;
    end;
    if ~exist(parms.fname_V0,'file')
      fprintf('%s: ERROR: DT V0 file %s not found\n',...
        mfilename,parms.fname_V0);
      errcode = 1;
      return;
    end;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_regT1(parms)
  errcode = 0;
  if ~parms.regFA_flag && (~parms.locflag || parms.DTIflag)
    % load registration info
%    fprintf('%s: loading registration info...\n',mfilename);
    args = mmil_parms2args(parms,parms.reg_tags);
    [RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(parms.ContainerPath,args{:});
    if errcode
      fprintf('%s: ERROR: registration file *_regT1.mat file not found\n',...
        mfilename);
      errcode = 1;
      return;
    end;
    parms.M_T1_to_DTI = RegInfo.M_T1_to_T2;

    % check that T1 file is still there
    parms.fname_T1 = RegInfo.fname_T1;
    if ~exist(parms.fname_T1,'file')
      % maybe container was moved
      [tmp_path,tmp_stem,tmp_ext] = fileparts(parms.fname_T1);
      if ~isempty(parms.FSContainerPath)
        parms.fname_T1 = sprintf('%s/mri/%s%s',parms.FSContainerPath,tmp_stem,tmp_ext);
      else
        parms.fname_T1 = sprintf('%s/%s%s',parms.ContainerPath,tmp_stem,tmp_ext);
      end;
    end;
    if ~exist(parms.fname_T1,'file')
      fprintf('%s: ERROR: T1 file %s not found\n',mfilename,parms.fname_T1);
      errcode = 1;
      return;
    end;
  else
    parms.M_T1_to_DTI = eye(4);
    parms.fname_T1 = [];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = set_track_parms(parms)
  errcode = 0;
  if ~parms.thresh_FA_flag
    parms.thresh_FA = 0;
  end;
  if isempty(parms.thresh_prob)
    if parms.countflag
      if parms.locflag
        if parms.thresh_FA_flag
          parms.thresh_prob = 0.08;
        else
          parms.thresh_prob = 0.10;
        end;
      else
        if parms.thresh_FA_flag
          parms.thresh_prob = 0.07;
        else
          parms.thresh_prob = 0.08;
        end;
      end;
    else
      if parms.locflag
        if parms.thresh_FA_flag
          parms.thresh_prob = 0.28;
        else
          parms.thresh_prob = 0.31;
        end;
      else
        if parms.thresh_FA_flag
          parms.thresh_prob = 0.22;
        else
          parms.thresh_prob = 0.24;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = track_fibers(parms)
  errcode = 0;
  % select reference image
  if ~parms.regFA_flag
    args = mmil_parms2args(parms,parms.T1_tags);
    [parms.fname_T1,errcode] = MMIL_Choose_T1(parms.ContainerPath,args{:});
    if errcode, return; end;
  end;
  % register to atlas
  args = mmil_parms2args(parms,parms.track_tags);
  try
    dti_AtlasTrack(args{:});
  catch
    errcode = 1;
    fprintf('\n%s: WARNING: dti_AtlasTrack failed for %s:\n%s\n\n',...
      mfilename,parms.ContainerPath,lasterr);
    return;
  end;
  % set fiber_dir for transform_fibers
  parms.fiber_dir = [parms.outdir '/' parms.mapsdir];
  if ~isempty(parms.fiber_atlasname)
    parms.fiber_dir = [parms.fiber_dir '_' parms.fiber_atlasname];
  end;
  if ~exist(parms.fiber_dir,'dir')
    errcode = 1;
    fprintf('\n%s: WARNING: dti_AtlasTrack output dir %s not found\n',...
      mfilename,parms.fiber_dir);
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = transform_fibers(parms,resT1flag)
  errcode = 0;
  if ~exist('resT1flag','var'), resT1flag = 0; end;
  parms.resT1flag = resT1flag;
  parms.fiber_outdir = parms.mapsdir;
  parms.paths_outdir = parms.pathsdir;
  if parms.resT1flag
    parms.fiber_outdir = [parms.fiber_outdir '_resT1'];
    parms.paths_outdir = [parms.paths_outdir '_resT1'];
    parms.fname_fseg = [parms.fseg_fstem '_resT1'];
    parms.create_paths_flag = 0; % not currently supported
  else
    parms.fname_fseg = [parms.fseg_fstem '_resDTI'];
  end;
  if parms.xcg_flag
    parms.fname_fseg = [parms.fname_fseg '_' parms.xcg_suffix];
  end;
  parms.fname_fseg = [parms.fname_fseg '.mgz'];
  if parms.atlas_flag && ~isempty(parms.fiber_atlasname)
    parms.fiber_outdir = sprintf('%s_%s',...
      parms.fiber_outdir,parms.fiber_atlasname);
    parms.paths_outdir = sprintf('%s_%s',...
      parms.paths_outdir,parms.fiber_atlasname);
  end;
  parms.fibers = parms.trans_fibers;
  args = mmil_parms2args(parms,parms.transform_tags);
  try
    dti_transform_fibers(parms.fiber_dir,args{:})
  catch
    errcode = 1;
    fprintf('\n%s: WARNING: dti_transform_fibers failed for %s:\n%s\n\n',...
      mfilename,parms.ContainerPath,lasterr);
    return;
  end;
return;

