function dti_AtlasTrack(varargin)
%function dti_AtlasTrack([options])
%
% Purpose: warp subject to atlas and generate white matter
%   tract ROIs in subject space
%
% Required Input:
%   At least one of fname_T1 or fname_FA must be specified (see below)
%
% Optional Parameters to specify input data:
%   'fname_T1': full path of T1 volume (T1-weighted structural)
%     if supplied, will use for nonlinear reg to atlas
%     {default = []}
%   'fname_FA': full path of FA volume (fractional anisotropy)
%     if supplied, will use for nonlinear reg to atlas
%     {default = []}
%   'fname_V0': full path of V0 volume (principal Eigenvector)
%     if supplied, will use to refine fiber tract probability maps
%     {default = []}
%   'M_T1_to_DTI': transformation matrix from T1 to DTI (FA)
%     ignored if fname_T1 not supplied
%     if empty, will use identity matrix
%     {default = []}
%
% Optional Parameters to control what is done:
%   'locflag': [0|1] toggle use of location information alone
%     {default = 0}
%   'fibers': vector of fiber numbers to generate
%     {default = [101:110,115:123,133:138,141:150]}
%   'subdiv_fibers': vector of fiber subdivision numbers
%     {default = [1014,1024,1231,1232]}
%   'divide_fibers_flag': [0|1] divide fibers into subdivisions
%     {default = 1}
%   'combine_fibers_flag': [0|1] combine fibers (e.g. left and right hemi)
%     {default = 1}
%   'create_paths_flag': [0|1] generate DTI Studio format fiber paths
%     {default = 1}
%   'save_mgz_flag': [0|1] save fiber in mgz format in addition to sparse
%     {default = 0}
%   'forceflag': overwrite existing output files
%     {default = 0}
%
% Optional Parameters to control creation of fiber paths:
%   'thresh_FA': FA threshold for fiber path streamline generation
%     {default = 0}
%   'thresh_prob': fiber probability threshold for fiber paths
%     {default = 0.08}
%   'min_fiberlen': minimum fiber length for fiber path generation
%     (ignored if create_paths_flag=0)
%     {default = 12}
%   'thresh_angle': maximum turning angle for fiber path generation
%     (ignored if create_paths_flag=0)
%     {default = 70}
%   'path_suffix': string attached to fiber paths
%     if empty, will attach a long string including information about parameters
%       e.g. pthresh0.02_threshFA0.00_minflen10_loc_countatlas_path
%     {default = []}
%
% Optional Parameters to specify atlases:
%   'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%   'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%   'fiber_atlasdir': full path containing fiber atlas files
%     {default = [getenv('MMPS_DIR') '/atlases/DTI_Atlas/AllSubjects']}
%   'fiber_atlasname': name of fiber atlas appended to mapsdir and pathsdir
%     {default = []}
%   'fiber_subdiv_dir': fiber subdivision directory
%     full path or relative to fiber_atlasdir
%     {default = 'fiber_subdivisions'}
%
% Optional Parameters to specify output directories:
%   'outdir': output directory
%     may be full path, otherwise relative to current directory
%     {default = 'AtlasTrack'}
%   'mapsdir': output directory for fiber masks
%     may be full path, otherwise relative to outdir
%     {default = 'fiber_maps'}
%   'pathsdir': output directory fiber paths
%     may be full path, otherwise relative to outdir
%     {default = 'fiber_paths'}
%
% Created:  04/04/12 by Don Hagler
% Last Mod: 01/04/16 by Don Hagler
%

%% todo: options to control fiber combinations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters and check for problems
parms = check_input(varargin);

% set parameters based on input
parms = set_parms(parms);

% create output directories and copy files
parms = init_output(parms);

% nonlinear registration of T1 to atlas
parms = reg2atlas(parms);

% generate fibers from atlas
track_fibers(parms);

% generate fiber subdivisions from atlas (e.g. left and right corpus callosum)
if parms.divide_fibers_flag
  divide_fibers(parms);
end;

% recombine fibers (e.g. left and right hemisphere)
if parms.combine_fibers_flag
  combine_fibers(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(args)
  parms_filter = {...
    'fname_T1',[],[],...
    'fname_FA',[],[],...
    'fname_V0',[],[],...
    'M_T1_to_DTI',eye(4),[],...
...    
    'locflag',false,[false true],...
    'fibers',[101:110,115:123,133:138,141:150],[],...
    'subdiv_fibers',[1014,1024,1231,1232],[],...
    'divide_fibers_flag',true,[false true],...
    'combine_fibers_flag',true,[false true],...
    'create_paths_flag',true,[false true],...
    'save_mgz_flag',false,[false true],...
    'forceflag',false,[false true],...
... % control creation of fiber paths:
    'thresh_FA',0,[],...
    'thresh_prob',0.08,[],...
    'min_fiberlen',12,[1 1000],...
    'thresh_angle',70,[30 90],...
    'path_suffix',[],[],...
... % specify atlases
    'atlasdir',[],[],...
    'atlasname','T1_Atlas/T1_atlas',[],...    
    'fiber_atlasdir',[],[],...
    'fiber_atlasname',[],[],...
    'fiber_subdiv_dir','fiber_subdivisions',[],...
... % specify output directories
    'outdir','AtlasTrack',[],...
    'mapsdir','fiber_maps',[],...
    'pathsdir','fiber_paths',[],...
... % hidden parameters
    'first_only_flag',1,[],... % tensor atlas generated using first eigen vectors only
    'tensor_smooth_sigma',5,[],... % smoothing applied in generation of atlas
    'countflag',true,[false true],...
    'orient_ref','LPS',[],...
    'prob_exponent',1,[],...
    'regFA_flag',false,[false true],...
    'DTIflag',false,[false true],...
    'ext','.mgh',{'.mgh','.mgz'},...
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
    'warp_tags',{'atlasdir','atlasname'...
                 'smoothflag','sampling','nK','tstep','astep',...
                 'scales','ns','sf','thresh','stdflag','maskflag',...
                 'stdbgval'},[],...
    'track_tags',{'prob_exponent','min_fiberlen','thresh_prob',...
                  'thresh_FA','thresh_angle'},[],...
    'infix_tags',{'atlas_flag','thresh_FA','thresh_prob'},[],...
  };
  parms = mmil_args2parms(args,parms_filter);
  
  if isempty(parms.fiber_atlasdir)
    parms.fiber_atlasdir = [getenv('MMPS_DIR') '/atlases/DTI_Atlas/AllSubjects'];  
  end;
  if ~exist(parms.fiber_atlasdir,'dir')
    error('fiber atlas directory %s not found\n',parms.fiber_atlasdir);
  end;
  if mmil_isrelative(parms.fiber_subdiv_dir)
    parms.fiber_subdiv_dir = [parms.fiber_atlasdir '/' parms.fiber_subdiv_dir];
  end;

  if isempty(parms.atlasdir)
    parms.atlasdir = [getenv('MMPS_DIR') '/atlases'];
  end;
  if ~exist(parms.atlasdir,'dir')
    error('dct morph atlas directory %s not found\n',parms.atlasdir);
  end;

  % check data
  if ~isempty(parms.fname_T1)
    if ~exist(parms.fname_T1,'file')
      error('T1 file %s not found',parms.fname_T1);
    end;
  end;
  if ~isempty(parms.fname_FA)
    if ~exist(parms.fname_FA,'file')
      error('FA file %s not found',parms.fname_FA);
    end;
  end;
  if ~isempty(parms.fname_V0)
    if ~exist(parms.fname_V0,'file')
      error('V0 file %s not found',parms.fname_V0);
    end;
  end;

  %% todo: add these to parms to be set by user
  % set fibers to be combined
  fibers_global = [101:110,115:120,133,134,141,142,147:150,123]; % includes CC
  fibersR = [101,103,105,107,109,115,117,119,133,141,147,149]; % no CC
  fibersL = [102,104,106,108,110,116,118,120,134,142,148,150]; % no CC
  fibersR_withCC = [101,103,105,107,109,115,117,119,133,141,147,149,1231];
  fibersL_withCC = [102,104,106,108,110,116,118,120,134,142,148,150,1232];
  parms.fiber_combinations = ...
    {fibers_global, fibersR, fibersL, fibersR_withCC,fibersL_withCC};
  parms.combined_fibers = [2000,2001,2002,2003,2004];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_parms(parms)
  % get info about input data and set other parameters
  if ~isempty(parms.fname_T1)
    parms.regFA_flag = 0;
  else
    parms.regFA_flag = 1;
  end;
  if isempty(parms.fname_FA) && parms.thresh_FA>0
    fprintf('%s: WARNING: no FA threshold will be used because fname_FA not specified\n',...
      mfilename);
    parms.thresh_FA = 0;
  end;
  if ~parms.locflag && isempty(parms.fname_V0)
    fprintf('%s: WARNING: setting locflag=1 because fname_V0 not specified\n',...
      mfilename);
    parms.locflag = 1;
  end;
  if isempty(parms.fname_FA) && isempty(parms.fname_T1)
    error('either fname_T1 or fname_FA must be specified');
  end;
  if isempty(parms.fname_FA) && isempty(parms.fname_V0)
    parms.DTIflag = 0;
    parms.M_T1_to_DTI = eye(4);
  else
    parms.DTIflag = 1;
  end;

  % set parameters for registering to FA atlas
  if parms.regFA_flag
    if parms.thresh>0.5, parms.thresh = 0; end;
    parms.maskflag = false;
    parms.stdflag = false;
  end;

  % set fiber_infix for paths
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
  args = mmil_parms2args(parms,parms.infix_tags);
  parms.fiber_infix = dti_set_fiber_infix(args{:});  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = init_output(parms)
  % set output directories
  if mmil_isrelative(parms.outdir)
    parms.outdir = [pwd '/' parms.outdir];
  end;
  if mmil_isrelative(parms.mapsdir)
    parms.mapsdir = [parms.outdir '/' parms.mapsdir];
  end;
  if mmil_isrelative(parms.pathsdir)
    parms.pathsdir = [parms.outdir '/' parms.pathsdir];
  end;
  if ~isempty(parms.fiber_atlasname)
    parms.mapsdir = [parms.mapsdir '_' parms.fiber_atlasname];
    parms.pathsdir = [parms.pathsdir '_' parms.fiber_atlasname];
  end;
  
  % create output directories
  mmil_mkdir(parms.outdir);
  mmil_mkdir(parms.mapsdir);
  if parms.create_paths_flag
    mmil_mkdir(parms.pathsdir);
  end;

  % make local copies of T1, FA, and V0, get info for each
  if ~isempty(parms.fname_T1)
    parms.fname_T1_orig = parms.fname_T1;
    parms.fname_T1 = sprintf('%s/T1%s',parms.outdir,parms.ext);
    if ~exist(parms.fname_T1,'file')
      fprintf('%s: copying %s to %s...\n',...
        mfilename,parms.fname_T1_orig,parms.fname_T1);
      [vol,parms.M_T1,mrp,parms.volsz_T1] = fs_load_mgh(parms.fname_T1_orig);
      fs_save_mgh(vol,parms.fname_T1,parms.M_T1,mrp);
    else
      [parms.M_T1,parms.volsz_T1] = mmil_load_mgh_info(parms.fname_T1);
    end;
  end;
  if ~isempty(parms.fname_FA)
    parms.fname_FA_orig = parms.fname_FA;
    parms.fname_FA = sprintf('%s/FA%s',parms.outdir,parms.ext);
    if ~exist(parms.fname_FA,'file')
      fprintf('%s: copying %s to %s...\n',...
        mfilename,parms.fname_FA_orig,parms.fname_FA);
      [vol,parms.M_FA,mrp,parms.volsz_FA] = fs_load_mgh(parms.fname_FA_orig);
      fs_save_mgh(vol,parms.fname_FA,parms.M_FA,mrp);
    else
      [parms.M_FA,parms.volsz_FA] = mmil_load_mgh_info(parms.fname_FA);
    end;
  end;
  if ~isempty(parms.fname_V0)
    parms.fname_V0_orig = parms.fname_V0;
    parms.fname_V0 = sprintf('%s/V0%s',parms.outdir,parms.ext);
    if ~exist(parms.fname_V0,'file')
      fprintf('%s: copying %s to %s...\n',...
        mfilename,parms.fname_V0_orig,parms.fname_V0);
      [vol,parms.M_V0,mrp,parms.volsz_V0] = fs_load_mgh(parms.fname_V0_orig);
      fs_save_mgh(vol,parms.fname_V0,parms.M_V0,mrp);
    else
      [parms.M_V0,parms.volsz_V0] = mmil_load_mgh_info(parms.fname_V0);
    end;
  end;

  if ~isempty(parms.fname_FA) && ~isempty(parms.fname_V0)
    if any(parms.volsz_FA(1:3) ~= parms.volsz_V0(1:3))
      error('FA and V0 volume sizes do not match');
    end;
  end;
  if ~isempty(parms.fname_FA)
    parms.M_DTI = parms.M_FA;
    parms.volsz_FA = parms.volsz_FA(1:3);
    parms.volsz_DTI = parms.volsz_FA;
  elseif ~isempty(parms.fname_V0)
    parms.M_DTI = parms.M_V0;
    parms.volsz_DTI = parms.volsz_V0(1:3);
  else
    parms.M_DTI = parms.M_T1;
    parms.volsz_T1 = parms.volsz_T1(1:3);
    parms.volsz_DTI = parms.volsz_T1;
  end;

  if ~isempty(parms.fname_T1)
    if any(parms.M_DTI(:) ~= parms.M_T1(:)) ||...
       any(parms.volsz_DTI(:)~=parms.volsz_T1(:))
      parms.fname_T1_resDTI = sprintf('%s/T1_resDTI%s',parms.outdir,parms.ext);
      if ~exist(parms.fname_T1_resDTI,'file') || parms.forceflag
        fprintf('%s: resampling %s to %s...\n',...
          mfilename,parms.fname_T1,parms.fname_T1_resDTI);
        ctx_vol_T1 = ctx_load_mgh(parms.fname_T1_orig);
        ctx_vol_DTI = ctx_mgh2ctx(zeros(parms.volsz_DTI),parms.M_DTI);
        ctx_vol_T1_resDTI = ...
          vol_resample_pad(ctx_vol_T1,ctx_vol_DTI,inv(parms.M_T1_to_DTI));
        ctx_save_mgh(ctx_vol_T1_resDTI,parms.fname_T1_resDTI);
        clear ctx_vol_T1 ctx_vol_DTI ctx_vol_T1_resDTI;
      end;
    end;
  end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = reg2atlas(parms)
  % select reference image
  if parms.regFA_flag % FA
    parms.fname_ref = parms.fname_FA;
  else % T1-weighted MRI volume file
    parms.fname_ref = parms.fname_T1;
  end;
  % register to atlas
  args = mmil_parms2args(parms,parms.warp_tags);
  [fname_atl,fname_reg,parms.fname_vxl]=...
    mmil_warp_to_atlas(parms.fname_ref,'vxlmap_flag',1,args{:});
  clear M_Atlas_to_Subj;
  [fpath,fstem,fext] = fileparts(fname_reg);
  fname_reg_rigid = sprintf('%s/%s_rigid%s',fpath,fstem,fext);
  if ~exist(fname_reg_rigid,'file') || parms.forceflag
    load(fname_reg);
    save(fname_reg_rigid,'M_Atlas_to_Subj');
  else
    load(fname_reg_rigid)
  end;
  parms.M_Atlas_to_Subj = M_Atlas_to_Subj;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function track_fibers(parms)
  % morph atlas fiber data to subject for each fiber
  vol_V0 = [];
  for f=parms.fibers
    [fname_fiber_atlas,fname_tensor_atlas] = set_atlas_fnames(parms,f);
    if ~exist(fname_fiber_atlas,'file')
      fprintf('%s: WARNING: fiber atlas probability file %s not found\n',...
        mfilename,fname_fiber_atlas);
      continue;
    end;
    if ~exist(fname_tensor_atlas,'file')
      fprintf('%s: WARNING: fiber atlas tensor file %s not found\n',...
        mfilename,fname_tensor_atlas);
      continue;
    end;

    % set output file names
    [fname_fiber,fname_tensor,fname_prob_dir,fname_prob] =...
       set_output_fnames(parms,f);

    % warp fiber probability from atlas to subject
    warp_fiber_from_atlas(parms,fname_fiber_atlas,fname_fiber);

    % warp atlas tensor to subj T1
    if ~parms.locflag || parms.create_paths_flag
      warp_tensor_from_atlas(parms,fname_tensor_atlas,fname_tensor);
    end;

    if ~parms.locflag && ~exist(fname_prob,'file') || parms.forceflag
      % load principal Eigenvector
      if isempty(vol_V0)
        fprintf('%s: loading V0 volume %s...\n',mfilename,parms.fname_V0);
        vol_V0 = fs_load_mgh(parms.fname_V0);
      end
      % calculate probability from direction
      calc_prob_dir(parms,vol_V0,fname_fiber,fname_tensor,...
        fname_prob_dir,fname_prob);
    end

    % generate FA and V0 from tensor
    if parms.create_paths_flag && ...
        (~parms.DTIflag || isempty(parms.fname_V0))
      calc_DTmeas_from_atlas(parms,f,fname_tensor);
    end;

    % generate streamlines within atlas-derived masks
    if parms.create_paths_flag
      if parms.locflag
        create_paths(parms,fname_fiber,f);
      else
        create_paths(parms,fname_prob,f);
      end;
    end;

    % make copies of fibers in mgz format
    if parms.save_mgz_flag
      mmil_sparse2mgh(fname_fiber,...
        regexprep(fname_fiber,'.mat','.mgz'),parms.forceflag);
      if ~parms.locflag
        mmil_sparse2mgh(fname_prob_dir,...
          regexprep(fname_prob_dir,'.mat','.mgz'),parms.forceflag);
        mmil_sparse2mgh(fname_prob,...
          regexprep(fname_prob,'.mat','.mgz'),parms.forceflag);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function warp_fiber_from_atlas(parms,fname_fiber_atlas,fname_fiber,maskflag)
  if ~exist('maskflag','var') || isempty(maskflag), maskflag = 0; end;
  if ~exist(fname_fiber,'file') || parms.forceflag
    [vol_fiber_atlas,M_fiber_atlas]=mmil_load_sparse(fname_fiber_atlas);
    ctx_vol_fiber_atlas = ctx_mgh2ctx(vol_fiber_atlas,M_fiber_atlas);
    if ~exist('vxlmap','var') | isempty(vxlmap), load(parms.fname_vxl); end;
    fprintf('%s: warping atlas fiber %s to subject...\n',...
      mfilename,fname_fiber_atlas);
    % warp from atlas to T1
    ctx_vol_fiber_atlas.maxI = 1;
    ctx_vol_fiber_atlas.minI = 0;
    interpm = 1; padding = 1; bclamp = 1;
    ctx_vol_prob_fiber = getVolFromVxlMap_amd(ctx_vol_fiber_atlas,vxlmap,1,...
      interpm,padding,bclamp);
    % resample from T1 atlas to DTI (or T1 subj)
    ctx_vol_DTI = ctx_mgh2ctx(zeros(parms.volsz_DTI),parms.M_DTI);
    ctx_vol_prob_fiber = vol_resample_pad(ctx_vol_prob_fiber,ctx_vol_DTI,...
      inv(parms.M_T1_to_DTI),interpm,bclamp);
    ctx_vol_prob_fiber.imgs(find(isnan(ctx_vol_prob_fiber.imgs))) = 0;
    [vol_prob_fiber,M] = ctx_ctx2mgh(ctx_vol_prob_fiber);
    if maskflag
       vol_prob_fiber = 1.0 * (vol_prob_fiber > eps); % make it a binary mask
    end;
    % save unwarped, resampled fiber probability distribution
    mmil_save_sparse(vol_prob_fiber,fname_fiber,M);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function warp_tensor_from_atlas(parms,fname_tensor_atlas,fname_tensor)
  if ~exist(fname_tensor,'file') || parms.forceflag
    [vol_tensor_atlas,M_atlas] = mmil_load_sparse(fname_tensor_atlas);
    if ~exist('vxlmap','var') | isempty(vxlmap), load(parms.fname_vxl); end;
    fprintf('%s: warping atlas tensor %s to subject...\n',...
      mfilename,fname_tensor_atlas);
    vol_tensor_subj = [];
    num_tensor_components = size(vol_tensor_atlas,4);
    vol_tensor_subj = zeros([parms.volsz_DTI,num_tensor_components]);
    for t=1:num_tensor_components
      ctx_vol = ctx_mgh2ctx(squeeze(vol_tensor_atlas(:,:,:,t)),M_atlas);
      % warp from atlas to T1
      ctx_vol.maxI = 1;
      ctx_vol.minI = -1;
      interpm = 1; padding = 1; bclamp = 0;
      ctx_vol_subj = getVolFromVxlMap_amd(ctx_vol,vxlmap,1,...
        interpm,padding,bclamp);
      % resample from T1 atlas to DTI (or T1 subj)
      ctx_vol_DTI = ctx_mgh2ctx(zeros(parms.volsz_DTI),parms.M_DTI);
      ctx_vol_subj = vol_resample_pad(ctx_vol_subj,ctx_vol_DTI,...
        inv(parms.M_T1_to_DTI),interpm,bclamp);
      ctx_vol_subj.imgs(find(isnan(ctx_vol_subj.imgs))) = 0;
      vol_tensor_subj(:,:,:,t) = ctx_vol_subj.imgs;
    end;
    % rotate coordinate system for tensors from atlas to subj (M_Atlas_to_DTI)
    M_Atlas_to_DTI = parms.M_T1_to_DTI*parms.M_Atlas_to_Subj;
    vol_tensor_subj = dti_rotate_tensors(vol_tensor_subj,M_Atlas_to_DTI);
    % save unwarped, resampled tensor probability distribution
    mmil_save_sparse(vol_tensor_subj,fname_tensor,parms.M_DTI);
    clear ctx_vol_subj ctx_vol;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_prob_dir(parms,vol_V0,fname_fiber,fname_tensor,...
                       fname_prob_dir,fname_prob)
  clear vol_prob_dir;
  if ~exist(fname_prob_dir,'file') || parms.forceflag
    fprintf('%s: calculating probability from direction to make %s...\n',...
        mfilename,fname_prob_dir);
    % load fiber location probability
    vol_prob_fiber = mmil_load_sparse(fname_fiber);
    % load fiber tensor
    vol_tensor_subj = mmil_load_sparse(fname_tensor);
    % initialize variables
    volsz = parms.volsz_DTI;
    vol_prob_dir = zeros(volsz);
    vol_tensor_std = std(vol_tensor_subj,0,4);
    voxels = find(vol_tensor_std>eps);
    nvox_ROI = length(voxels);
    fprintf('%s: %d selected voxels\n',mfilename,nvox_ROI);
    for m=1:nvox_ROI
      ind=voxels(m);
      [i,j,k] = ind2sub(volsz,ind);
      T = dti_components2tensor(squeeze(vol_tensor_subj(i,j,k,:)));
      v = squeeze(vol_V0(i,j,k,:)); 
      vol_prob_dir(i,j,k) = abs(v'*T*v);
    end;
    mmil_save_sparse(vol_prob_dir,fname_prob_dir,parms.M_DTI);
    clear vol_tensor_std;
  end;
  % calculate combined probability
  if ~exist(fname_prob,'file') || parms.forceflag
    fprintf('%s: calculating combined probability to make %s...\n',...
        mfilename,fname_prob);
    if ~exist('vol_prob_dir','var')
      vol_prob_dir = mmil_load_sparse(fname_prob_dir);
    end;
    if ~exist('vol_prob_fiber','var')
      vol_prob_fiber = mmil_load_sparse(fname_fiber);
    end;
    vol_prob_comb = vol_prob_fiber.*vol_prob_dir; % could weight prob_fiber and prob_dir with exponents
    mmil_save_sparse(vol_prob_comb,fname_prob,parms.M_DTI);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_DTmeas_from_atlas(parms,f,fname_tensor);
  % generate FA and V0 from tensor
  [fname_FA_atlas,fname_V0_atlas] = set_DTmeas_atlas_fnames(parms,f);
  if ~exist(fname_FA_atlas,'file') ||...
     ~exist(fname_V0_atlas,'file') || parms.forceflag
    vol_tensor_subj = mmil_load_sparse(fname_tensor);
    DTmeas = dti_calc_eig(vol_tensor_subj);
    if ~exist(fname_FA_atlas,'file') || parms.forceflag
      mmil_save_sparse(DTmeas.volFA,fname_FA_atlas,parms.M_DTI);
    end;
    if ~exist(fname_V0_atlas,'file') || parms.forceflag
      mmil_save_sparse(squeeze(DTmeas.volV(:,:,:,1,:)),...
        fname_V0_atlas,parms.M_DTI);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function divide_fibers(parms)
  for f=1:length(parms.subdiv_fibers)
    subdiv_fibernum = parms.subdiv_fibers(f);
    fibernum = floor(subdiv_fibernum/10);
    % corpus callosum (123) has extra set of subdivisions 1240 - 1248
    if fibernum == 124, fibernum = 123; end;
    % do not divide fibers if not included in fibers
    if ~ismember(fibernum,parms.fibers), continue; end;    
    % output filenames from track_fibers (used as input for this function)
    [fname_fiber,fname_tensor,fname_prob_dir,fname_prob] = ...
      set_output_fnames(parms,fibernum);
    % output filenames
    [fname_fiber_split,fname_tensor_split,...
     fname_prob_dir_split,fname_prob_split] = ...
      set_output_fnames(parms,subdiv_fibernum);
    if ~parms.locflag
      fname_in = fname_prob;
      fname_out = fname_prob_split;
    else
      fname_in = fname_fiber;
      fname_out = fname_fiber_split;
    end;
    if ~exist(fname_out,'file') || parms.forceflag
      % warp atlas-space subdivision to subject
      fname_fibermask_atl = sprintf('%s/fiber_%04d_mask.mat',...
        parms.fiber_subdiv_dir,subdiv_fibernum);
      fname_fibermask_subj = sprintf('%s/fiber_%04d_mask.mat',...
        parms.mapsdir,subdiv_fibernum);
      warp_fiber_from_atlas(parms,fname_fibermask_atl,fname_fibermask_subj,1);
      % mask fname_prob or fname_fiber using the fiber mask
      vol_in = mmil_load_sparse(fname_in); 
      vol_mask = mmil_load_sparse(fname_fibermask_subj); 
      vol_out = vol_in .* vol_mask;
      clear vol_in vol_mask;
      mmil_save_sparse(vol_out,fname_out,parms.M_DTI);
      clear vol_out;
    end;
    % make copies of fibers in mgz format
    if parms.save_mgz_flag
      mmil_sparse2mgh(fname_out,...
        regexprep(fname_out,'.mat','.mgz'),parms.forceflag);
    end;
    % generate streamlines within atlas-derived masks
    if parms.create_paths_flag
      create_paths(parms,fname_out,subdiv_fibernum);
    end; 
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combine_fibers(parms)
  % combining left and right hemisphere fibers
  for i = 1:length(parms.fiber_combinations)
    fiber_list = parms.fiber_combinations{i};
    nfibers = length(fiber_list);
    % do not combine fibers if they are not all in fibers or subdiv_fibers
    if parms.divide_fibers_flag
      all_fibers = union(parms.fibers,parms.subdiv_fibers);
    else
      all_fibers = parms.fibers;
    end;
    if ~all(ismember(fiber_list,all_fibers)), continue; end;
    fibernum = parms.combined_fibers(i);
    [fname_fiber_comb,fname_tensor_comb,fname_prob_dir_comb,fname_prob_comb] = ...
      set_output_fnames(parms,fibernum);
    if ~parms.locflag
      fname_out = fname_prob_comb;
    else
      fname_out = fname_fiber_comb;
    end;
    if ~exist(fname_out,'file') || parms.forceflag
      vol_comb = zeros([parms.volsz_DTI,nfibers]);
      for f = 1:length(fiber_list)
        fnum = fiber_list(f);
        [fname_fiber,fname_tensor,fname_prob_dir,fname_prob] = ...
            set_output_fnames(parms,fnum);
        if ~parms.locflag
          fname_in = fname_prob;
        else
          fname_in = fname_fiber;
        end;
        vol = mmil_load_sparse(fname_in);
        vol_comb(:,:,:,f) = vol;
        clear vol;
      end;
      vol_comb = max(vol_comb,[],4);
      mmil_save_sparse(vol_comb,fname_out,parms.M_DTI);
      clear vol_comb;
    end;
    % generate streamlines within atlas-derived masks
    if parms.create_paths_flag
      create_paths(parms,fname_out,fibernum);
    end; 
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_paths(parms,fname_prob,f)
  [tmp_path,tmp_fname,tmp_ext] = fileparts(fname_prob);
  local_fname_prob = [tmp_fname,tmp_ext];
  if isempty(parms.path_suffix)
    outstem = sprintf('%s/fiber_%02d_%s_minflen%d',...
      parms.pathsdir,f,parms.fiber_infix,parms.min_fiberlen);
    outstem = [outstem '_path'];
  else
    outstem = sprintf('%s/fiber_%02d_%s',parms.pathsdir,f,parms.path_suffix);
  end;    
  fname_out = sprintf('%s.grp',outstem);
  if ~exist(fname_out,'file') || parms.forceflag
    args = mmil_parms2args(parms,parms.track_tags);
    orient = fs_read_orient([],parms.M_DTI);
    [permvec,flipvec] = fs_compare_orient(orient,parms.orient_ref);
    flipflags = flipvec<0;
    if ~parms.DTIflag
      [fname_FA,fname_V0] = set_DTmeas_atlas_fnames(parms,f);
    elseif isempty(parms.fname_V0)
      [fname_FA,fname_V0] = set_DTmeas_atlas_fnames(parms,f);
      fname_FA = parms.fname_FA;
    else
      fname_FA = parms.fname_FA;
      fname_V0 = parms.fname_V0;
    end;
    dti_track_fibers(fname_FA,fname_V0,...
      'fname_prob',fname_prob,...
      'outstem',outstem,...
      'flipx_flag',flipflags(1),...
      'flipy_flag',flipflags(2),...
      'flipz_flag',flipflags(3),...
      args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_fiber_atlas,fname_tensor_atlas] = set_atlas_fnames(parms,f)
  if parms.countflag
    fname_fiber_atlas = sprintf('%s/fiber_%03d_mean_countatlas.mat',...
      parms.fiber_atlasdir,f);
  else
    fname_fiber_atlas = sprintf('%s/fiber_%03d_mean_maskatlas.mat',...
      parms.fiber_atlasdir,f);
  end;
  if parms.first_only_flag
    fname_tensor_atlas = sprintf('%s/fiber_%03d_mean_tensorV0_atlas_sm%0.2f_norm.mat',...
      parms.fiber_atlasdir,f,parms.tensor_smooth_sigma);
  else
    fname_tensor_atlas = sprintf('%s/fiber_%03d_mean_tensor_atlas_sm%0.2f_norm.mat',...
      parms.fiber_atlasdir,f,parms.tensor_smooth_sigma);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_fiber,fname_tensor,fname_prob_dir,fname_prob] = set_output_fnames(parms,f)
  if parms.countflag
    suffix = 'countatlas';
  else
    suffix = 'maskatlas';
  end;
  fname_fiber = sprintf('%s/fiber_%03d_%s.mat',parms.mapsdir,f,suffix);
  fname_prob = sprintf('%s/fiber_%03d_prob_%s.mat',parms.mapsdir,f,suffix);
  fname_tensor = sprintf('%s/fiber_%03d_tensor.mat',parms.mapsdir,f);
  fname_prob_dir = sprintf('%s/fiber_%03d_prob_dir.mat',parms.mapsdir,f);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_FA_atlas,fname_V0_atlas] = set_DTmeas_atlas_fnames(parms,f)
  fname_FA_atlas = sprintf('%s/fiber_%03d_FA_atlas.mat',parms.mapsdir,f);
  fname_V0_atlas = sprintf('%s/fiber_%03d_V0_atlas.mat',parms.mapsdir,f);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

