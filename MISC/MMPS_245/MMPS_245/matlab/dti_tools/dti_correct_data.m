function [fname_out,fname_B0dx,fname_qmat] = dti_correct_data(fname,varargin)
%function [fname_out,fname_B0dx,fname_qmat] = dti_correct_data(fname,[options])
%
% Purpose: correct diffusion images for distortions from eddy currents,
%   head motion, B0 inhomogeneity, and grad warp
%
% Required Parameters:
%   fname: full path name of mgh/mgz file containing 4D DTI volume
%
% Optional Control Parameters:
%   'B0unwarp_flag': [0|1|2] whether to correct for B0 distortions
%     if 2, processing aborted if B0 unwarping fails (e.g. no rev scans)
%     {default = 1}
%   'ecc_flag': [0|1] whether to perform eddy current correction
%     {default = 1}
%   'censor_flag': [0|1] automatically reject bad slices based on tensor fit
%      values for bad slices will be replaced with values
%      synthesized from tensor fit
%     {default = 1}
%   'motion_B0uw_flag': [0|1] whether to estimate head motion
%     and apply to B0dx field
%     {default = 1}
%   'mc_flag': [0|1] whether to correct for head motion
%     {default = 1}
%   'gradunwarp_flag': [0|1|2] whether to correct for grad warp
%     {default = 1}
%
% Optional Parameters:
%   'fname_out': output file fname for corrected data
%     If empty, will append '_corr.mgz' to file stem of fname
%     {default = []}
%   'fname_qmat': output file name for adjusted qmat
%     also contains motion estimates
%     If empty, will append '_corr_qmat.mat' to file stem of fname
%     {default = []}
%   'qmat': matrix of diffusion direction vectors (ndirs x 3)
%     required for tensor-based eddy current and motion correction
%     {default = []}
%   'bvals': vector of b values (one for each diffusion direction)
%     If single value supplied, will use same for all
%     {default = [1000]}
%   'fname_censor': text file containing frame numbers to exclude (1-based)
%     {default = []}
%   'censor_thresh': error threshold for censoring bad frames
%     normalized to median error for each slice
%     higher values mean less censoring
%     {default = 3.2}
%   'censor_min_ndirs': minimum number of diffusion directions (not including
%     b=0 images) required for tensor fit after censoring
%     will not do censoring if remaining ndirs < min_ndirs
%     {default = 6}
%   'nonlin_flag': [0|1] use nonlinear optimization for tensor fits
%     {default = 0}
%   'b0_thresh': threshold used for considering a b-value to be 0
%     {default = 10}
%   'interpm': interpolation method for motion correction and grad unwarp
%      0 = nearest neighbor, 1 = linear, 2 = cubic
%      3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%     { default = 2 }
%   'maskoutput': apply dilated brain mask to set background to zero
%     {default = 0}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': overwrite existing output files
%     {default = 0}
%
% Optional Parameters for B0 unwarping:
%   'optimize_B0uw_flag': [0|1] search for optimal B0 unwarp parameters
%     kernelWidthMax and lambda2
%     {default = 0}
%   'inorm_B0uw_flag': [0|1] intensity normalize 'forward' and 'reverse'
%     polarity images {default = 1}
%   'fname_B0dx'; name of input/output B0dx file with estimated displacements
%     if not found, will create using fname_for and fname_rev
%     if empty, will not do B0 unwarping
%     {default = []}
%   'fname_for'; name of input b=0 image with "forward" phase-encode polarity
%     required for B0 unwarping if fname_B0dx not supplied or does not exist
%     {default = []}
%   'fname_rev'; name of input b=0 image with "reverse" phase-encode polarity
%     required for B0 unwarping if fname_B0dx not supplied or does not exist
%     {default = []}
%   'fname_B0uw_ref': full path file name of reference volume (mgh/mgz format)
%     for registering to fname_B0dx
%     {default = []}
%   'fname_B0uw_mask': name of mask volume for registration to fname_B0uw_ref
%     If empty, will generate one from fname_ref using mmil_quick_brainmask
%     {default = []}
%   'fname_B0uw_reg': output file name for registration to fname_B0uw_ref
%     {default = []}
%   'PhaseDir': phase encode direction (obtained from dicom header)
%     'COL' or 'ROW'
%     {default = 'COL'}
%   'revflag': [0|1] whether to treat input data as "reverse" 
%     phase-encode polarity for applying B0 unwarping
%     {default = 0}
%
% Optional Parameters for eddy current correction
%   'driftcorr': [0|1] estimate drift correction with eddy current correction
%     {default = 0}
%
% Optional Parameters for estimating motion as part of B0uw:
%   'motion_B0uw_iters': number of iterations to estimate motion
%     and B0 displacement
%     {default = 2}
%   'min_trans': minimum translation (mm) for estimating motion for B0uw
%     {default = 0.05}
%   'min_rot': minimum rotation (degrees) for estimating motion for B0uw
%     {default = 0.05}
%   'fname_motion_B0uw': output file name of mat file containing
%     estimated motion; if empty, will be temporary
%     {default = []}
%
% Optional Parameters for grad unwarping:
%   'gruw_type': gradient type number
%     0:  Siemens Sonata
%     1:  Siemens Allegra
%     2:  GE BRM
%     3:  GE CRM
%     4:  Siemens Avanto
%     5:  Siemens AXXess/Espree
%     6:  Siemens Quantum/Symphony
%     7:  GE Twin Speed Whole Body
%     8:  GE Twin Speed Zoom
%     9:  GE mr450 or mr750
%     10: GE MR750W
%     11: Siemens Skyra
%     12: Siemens Connectome Skyra
%     13: Siemens Prisma
%     {default = 8}
%   'gruw_unwarpflag' [0|1|2]
%     0: unwarp 3D
%     1: unwarp through plan only
%     2: unwarp inplane only
%     {default = 0}
%   'gruw_isoctrflag': [0|1] whether to adjust for isocenter coordinates
%     {default = 1}
%
% Optional Parameters for motion correction:
%   'fname_ref': full path file name of reference volume (mgh/mgz format)
%     for between-scan motion correction
%     {default = []}
%   'fname_mask': name of mask volume
%     If empty, will generate one from fname_ref using mmil_quick_brainmask
%     {default = []}
%   'fname_reg': output file name for between-scan registration matrix
%     {default = []}
%
% Output:
%   fname_out: mgz file containing 4D diffusion volume corrected for
%     eddy current distortion, head motion, B0 distortion, and gradient warping
%   fname_qmat: mat file containing:
%    'qmat': diffusion direction matrix (ndirs x 3) adjusted for rotation
%      from motion correction
%    'M_motion': cell array containing registration matrices M_ref_to_orig
%      for each frame
%   fname_B0dx: mgz file containing volume of estimated displacement
%     along phase encode direction
%
% Created:  02/23/10 by Don Hagler
% Prev Mod: 03/23/16 by Don Hagler
% Last Mod: 01/08/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_B0dx = []; fname_qmat = [];

if ~mmil_check_nargs(nargin, 1), return; end;

parms = check_input(fname,varargin);

fname_out = parms.fname_out;
fname_B0dx = parms.fname_B0dx;
fname_qmat = parms.fname_qmat;

if output_exists(parms), return; end;
mmil_mkdir(parms.tmpdir);

% estimate distortions due to magnetic field susceptibility inhomogeneities
if parms.B0unwarp_flag
  estimate_B0_distortion(parms);
end;

if ~isempty(parms.fname_censor)
  parms = censor_frames(parms);
end;

if parms.DT_flag
  if parms.ecc_flag
    % correct for distortions due to eddy currents
    parms = correct_eddycurr(parms);
  elseif parms.censor_flag
    % censor slices using tensor fit
    parms = censor_slices(parms);
  end;
  % load DTfit and set updated parameters
  load(parms.fname_DT);
  parms.DTfit = DTfit;
  parms.qmat = DTfit.qmat;
  parms.bvals = DTfit.bvals;
  parms.censor_mat = DTfit.censor_mat;
  if parms.censor_flag
    save_censor_result(parms);
  end;
end;

% correct for head motion
if parms.mc_flag && ~parms.motion_B0uw_flag
  parms = correct_motion(parms,1);
end;

% correct for distortion due to magnetic field susceptibility inhomogeneities
if parms.B0unwarp_flag
  parms = correct_B0_distortion(parms);
end;

% replace censored slices
if parms.censor_flag && (parms.motion_B0uw_flag || ~parms.mc_flag)
  parms = replace_censored_slices(parms);
end;

% correct for distortion due to gradient coil inhomogeneities
if parms.gradunwarp_flag
  parms = correct_gradwarp(parms);
end;

% correct for head motion
if parms.mc_flag || ~isempty(parms.fname_ref)
  parms = correct_motion(parms);
end;

% copy last file to fname_out
fs_copy_mgh(parms.fname_in,parms.fname_out);

% remove temporary files
cleanup_tmpdir(parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
    'fname_orig',fname,[],...
  ...
    'B0unwarp_flag',1,[0:2],...
    'ecc_flag',true,[false true],...
    'censor_flag',true,[false,true],...
    'motion_B0uw_flag',true,[false true],...
    'mc_flag',true,[false true],...
    'gradunwarp_flag',1,[0:2],...
  ...
    'fname_out',[],[],...
    'fname_qmat',[],[],...
    'qmat',[],[],...
    'bvals',1000,[],...
    'fname_censor',[],[],...
    'censor_min_ndirs',6,[],...
    'censor_thresh',3.2,[],...
    'nonlin_flag',false,[false true],...
    'b0_thresh',10,[0,500],...
    'interpm',2,[0:5],...
    'maskoutput',false,[false true],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'optimize_B0uw_flag',false,[false true],...
    'inorm_B0uw_flag',true,[false true],...
    'fname_B0dx',[],[],...
    'fname_for',[],[],...
    'fname_rev',[],[],...
    'fname_B0uw_ref',[],[],...
    'fname_B0uw_mask',[],[],...
    'fname_B0uw_reg',[],[],...
    'PhaseDir','COL',{'ROW','COL'},...
    'revflag',false,[false,true],...
  ...
    'driftcorr',false,[false true],...
  ...
    'motion_B0uw_iters',2,[1:10],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'fname_motion_B0uw',[],[],...
  ...
    'gruw_type',8,[0:13],...
    'gruw_unwarpflag',0,[],...
    'gruw_isoctrflag',true,[false true],...
  ...
    'fname_ref',[],[],...
    'fname_mask',[],[],...
    'fname_reg',[],[],...
  ... % more parameters for rbreg_vol2vol_icorr_mask
    'mstep',4,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
  ... % censoring
    'ecc_censor_niter',3,[0,100],...
    'mc_censor_niter',1,[0,100],...
  ... % grad unwarp
    'gruw_jacobian_flag',false,[false true],...
  ... % B0uw optimization
    'kernelWidthMax',25,[1:100],...
    'lambda2',1100,[1:10000],...
    'kernelWidthMax_vec',[25,31,35],[1:100],...
    'lambda2_vec',[1100,1500,1900],[1:10000],...
    'multi_opt_flag',false,[false true],...
  ... % other
    'outfix','corr',[],...
    'outext','.mgz',{'.mgh','.mgz'},...
    'tmpext','.mgh',{'.mgh','.mgz'},...
    'cleanupflag',true,[false true],...
  ...
    'unwarpB0_tags',{'fname_out','verbose','forceflag',...
                     'optimize_B0uw_flag','inorm_B0uw_flag',...
                     'fname_B0dx','fname_for','fname_rev',...
                     'revflag','PhaseDir',...
                     'motion_B0uw_flag','motion_B0uw_iters','qmat',...
                     'fname_ref','fname_mask','fname_reg','censor_niter',...
                     'censor_min_ndirs','censor_thresh','censor_mat',...
                     'bvals','min_trans','min_rot','fname_motion_B0uw',...
                     'kernelWidthMax','lambda2','kernelWidthMax_vec',...
                     'lambda2_vec','multi_opt_flag','outfix',...
                     'outext','tmpext','cleanupflag'},[],...
    'tensor_tags',{'bvals','M','censor_niter','censor_min_ndirs',...
                   'censor_thresh','censor_mat','nonlin_flag','b0_thresh'},[],...
    'eddy_tags',{'fname_in','fname_out','fname_DT',...
                 'niter','censor_niter','censor_min_ndirs','censor_thresh',...
                 'nonlin_flag','b0_thresh',...
                 'bvals','nparams','phasedir','driftcorr','smf'},[],...
    'motion_tags',{'fname_out','fname_qmat','interpm','verbose',...
                   'forceflag','intra_flag','qmat',...
                   'fname_ref','fname_mask','fname_reg',...
                   'censor_niter','censor_min_ndirs','censor_thresh',...
                   'nonlin_flag','b0_thresh',...
                   'censor_mat','bvals','min_trans','min_rot',...
                   'mstep','scales','smf','motion_radius'},[],...
    'censor_tags',{'volmask','censor_min_ndirs','nb0'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);

  % check input file
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;

  % create output file names
  [tpath,parms.fstem,text] = fileparts(parms.fname_in);
  if isempty(parms.fname_out)
    parms.fname_out = [tpath '/' parms.fstem '_' parms.outfix parms.outext];
  end;
  parms.outdir = fileparts(parms.fname_out);
  parms.tmpdir = [parms.outdir '/tmp_DTI_corr'];
  if isempty(parms.fname_qmat)
    parms.fname_qmat = [tpath '/' parms.fstem '_' parms.outfix '_qmat.mat'];
  end;

  % check input files for B0 distortion correction
  if parms.B0unwarp_flag
    if isempty(parms.fname_B0dx) && ...
       ~isempty(parms.fname_for) && ~isempty(parms.fname_rev)
      [tpath,fstem_rev,text] = fileparts(parms.fname_rev);
      [tpath,fstem_for,text] = fileparts(parms.fname_for);
      parms.fname_B0dx = ...
        [parms.outdir '/' fstem_for 'VS' fstem_rev '_B0dx' parms.outext];
    end;
    if isempty(parms.fname_B0dx)
      if parms.B0unwarp_flag==2
        error('files required for B0 distortion correction not supplied');
      else
        fprintf('%s: WARNING: skipping B0 distortion correction\n',mfilename);
        parms.B0unwarp_flag = 0;
      end;
    else
      % check that fname_for and fname_rev exist
      if ~exist(parms.fname_B0dx,'file')
        if ~exist(parms.fname_for,'file')
          error('file %s not found',parms.fname_for);
        end;
        if ~exist(parms.fname_rev,'file')
          error('file %s not found',parms.fname_rev);
        end;
      end;
    end;
  end;

  % whether phase-encode direction is along row or column (for eddycurr)
  if strcmp(parms.PhaseDir,'ROW')
    parms.phasedir = 1;
  elseif strcmp(parms.PhaseDir,'COL')
    parms.phasedir = 2;
  end;

  % whether diffusion tensor fitting is possible/required
  if isempty(parms.qmat)
    parms.DT_flag = 0;
    parms.censor_flag = 0;
  elseif ~parms.ecc_flag && ~parms.censor_flag
    parms.DT_flag = 0;
  else
    parms.DT_flag = 1;
  end;

  % get information about input file
  if parms.DT_flag || parms.ecc_flag
    [M,volsz] = ...
      mmil_load_mgh_info(parms.fname_in,parms.forceflag);
    nframes = volsz(4);
  end;

  % check for agreement between nframes and qmat ndirs
  if parms.DT_flag
    if nframes~=size(parms.qmat,1)
      fprintf('%s: WARNING: number of frames (%d) in %s does not match length of qmat (%d)\n',...
        mfilename,nframes,parms.fname_in,size(parms.qmat,1));
      parms.DT_flag = 0;
    else
      parms.DT_flag = 1;
    end;
  end;

  % warn if unable to perform eddy current correction
  if ~parms.DT_flag && parms.ecc_flag
    if nframes>1
      fprintf('%s: WARNING: unable to perform eddy current correction for %s\n',...
        mfilename,parms.fname_in);
    end;
    parms.ecc_flag = 0;
  end;

  if parms.DT_flag
    parms.fname_DT = [parms.tmpdir '/' parms.fstem '_DT.mat'];
  else
    parms.fname_DT = [];
  end;

  if parms.verbose
    fprintf('%s: input file %s\n',mfilename,parms.fname_in);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exist_flag = output_exists(parms)
  if parms.forceflag
    exist_flag = 0;
    return;
  else
    exist_flag = 1;
  end;
  if ~isempty(parms.fname_out) && ~exist(parms.fname_out,'file')
    exist_flag = 0;
  end;
  if parms.B0unwarp_flag && ~exist(parms.fname_B0dx,'file')
    exist_flag = 0;
  end;
  if ~isempty(parms.qmat) && ~exist(parms.fname_qmat,'file')
    exist_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estimate_B0_distortion(parms)
  args = mmil_parms2args(parms,parms.unwarpB0_tags);
  dti_unwarpB0([],args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = censor_frames(parms)
  [tpath,fstem,text] = fileparts(parms.fname_in);
  fname_out = [parms.tmpdir '/' fstem '_cf' parms.tmpext];
  fname_mat = [parms.outdir '/' fstem '_cf.mat'];
  if ~exist(fname_out,'file') ||...
     ~exist(parms.fname_mat,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: censoring frames...\n',mfilename);
      tic
    end;
    % load data
    [vol,M] = fs_load_mgh(parms.fname_in);
    % censor bad frames from vol, qmat, and bvals
    [vol,qmat,bvals,ind_censored_frames] = ...
      dti_censor_frames(parms.fname_censor,vol,parms.qmat,parms.bvals);
    % save censored volume
    fs_save_mgh(vol,fname_out,M);
    % save qmat, bvals, ind_censored_frames
    save(fname_mat,'qmat','bvals','ind_censored_frames');
    if parms.verbose, toc; end;
  end;
  parms.fname_in = fname_out;
  load(fname_mat);
  parms.qmat = qmat;
  parms.bvals = bvals;
  parms.ind_censored_frames = ind_censored_frames;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = correct_eddycurr(parms)
  [tpath,fstem,text] = fileparts(parms.fname_in);
  tparms = parms;
  tparms.fname_out = [parms.tmpdir '/' fstem '_ecc' parms.tmpext];
  if ~exist(tparms.fname_out,'file') ||...
     ~exist(tparms.fname_DT,'file') || tparms.forceflag
    if parms.verbose
      fprintf('%s: correcting for eddy current distortion...\n',mfilename);
      tic
    end;
    if parms.censor_flag
      tparms.censor_niter = parms.ecc_censor_niter;
    else
      tparms.censor_niter = 0;
    end;
    args = mmil_parms2args(tparms,parms.eddy_tags);
    dti_eddycurr_corr([],parms.qmat,args{:});
    if parms.verbose, toc; end;
  end;
  parms.fname_in = tparms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = censor_slices(parms)
  if ~exist(parms.fname_DT,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: censoring slices...\n',mfilename);
      tic
    end;
    if parms.censor_flag
      tparms.censor_niter = parms.ecc_censor_niter;
    else
      tparms.censor_niter = 0;
    end;
    % load data
    [vol,M] = fs_load_mgh(parms.fname_in);
    % fit tensor, censoring bad frames for each slice
    args = mmil_parms2args(tparms,parms.tensor_tags);
    DTfit = dti_fit_tensor(vol,parms.qmat,args{:});
    % save tensor fit
    save(parms.fname_DT,'DTfit');
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_censor_result(parms)
  fname_mat = [parms.outdir '/' parms.fstem '_censor.mat'];
  if ~exist(fname_mat,'file') || parms.forceflag
    censor_mat = parms.DTfit.censor_mat;
    censor_err = parms.DTfit.censor_err;
    save(fname_mat,'censor_mat','censor_err');
  end;
  fname_plot_mat = [parms.outdir '/' parms.fstem '_censor_mat.tif'];
  fname_plot_err = [parms.outdir '/' parms.fstem '_censor_err.tif'];
  if ~exist(fname_plot_mat,'file') ||...
     ~exist(fname_plot_err,'file') ||...
     parms.forceflag
    try
      % save censor_mat and censor_err images
      figure; set(gcf,'Visible','off');
      imagesc(parms.DTfit.censor_mat); % green if nothing, blue if anything, red if bad
      print(gcf,'-dtiff',fname_plot_mat);
      close(gcf);
      figure; set(gcf,'Visible','off');
      imagesc(parms.DTfit.censor_err,[0,10]);
      print(gcf,'-dtiff',fname_plot_err);
      close(gcf);
    catch
      fprintf('%: WARNING: failed to save censor_mat and censor_err images:\n%s\n',...
        mfilename,lasterr);
    end;
  end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = correct_B0_distortion(parms)
  [tpath,fstem,text] = fileparts(parms.fname_in);
  tparms = parms;
  tparms.fname_out = [parms.tmpdir '/' fstem '_B0uw' parms.tmpext];
  if ~exist(tparms.fname_out,'file')
    if parms.verbose
      fprintf('%s: correcting for B0 distortions...\n',mfilename);
    end;
    if parms.motion_B0uw_flag
      if isempty(parms.fname_motion_B0uw)
        tparms.fname_motion_B0uw = [parms.tmpdir '/' fstem '_motion_B0uw.mat'];
      end;
      % register to scans used to estimate B0 distortion
      [tparms.fname_ref,tparms.fname_mask,tparms.fname_reg] = ...
        set_fname_ref(parms);
    end;
    % do not do between-scan registration for same scan
    if strcmp(tparms.fname_ref,parms.fname_orig)
      tparms.fname_ref = tparms.fname_in;
    end;
    if parms.censor_flag
      tparms.censor_niter = parms.mc_censor_niter;
    else
      tparms.censor_niter = 0;
    end;
    args = mmil_parms2args(tparms,parms.unwarpB0_tags);
    dti_unwarpB0(parms.fname_in,args{:});
  end;
  parms.fname_in = tparms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_ref,fname_mask,fname_reg] = set_fname_ref(parms)
  if ~isempty(parms.fname_B0uw_ref)
    fname_ref = parms.fname_B0uw_ref;
    fname_mask = parms.fname_B0uw_mask;
    fname_reg = parms.fname_B0uw_reg;
  else
    if ~parms.revflag
      fname_ref = parms.fname_for;
    else
      fname_ref = parms.fname_rev;
    end;
    fname_mask = [];
    fname_reg = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = replace_censored_slices(parms)
  [tpath,fstem,text] = fileparts(parms.fname_in);
  fname_out = [parms.tmpdir '/' fstem '_cens' parms.tmpext];
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: replacing censored slices...\n',mfilename);
      tic;
    end;
    % load fname_in
    [vol,M] = fs_load_mgh(parms.fname_in);
    % fit tensor, censoring bad slices
    parms.censor_niter = parms.mc_censor_niter;
    args = mmil_parms2args(parms,parms.tensor_tags);
    DTfit = dti_fit_tensor(vol,parms.qmat,args{:});
    % synthesize volume using tensor fit
    vol_synth = dti_synth_vol(DTfit);
    % replace censored slices in vol
    tparms = parms;
    tparms.nb0 = length(DTfit.i_b0);
    tparms.volmask = DTfit.volmask;
    args = mmil_parms2args(tparms,parms.censor_tags);
    vol = dti_censor_vol(vol,vol_synth,DTfit.censor_mat,args{:});
    % save censored volume
    fs_save_mgh(vol,fname_out,M);
    if parms.verbose, toc; end;
  end;
  parms.mc_censor_niter = 0; % no need to censor again later
  parms.fname_in = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = correct_motion(parms,preB0uw_flag)
  if ~exist('preB0uw_flag','var') || isempty(preB0uw_flag)
    preB0uw_flag = 0;
  end;
  [tpath,fstem,text] = fileparts(parms.fname_in);
  tparms = parms;
  tparms.fname_out = [parms.tmpdir '/' fstem '_mc' parms.tmpext];
  if ~exist(tparms.fname_out,'file') ||...
     ~exist(parms.fname_qmat,'file') || parms.forceflag
    if parms.censor_flag
      tparms.censor_niter = parms.mc_censor_niter;
    else
      tparms.censor_niter = 0;
    end;
    if preB0uw_flag
      % register to scans used for estimating B0 distortion
      [tparms.fname_ref,tparms.fname_mask,tparms.fname_reg] = ...
        set_fname_ref(parms);
      tparms.fname_qmat = [parms.outdir '/' fstem '_mc_qmat.mat'];
    end;
    % do not do between-scan registration for same scan
    if strcmp(tparms.fname_ref,parms.fname_orig)
      tparms.fname_ref = tparms.fname_in;
    end;
    if ~parms.mc_flag || (~parms.motion_B0uw_flag && ~preB0uw_flag)
      % between-scan registration and resampling only
      tparms.intra_flag = 0;
    else
      tparms.intra_flag = 1;
    end;
    if parms.verbose
      if parms.DT_flag && tparms.intra_flag
        fprintf('%s: correcting for head motion...\n',mfilename);
      else
        fprintf('%s: correcting for head motion (between scan only)...\n',...
          mfilename);
      end;
      tic;
    end;
    args = mmil_parms2args(tparms,parms.motion_tags);
    [fname_out,fname_qmat] = ...
      dti_motion_corr(parms.fname_in,args{:});
    % load qmat adjusted for within-scan motion
    if ~isempty(parms.qmat) && preB0uw_flag
      load(fname_qmat)
      parms.qmat = qmat;
    end;
    if parms.verbose, toc; end;
  end;
  parms.fname_in = tparms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = correct_gradwarp(parms)
  [tpath,fstem,text] = fileparts(parms.fname_in);
  fname_out = [parms.tmpdir '/' fstem '_gruw' parms.tmpext];
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: correcting for gradient warping...\n',mfilename);
      tic
    end;
    epi_gradunwarp(parms.fname_in,...
      'fname_out',fname_out,...
      'gwtype',parms.gruw_type,...
      'unwarpflag',parms.gruw_unwarpflag,...
      'isoctrflag',parms.gruw_isoctrflag,...
      'jacobian_flag',parms.gruw_jacobian_flag,...
      'interpm',parms.interpm,...
      'forceflag',parms.forceflag);
    if parms.verbose, toc; end;
  end;
  parms.fname_in = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if parms.cleanupflag && exist(parms.tmpdir,'dir')
    % delete tmp dir
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,result] = unix(cmd);
    if status
      error('deletion of temporary directory %s failed:\n%s\n',...
        parms.tmpdir,result);
    end;
  end;
return;

