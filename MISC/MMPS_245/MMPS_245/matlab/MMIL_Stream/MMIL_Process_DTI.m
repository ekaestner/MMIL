function [infix,errcode] = MMIL_Process_DTI(ContainerPath,varargin)
%function [infix,errcode] = MMIL_Process_DTI(ContainerPath,[options])
%
% Required Input:
%   ContainerPath: full path of DTIPROC Container
%
% Optional Control Parameters:
%   'B0unwarp_flag': [0|1|2] whether to correct for B0 distortions
%     if 2, processing aborted if B0 unwarping fails (e.g. no rev scans)
%     {default = 1}
%   'ecc_flag': [0|1] whether to perform eddy current correction
%     {default=1}
%   'censor_flag': [0|1] reject bad slices based on tensor fit
%     {default = 1}
%   'motion_B0uw_flag': [0|1] whether to estimate head motion
%     and apply to B0dx field
%     {default = 1}
%   'mc_flag': [0|1] whether to perform motion correction
%     {default=1}
%   'gradunwarp_flag': [0|1|2] whether to correct for grad warp
%     {default = 1}
%   'resample_flag': [0|1] whether to resample diffusion data
%     e.g. to make isotropic, to correct inter-scan motion, to register to T1
%     {default = 1}
%
% Optional Parameters for Resampling:
%   'native_flag': resample but keep original nvoxels and resolution
%     {default = 0}
%   'deoblique_flag': [0|1] whether to resample oblique slices to on-axis
%     ignored if native_flag = 1
%     {default = 0}
%   'nvoxels': vector of number of resampled voxels [nx,ny,nz]
%     {default = [120 120 70]}
%   'resolution': desired voxel sizes of resampled data (x y z)
%     {default = [2 2 2]}
%   'std_orient' : specify the resampled output orientation 
%     {default = []}  
%   'rot': rotation applied to processed diffusion data (x,y,z deg)
%     order of rotations: 'x','y','z'
%     if regT1flag=2, suggested values are [-15 0 0]
%     {default: [0,0,0]}
%   'trans': translation applied to processed diffusion data (x,y,z mm)
%     {default = [0,0,0]}
%     if regT1flag=2, suggested values are [0 -15 -25]
%   'smooth': isotropic smoothing kernel sigma (in voxels)
%     {default = 0}
%
% Other Optional Parameters:
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     0: use only forward phase-encode polarity data
%     1: use only reverse phase-encode polarity data
%     2: use both forward and reverse data
%     {default = 2}
%   'optimize_B0uw_flag': [0|1] search for optimal B0 unwarp parameters
%     kernelWidthMax and lambda2
%     {default = 0}
%   'censor_min_ndirs': minimum number of diffusion directions (not including
%     b=0 images) required for tensor fit after censoring
%     will not do censoring if remaining ndirs < min_ndirs
%     {default = 6}
%   'censor_thresh': error threshold for censoring bad frames
%     normalized to median error for each slice
%     higher values mean less censoring
%     {default=3.2}
%   'driftcorr': [0|1] estimate drift correction with eddy current correction
%     {default = 0}
%   'motion_B0uw_iters': number of iterations to estimate motion
%     and B0 displacement
%     {default = 2}
%   'min_trans': minimum translation (mm) for estimating motion
%     {default = 0.05}
%   'min_rot': minimum rotation (degrees) for estimating motion
%     {default = 0.05}
%   'regT1flag': [0,1,2] whether to register diffusion data to T1
%     requires T1ContainerPath or FSContainerPath
%     0=do not register, 1=register only, 2=apply registration to diffusion data
%     {default = 2}
%   'bbregflag': [0|1] use FreeSurfer's bbregister (irrelevant if regT1flag = 0)
%     {default = 0}
%   'T1type': which T1 series to reg the DTI to
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%   'T1ContainerPath': full path of MRIPROC Container for structural
%      Used if regT1flag ~= 0
%     {default = []}
%   'FSContainerPath': full path of directory containing freesurfer recon
%     Used if bbregflag==1 instead of MPR/hiFA in ContainerPath
%        or if regT1flag = 1 and T1ContainerPath is empty
%     {default = []}
%   'outfix': string added to DTI file name after correction
%     {default = 'corr'}
%   'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%   'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default = 'T1_Atlas/T1_atlas'}
%   'forceflag': [0|1] overwrite existing files
%     {default = 0}
%
% Output:
%   infix: string applied to end of output file names
%   errcode: 0 if successful, 1 if error
%
% Created:  01/28/09 by Don Hagler
% Last Mod: 07/07/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infix = []; errcode = [];

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ContainerPath,varargin);

% get DTI scan info, determine valid scans, reference scans
[parms,errcode] = get_sess_info(parms);
if errcode, return; end;

% correct reference scan for distortions, register to T1, resample
if parms.resample_flag
  [parms,errcode] = process_reference(parms);
  if errcode, return; end;
end;

% correct each scan for distortions, correct for motion, resample
[parms,errcode] = process_data(parms);
if errcode, return; end;

infix = parms.outfix;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ContainerPath,options)
  parms_filter = {...
    'ContainerPath',ContainerPath,[],...
  ... % control flags
    'B0unwarp_flag',1,[0:2],...
    'ecc_flag',true,[false true],...
    'censor_flag',true,[false,true],...
    'motion_B0uw_flag',true,[false true],...
    'mc_flag',true,[false true],...
    'gradunwarp_flag',1,[0:2],...
    'regT1flag',2,[0 1 2],...
    'resample_flag',true,[false true],...
  ... % for resampling
    'native_flag',false,[false true],...
    'deoblique_flag',false,[false true],...
    'nvoxels',[120,120,70],[],...
    'resolution',[2 2 2],[],...
    'std_orient',[],[],...
    'rot',[0 0 0],[],...
    'trans',[0 0 0],[],...
    'smooth',0,[0,100],...
  ... % other
    'revflag',2,[0,1,2],...
    'optimize_B0uw_flag',false,[false true],...
    'inorm_B0uw_flag',true,[false true],...
    'bbregflag',false,[false true],...
    'T1type',2,[0,1,2,3],... 
    'T1ContainerPath',[],[],...
    'FSContainerPath',[],[],...
    'censor_min_ndirs',6,[],...
    'censor_thresh',3.2,[],...
    'driftcorr',false,[false true],...
    'motion_B0uw_iters',2,[1:10],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'outfix','corr',[],...
    'cleanupflag',true,[false true],...
    'forceflag',false,[false true],...
  ... % B0uw optimization
    'kernelWidthMax',25,[1:100],...
    'lambda2',1100,[1:10000],...
    'kernelWidthMax_vec',[25,31,35],[1:100],...
    'lambda2_vec',[1100,1500,1900],[1:10000],...
    'multi_opt_flag',false,[false true],...
  ... % hidden
    'fnamestem','DTI',[],...
    'ext','.mgz',{'.mgh','.mgz'},...
    'outext','.mgz',{'.mgh','.mgz'},...
		'vox2ras_maxdiff',1e-5',[1e-16,1],...
  ... % T1 atlas
    'atlasdir',[],[],...
    'atlasname','T1_Atlas/T1_atlas',[],... 
  ...
    'info_tags',{'revflag'},[],...
    'correct_tags',{'B0unwarp_flag','ecc_flag','censor_flag',...
                    'motion_B0uw_flag','mc_flag','gradunwarp_flag',...
                    'fname_out','fname_qmat','qmat','bvals','fname_censor',...
                    'censor_min_ndirs','censor_thresh','interpm',...
                    'maskoutput','verbose','forceflag',...
                    'optimize_B0uw_flag','inorm_B0uw_flag','fname_B0dx',...
                    'fname_for','fname_rev','fname_B0uw_ref',...
                    'fname_B0uw_mask','fname_B0uw_reg','PhaseDir','revflag',...
                    'driftcorr','motion_B0uw_iters','fname_ref','fname_mask',...
                    'fname_reg','min_trans','min_rot','fname_motion_B0uw',...
                    'gruw_type','gruw_unwarpflag','gruw_isoctrflag',...
                    'mstep','scales','ecc_censor_niter',...
                    'mc_censor_niter','gruw_jacobian_flag','kernelWidthMax',...
                    'lambda2','kernelWidthMax_vec','lambda2_vec',...
                    'multi_opt_flag','outfix','outext','tmpext',...
                    'cleanupflag'},[],...
    'regT1_tags',{'fname_T1','T1ContainerPath','T1type','FSContainerPath',...
                  'fname_reg','bbregflag','T2_type','interpm','ext',...
                  'forceflag','atlasdir','atlasname','cleanup_flag',...
                  'mask_smooth'},[],...
    'resample_tags',{'fname_out','save_reg_flag','fname_reg','qmat',...
                     'fname_qmat','native_flag','nvoxels','resolution',...
                     'deoblique_flag','std_orient','M_reg','M_ref',...
                     'volsz_ref','regT1flag','fname_T1','M_T1_to_EPI',...
                     'smooth','rot','trans','EPI_type','ext','forceflag',...
                     'interpm','smf'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = get_sess_info(parms)
  errcode = 0;
  args = mmil_parms2args(parms,parms.info_tags);
  [parms.ScanInfo,parms.SessInfo,errcode] = ...
    DTI_MMIL_Get_ScanInfo(parms.ContainerPath,args{:});
  if errcode~=0 || isempty(parms.ScanInfo), errcode = 1; return; end;
  % check for files required for B0 unwarp
  if parms.B0unwarp_flag
    if isempty(parms.SessInfo.B0uw_refs_for) ||...
       isempty(parms.SessInfo.B0uw_refs_rev)
      fprintf('%s: WARNING: missing scans required for B0unwarp in %s\n',...
        mfilename,parms.ContainerPath);
      if parms.B0unwarp_flag==2
        errcode = 1;
      end;
      parms.B0unwarp_flag = 0;
    end;
  end;
  fprintf('%s: processing DTI data for %s...\n',mfilename,parms.ContainerPath);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = process_reference(parms)
  errcode = 0;
  % correct reference scan for distortions
  [parms,errcode] = correct_reference(parms);
  if errcode, return; end;
  % register reference scan to T1-weighted image
  if parms.regT1flag
    [parms,errcode] = register_reference_to_T1(parms);
    if errcode, return; end;
  else
    parms.M_T1_to_DTI = [];
    parms.RegInfo = [];
  end;
  % resample reference to final space and resolution
  [parms,errcode] = resample_reference(parms);
  if errcode, return; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = correct_reference(parms)
  errcode = 0;
  tparms = parms;

  % get info about forward scan
  s_for = parms.SessInfo.reg_ref_for;
  [fname_for,M_for,volsz_for] = get_scan_info(parms,s_for,0);
  
  % get info about reverse scan
  s_rev = parms.SessInfo.reg_ref_rev;
  [fname_rev,M_rev,volsz_rev] = get_scan_info(parms,s_rev,1);

  % check that forward and reverse scans match
  if parms.B0unwarp_flag
    % check vox2ras matrices
    if any(M_for(:)-M_rev(:)>parms.vox2ras_maxdiff)
      fprintf('%s: WARNING: vox2ras matrices do not match for %s and %s\n',...
        mfilename,fname_for,fname_rev);
      errcode = 1;
      return;
    end;
    % check dimensions
    if any(volsz_for(1:3)~=volsz_rev(1:3))
      fprintf('%s: WARNING: volume dimensions do not match for %s and %s\n',...
        mfilename,fname_for,fname_rev);
      errcode = 1;
      return;
    end;
    % check PhaseDir
    if ~strcmp(parms.ScanInfo(s_for).PhaseDir,parms.ScanInfo(s_rev).PhaseDir)
      fprintf('%s: WARNING: phase-encoding direction does not match for %s and %s\n',...
        mfilename,fname_for,fname_rev);
      errcode = 1;
      return;
    end;
  end;

  % correct forward reference scan for distortions
  if ~isempty(fname_for)
    [tparms,fname_for,errcode] = set_correct_parms(parms,s_for,0,1);
    if errcode, return; end;
    [parms.fname_ref_for,errcode] = correct_scan(fname_for,tparms);
    if errcode, return; end;
  else
    parms.fname_ref_for = [];
  end;

  % correct reverse reference scan for distortions
  if ~isempty(fname_rev)
    [tparms,fname_rev,errcode] = set_correct_parms(parms,s_rev,1,1);
    if errcode, return; end;
    [parms.fname_ref_rev,errcode] = correct_scan(fname_rev,tparms);
    if errcode, return; end;
  else
    parms.fname_ref_rev = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname,fstem] = set_scan_fname(parms,s,revflag,infix)
  if ~exist('infix','var'), infix = []; end;
  fname = []; fstem = [];
  if ~isempty(s)
    fstem = sprintf('%s/%s%d',parms.ContainerPath,parms.fnamestem,s);
    if revflag, fstem = [fstem '_rev']; end;
    if ~isempty(infix), fstem = [fstem '_' infix]; end;
    fname = [fstem parms.ext];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname,M,volsz] = get_scan_info(parms,s,revflag,outdir)
  if ~exist('outdir','var') || isempty(outdir)
    outdir = parms.ContainerPath;
  end;
  fname = []; M = []; volsz = [];
  if ~isempty(s)
    fname = set_scan_fname(parms,s,revflag);
    [M,volsz] = mmil_load_mgh_info(fname,parms.forceflag,outdir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = copy_first_frame(fname_in,parms)
  [fpath,fstem,fext] = fileparts(fname_in);
  fname_out = sprintf('%s/%s_f0%s',fpath,fstem,fext);
  if ~exist(fname_out,'file') || parms.forceflag
    [vol,M]=fs_load_mgh(fname_in,[],1);
    fs_save_mgh(vol,fname_out,M);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,fname,errcode] = set_correct_parms(parms,s,revflag,refflag)
  fname = [];
  errcode = 0;
  if ~exist('refflag','var') || isempty(refflag), refflag = 0; end;
  parms.revflag = revflag;

  % set B0uw_ref file names
  sA = parms.ScanInfo(s).B0uw_ref_for;
  sB = parms.ScanInfo(s).B0uw_ref_rev;
  if ~isempty(sA) && ~isempty(sB)
    parms.fname_for = sprintf('%s/%s%d%s',...
      parms.ContainerPath,parms.fnamestem,sA,parms.ext);
    parms.fname_rev = sprintf('%s/%s%d_rev%s',...
      parms.ContainerPath,parms.fnamestem,sB,parms.ext);
    parms.PhaseDir = parms.ScanInfo(s).PhaseDir;
  else
    parms.fname_for = [];
    parms.fname_rev = [];
    parms.PhaseDir = [];
  end;

  % get gradwarp parameters from parms.ScanInfo
  if isfield(parms.ScanInfo(s).gradwarpinfo,'gwtype')
    parms.gruw_type = parms.ScanInfo(s).gradwarpinfo.gwtype;
  else
    parms.gruw_type = 0;
    if parms.gradunwarp_flag
      if isfield(parms.ScanInfo(s).gradwarpinfo,'ambiguousgwtype') &...
         parms.ScanInfo(s).gradwarpinfo.ambiguousgwtype
        fprintf('%s: WARNING: ambiguous gwtype... skipping grad unwarp\n',...
          mfilename);
      elseif ~mmil_Philips(parms.SessInfo)
        fprintf('%s: WARNING: missing gwtype... skipping grad unwarp\n',...
          mfilename);
      end;
      if parms.gradunwarp_flag==2
        errcode = 1;
      end;
      parms.gradunwarp_flag = 0;
    end;
  end;

  % set unwarpflag to 0 if no unwarpflag exists
  %  (so it will work with Philips, which has no gradwarp)
  parms.gruw_unwarpflag = ...
    mmil_getfield(parms.ScanInfo(s).gradwarpinfo,'unwarpflag',0);
  parms.gruw_isoctrflag = ...
    mmil_getfield(parms.ScanInfo(s).gradwarpinfo,'isoctrflag',0);

  % set parameters differently for reference scan vs others
  if refflag
    % save first frame of forward scan if multi-frame
    [fname,M,volsz] = get_scan_info(parms,s,revflag);
    if volsz(4)>1, fname = copy_first_frame(fname,parms); end;
    % no eddy current or motion correction for single-frame reference
    parms.ecc_flag = 0;
    parms.mc_flag = 0;
    parms.fname_B0uw_ref = fname;    
    parms.motion_B0uw_iters = 1;
  else
    % set input file name for scan
    [fname,fstem] = set_scan_fname(parms,s,revflag);
    % check for frame-censoring file name
    parms.fname_censor = [fstem '_censor.txt'];
    if ~exist(parms.fname_censor,'file')
      parms.fname_censor = [];
    end;
    % matrix of diffusion vectors (ndirs x 3)
    parms.qmat = set_qmat(parms.ScanInfo(s),revflag);
    parms.bvals = set_bvals(parms.ScanInfo(s),parms.qmat);
    % set parameters for resampling to reference
    if parms.resample_flag
      % set outfix to include '_regT1' or '_resDTI' depending on regT1flag
      if parms.regT1flag==2
        parms.outfix = [parms.outfix '_regT1'];
      else
        parms.outfix = [parms.outfix '_resDTI'];
      end;
      % set file names for between-scan motion correction
      if ~revflag
        s_ref = parms.SessInfo.reg_ref_for;
        parms.fname_ref = parms.fname_ref_for;
      else
        s_ref = parms.SessInfo.reg_ref_rev;
        parms.fname_ref = parms.fname_ref_rev;
      end;
      [fpath,fstem_ref] = fileparts(parms.fname_ref);
      fstem_ref = [fpath '/' fstem_ref];
      parms.fname_mask = sprintf('%s_mask%s',fstem_ref,parms.ext);
      parms.fname_reg = sprintf('%s_reg%d_%s.mat',fstem,s_ref,parms.outfix);
    end;
    % set reference for registering to scans used for estimating B0 distortion
    if ~revflag
      s_ref = parms.ScanInfo(s).B0uw_ref_for;
    else
      s_ref = parms.ScanInfo(s).B0uw_ref_rev;
    end;
    [parms.fname_B0uw_ref,fstem_B0uw_ref] = ...
      set_scan_fname(parms,s_ref,revflag);
    parms.fname_B0uw_mask = sprintf('%s_mask%s',fstem_B0uw_ref,parms.ext);
    parms.fname_B0uw_reg = sprintf('%s_reg%d.mat',fstem,s_ref);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,errcode] = correct_scan(fname,parms)
  errcode = 0;
  fname_out = [];
  args = mmil_parms2args(parms,parms.correct_tags);
  try
    [fname_out,fname_B0dx,fname_qmat] = dti_correct_data(fname,args{:});
  catch me
    fprintf('%s: WARNING: dti_correct_data failed for %s:\n%s\n',...
      mfilename,fname,me.message);
    errcode = 1;
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = register_reference_to_T1(parms)
  errcode = 0;
  if ~parms.SessInfo.revflag
    fname = parms.fname_ref_for;
  else
    fname = parms.fname_ref_rev;
  end;
  [fpath,fstem_ref] = fileparts(fname);
  tparms = parms;
  tparms.T2_type = ['T2_' parms.fnamestem];
  if ~parms.bbregflag && ~isempty(parms.T1ContainerPath)
    tparms.FSContainerPath = [];
  end;
  args = mmil_parms2args(tparms,parms.regT1_tags);
  [parms.M_T1_to_DTI,parms.RegInfo,errcode] = ...
    MMIL_Register_T2_to_T1(fname,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = resample_reference(parms)
  errcode = 0;
  % resample forward reference scan
  if ~isempty(parms.fname_ref_for)
    [tparms,fname] = set_resample_parms(parms,0);
    [parms.fname_ref_for,errcode] = resample_scan(fname,tparms);
    if errcode, return; end;
  end;
  % resample reverse reference scan
  if ~isempty(parms.fname_ref_rev)
    [tparms,fname] = set_resample_parms(parms,1);
    [parms.fname_ref_rev,errcode] = resample_scan(fname,tparms);
    if errcode, return; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,fname] = set_resample_parms(parms,revflag)
  if ~revflag
    fname = parms.fname_ref_for;
  else
    fname = parms.fname_ref_rev;
  end;
  [fpath,fstem,fext] = fileparts(fname);
  if parms.regT1flag==2
    outfix = 'regT1';
  else
    outfix = 'resDTI';
  end;
  % set fname_out
  parms.fname_out = [fpath '/' fstem '_' outfix parms.ext];
  % set fname_reg
  if parms.regT1flag>0
    parms.fname_reg = [fpath '/' fstem '_' outfix '_regT1.mat'];
  else
    parms.fname_reg = [];
  end;
  % set fname_T1
  if isempty(parms.RegInfo)
    parms.fname_T1 = [];
  else
    parms.fname_T1 = parms.RegInfo.fname_T1;
    if ~exist(parms.fname_T1)
      % look for file with same fstem in T1ContainerPath
      [fpath_T1,fstem_T1,fext_T1] = fileparts(parms.fname_T1);
      tmp_fname_T1 = [parms.T1ContainerPath '/' fstem_T1 fext_T1];
      if exist(tmp_fname_T1,'file')
        parms.fname_T1 = tmp_fname_T1;
      else
        fprintf('%s: WARNING: T1 file %s not found\n',...
          mfilename,parms.fname_T1);
      end;
    end;
  end;
  % set additional parameters
  parms.save_reg_flag = 1;
  parms.regT1flag = (parms.regT1flag==2);
  parms.M_T1_to_EPI = parms.M_T1_to_DTI;
  parms.EPI_type = parms.fnamestem;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,errcode] = resample_scan(fname,parms)
  errcode = 0;
  fname_out = [];
  args = mmil_parms2args(parms,parms.resample_tags);
  try
    epi_resample_data(fname,args{:});
  catch me
    fprintf('%s: WARNING: epi_resample_data failed for %s:\n%s\n',...
      mfilename,fname,me.message);
    errcode = 1;
    return;
  end;
  fname_out = parms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = process_data(parms)
  errcode = 0;
  for s=parms.SessInfo.snums_valid
    revflag_list = set_revflag_list(parms.ScanInfo(s),parms.revflag);
    for revflag = revflag_list
      [tparms,fname,errcode] = set_correct_parms(parms,s,revflag);
      if errcode, continue; end;
      [fname_out,errcode] = correct_scan(fname,tparms);
    end;
  end;
  parms.outfix = tparms.outfix;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function revflag_list = set_revflag_list(ScanInfo,revflag)
  revflag_list = [];
  if ismember(ScanInfo.pepolar,[0,1])
    revflag_list = ScanInfo.pepolar;
  elseif revflag==2
    revflag_list = [0,1];
  else
    revflag_list = revflag;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qmat = set_qmat(ScanInfo,revflag)
  if ((ScanInfo.pepolar == 2 && revflag == 1) || ...
      (ScanInfo.pepolar == 3 && revflag == 0)) && ...
      ~ismember(ScanInfo.DTI_Sequence_Type,[1:3])
    qmat = [];
  else
    qmat = ScanInfo.qmat;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bvals = set_bvals(ScanInfo,qmat)
  if isempty(qmat)
    bvals = ScanInfo.bval;
  elseif ScanInfo.DTI_Sequence_Type==6
    bvals = ScanInfo.bval*sum(qmat.^2,2);
  else
    bvals = ScanInfo.bval*ones(size(qmat,1),1);
  end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

