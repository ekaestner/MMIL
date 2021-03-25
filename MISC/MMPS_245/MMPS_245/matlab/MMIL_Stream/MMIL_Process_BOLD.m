function [infix,errcode]=MMIL_Process_BOLD(ContainerPath,varargin)
%function [infix,errcode]=MMIL_Process_BOLD(ContainerPath,varargin)
%
% Required Input:
%   ContainerPath: full path of BOLDPROC Container Directory
%
% Optional Control Parameters:
%   'snums' - vector of scan numbers to process (if empty, do all)
%     {default = []}
%   'outfix': string attached to BOLD file names after correction
%     {default = 'corr'}
%   'B0unwarp_flag': [0|1|2] whether to correct for B0 distortions
%     if 2, processing aborted if B0 unwarping fails (e.g. no rev scans)
%     {default = 1}
%   'optimize_B0uw_flag': [0|1] search for optimal B0 unwarp parameters
%     kernelWidthMax and lambda2
%     {default = 0}
%   'regref_B0uw_flag': [0|1] whether to register reference image (fname_ref)
%      to fname_for (if revflag = 0) or fname_rev (if revflag = 1)
%      and apply transformation to B0dx
%     {default = 1}
%   'motion_B0uw_flag': [0|1] whether to estimate head motion
%     and apply to B0dx field
%     {default = 0}
%   'tshift_flag': [0|1] correct for differences in slice timing
%     {default = 1}
%   'mc_flag' - [0|1] correct for within scan motion
%     {default = 1}
%   'gradunwarp_flag': [0|1|2] whether to correct for gradient distortions
%     {default = 1}
%   'regT1flag': [0,1,2] whether to register BOLD data to T1
%     requires T1ContainerPath or FSContainerPath
%     0=do not register, 1=register only, 2=apply registration to BOLD data
%     {default = 1}
%
% Other Optional Parameters:
%   'motion_B0uw_iters': number of iterations to estimate motion
%     and B0 displacement
%     {default = 2}
%   'min_trans': minimum translation (mm) for estimating motion for B0uw
%     {default = 0.05}
%   'min_rot': minimum rotation (degrees) for estimating motion for B0uw
%     {default = 0.05}
%   'tpattern': slice time pattern
%     Allowed values: {'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z'}
%     {default = 'alt+z'}
%   'skipTRs': number of TRs at beginning of scan to be ignored in time shifting
%     {default = 0}
%   'bbegflag': [0|1] use FreeSurfer's bbregister (irrelevant if regT1flag = 0)
%      {default = 0}
%   'T1type': which T1 series to reg the BOLD data to
%      0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%      {default=2}
%   'T1ContainerPath': full path of MRIPROC Container for structural
%      used if regT1flag ~= 0
%      {default = []}
%   'FSContainerPath': full path of directory containing freesurfer recon
%      used if bbregflag=1 or if regT1flag = 1 and T1ContainerPath is empty
%      {default = []}
%   'stdev_flag': [0|1] calculate inter- and intra-scan variability
%      {default = 1}
%   'stdev_plotflag': [0|1] plot images of inter- and intra-scan variability
%      irrelevant if stdev_flag = 0
%      {default = 1}
%   'cleanupflag': [0|1] remove intermediate files after processing
%      {default = 1}
%   'verbose': [0|1] display status messages
%      {default = 1}
%   'forceflag': [0|1] whether to overwrite existing output
%      {default = 0}
%
% Optional Parameters for Resampling:
%   'resample_flag': [0|1] whether to resample BOLD data
%     e.g. to make isotropic, to correct inter-scan motion, to register to T1
%     {default = 1}
%   'native_flag': resample but keep original nvoxels and resolution
%     {default = 0}
%   'nvoxels': vector of number of resampled voxels [nx,ny,nz]
%     ignored if native_flag = 1
%     {default = [120 120 70]}
%   'resolution': desired voxel sizes of resampled data (x y z)
%     ignored if native_flag = 1
%     {default = [2 2 2]}
%   'deoblique_flag': [0|1] whether to resample oblique slices to on-axis
%     ignored if native_flag = 1
%     {default = 1}
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
% Output:
%   infix: string applied to end of output file names
%   errcode: 0 if successful, 1 if error
%
% Created:  06/25/08 by Don Hagler
% Prev Mod: 07/24/17 by Don Hagler
% Last Mod: 08/15/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infix = []; errcode = [];

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ContainerPath,varargin);

% get BOLD scan info, determine valid scans, reference scans
[parms,errcode] = get_sess_info(parms);
if errcode, return; end;

% correct reference scan for distortions, register to T1, resample
if parms.resample_flag
  [parms,errcode] = process_reference(parms);
  if errcode, return; end;
end;

% correct each scan for distortions, slice timing, and motion
[parms,errcode] = process_data(parms);
if errcode, return; end;

% for BOLD_ape scans, recombine for and rev, also motion.1D files
[parms,errcode] = recombine_ape_scans(parms);
if errcode, return; end;

% calculate inter- and intra-scan variability
if parms.stdev_flag
  [parms,errcode] = calc_stdev(parms);
  if errcode, return; end;
end;

infix = parms.outfix;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ContainerPath,options)
  parms_filter = {...
    'ContainerPath',ContainerPath,[],...
  ...
    'snums',[],[],...
    'outfix','corr',[],...
    'B0unwarp_flag',1,[0:2],...
    'optimize_B0uw_flag',false,[false true],...
    'regref_B0uw_flag',true,[false true],...
    'mask_thresh',0.82,[0,1],...
    'motion_B0uw_flag',false,[false true],...
    'tshift_flag',true,[false true],...
    'mc_flag',true,[false true],...
    'gradunwarp_flag',1,[0:2],...
    'regT1flag',1,[0 1 2],...
    'motion_B0uw_iters',2,[1:10],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'tpattern','alt+z',{'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z'},...
    'skipTRs',0,[0,Inf],...
    'bbregflag',false,[false true],...
    'T1type',2,[0,1,2,3],... 
    'T1ContainerPath',[],[],...
    'FSContainerPath',[],[],...
    'stdev_flag',true,[false true],...
    'stdev_plotflag',true,[false true],...
    'cleanupflag',true,[false true],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % for resampling
    'resample_flag',true,[false true],...
    'native_flag',false,[false true],...
    'nvoxels',[120,120,70],[],...
    'resolution',[2 2 2],[],...
    'deoblique_flag',true,[false true],...
    'rot',[0 0 0],[],...
    'trans',[0 0 0],[],...
    'smooth',0,[0,100],...
  ... % B0uw optimization
    'kernelWidthMax',25,[1:100],...
    'lambda2',1100,[1:10000],...
    'kernelWidthMax_vec',[25,31,35],[1:100],...
    'lambda2_vec',[1100,1500,1900],[1:10000],...
    'multi_opt_flag',false,[false true],...
  ... % hidden
    'fnamestem','BOLD',[],...
    'ext','.mgz',{'.mgh','.mgz'},...
    'outext','.mgz',{'.mgh','.mgz'},...
  ...
    'info_tags',{'snums','revflag'},[],...
    'correct_tags',{'B0unwarp_flag','regref_B0uw_flag','mask_thresh',...
                    'motion_B0uw_flag','tshift_flag',...
                    'mc_flag','gradunwarp_flag','fname_out','verbose',...
                    'forceflag','optimize_B0uw_flag',...
                    'fname_B0dx','fname_for','fname_rev',...
                    'fname_B0uw_ref','fname_B0uw_mask','revflag',...
                    'EchoSpacing','AcquisitionColumns','AcquisitionRows',...
                    'PhaseDir','EchoSpacing_for','AcquisitionColumns_for',...
                    'AcquisitionRows_for','EchoSpacing_rev',...
                    'AcquisitionColumns_rev','AcquisitionRows_rev',...
                    'motion_B0uw_iters','min_trans','min_rot',...
                    'fname_motion_B0uw','fname_regref_B0uw',...
                    'tpattern','skipTRs','ts_suffix',...
                    'fname_ref','ref_frame','fname_motion_mat',...
                    'fname_motion_1D','mc_suffix','gruw_type',...
                    'gruw_unwarpflag','gruw_isoctrflag','interpm',...
                    'gruw_jacobian_flag','kernelWidthMax','lambda2',...
                    'kernelWidthMax_vec','lambda2_vec','multi_opt_flag',...
                    'outfix','outext','tmpext','cleanupflag'},[],...
    'regT1_tags',{'fname_T1','T1ContainerPath','T1type','FSContainerPath',...
                  'fname_reg','bbregflag','T2_type','interpm','ext',...
                  'forceflag','atlasdir','atlasname','cleanup_flag',...
                  'mask_smooth'},[],...
    'resample_tags',{'fname_out','save_reg_flag','fname_reg',...
                     'native_flag','nvoxels','resolution',...
                     'deoblique_flag','std_orient','M_reg','M_ref',...
                     'volsz_ref','regT1flag','fname_T1','M_T1_to_EPI',...
                     'smooth','rot','trans','EPI_type','ext','forceflag',...
                     'interpm','smf'},[],...
     'merge_tags',{'pepolar','mc_flag','infix','fname_out',...
                   'fname_motion_out','normflag','verbose','forceflag'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = get_sess_info(parms)
  errcode = 0;

  % get scan and sess info
  args = mmil_parms2args(parms,parms.info_tags);
  [parms.ScanInfo,parms.SessInfo,errcode] = ...
    BOLD_MMIL_Get_ScanInfo(parms.ContainerPath,args{:});
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

  % check for valid scans
  ind_valid = parms.SessInfo.snums_valid;
  if isempty(ind_valid)
    errcode = 1;
    return;
  end;
  
  fprintf('%s: processing BOLD data for %s...\n',mfilename,parms.ContainerPath);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = process_reference(parms)
  errcode = 0;
  % correct reference scan for distortions
  [parms,errcode] = correct_reference(parms);
  if errcode, return; end;
  % register SE reference to GE reference
  if parms.regref_B0uw_flag
    [parms,errcode] = register_reference(parms);
    if errcode, return; end;
  end;
  % register reference scan to T1-weighted image
  if parms.regT1flag
    [parms,errcode] = register_reference_to_T1(parms);
    if errcode, return; end;
  else
    parms.M_T1_to_BOLD = [];
    parms.RegInfo = [];
  end;
  % resample reference to final space and resolution
  [parms,errcode] = resample_reference(parms);
  if errcode, return; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = correct_reference(parms)
  errcode = 0;

  % correct regT1 forward reference scan for distortions
  parms.fname_ref_regT1 = [];
  s_for = parms.SessInfo.regT1_ref_for;
  if ~isempty(s_for)
    [tparms,fname_ref_regT1_for,errcode] = set_correct_parms(parms,s_for,0,1);
    if errcode, return; end;
    [parms.fname_ref_regT1_for,errcode] = correct_scan(fname_ref_regT1_for,tparms);
    if errcode, return; end;
  else
    parms.fname_ref_regT1_for = [];
  end;

  % correct regT1 reverse reference scan for distortions
  s_rev = parms.SessInfo.regT1_ref_rev;
  if ~isempty(s_rev)
    [tparms,fname_ref_regT1_rev,errcode] = set_correct_parms(parms,s_rev,1,1);
    if errcode, return; end;
    [parms.fname_ref_regT1_rev,errcode] = correct_scan(fname_ref_regT1_rev,tparms);
    if errcode, return; end;
  else
    parms.fname_ref_regT1_rev = [];
  end;

  if parms.SessInfo.revflag
    parms.fname_ref_regT1 = parms.fname_ref_regT1_rev;
  else
    parms.fname_ref_regT1 = parms.fname_ref_regT1_for;
  end;

  % correct forward reference scan for distortions
  s_for = parms.SessInfo.reg_ref_for;
  if ~isempty(s_for)
    [tparms,fname_ref_for,errcode] = set_correct_parms(parms,s_for,0,1);
    if errcode, return; end;
    [parms.fname_ref_for,errcode] = correct_scan(fname_ref_for,tparms);
    if errcode, return; end;
  else
    parms.fname_ref_for = [];
  end;
  
  % correct reverse reference scan for distortions
  s_rev = parms.SessInfo.reg_ref_rev;
  if ~isempty(s_rev)
    [tparms,fname_ref_rev,errcode] = set_correct_parms(parms,s_rev,1,1);
    if errcode, return; end;
    [parms.fname_ref_rev,errcode] = correct_scan(fname_ref_rev,tparms);
    if errcode, return; end;
  else
    parms.fname_ref_rev = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = register_reference(parms)
  errcode = 0;
  parms.M_GE_to_SE = [];
  % set file names of volumes to be registered
  revflag = parms.SessInfo.revflag;
  if revflag
    fname_GE = parms.fname_ref_rev;
    fname_SE = parms.fname_ref_regT1_rev;
  else
    fname_GE = parms.fname_ref_for;
    fname_SE = parms.fname_ref_regT1_for;
  end;
  [fpath,fstem_SE] = fileparts(fname_SE);
  [fpath,fstem_GE] = fileparts(fname_GE);
  if ~isempty(regexp(fname_GE,fname_SE))
    if parms.verbose
      fprintf('%s: skipping registration between %s and %s...\n',mfilename,fstem_SE,fstem_GE);
    end;
    return;
  end;
  if parms.verbose
    fprintf('%s: registering reference images %s and %s...\n',...
      mfilename,fstem_SE,fstem_GE);
  end;
  fname_regref_B0uw = sprintf('%s/%s_regref_B0uw.mat',...
    parms.ContainerPath,fstem_GE);
  if ~exist(fname_regref_B0uw,'file') || parms.forceflag
    tparms = [];
    tparms.outdir = sprintf('%s/tmp_BOLD_corr/reg_GE_%s_SE_%s',...
      parms.ContainerPath,fstem_GE,fstem_SE);
    tparms.mask_thresh = parms.mask_thresh;
    tparms.cleanup_flag = 0;
    tparms.forceflag = 0;
    args = mmil_parms2args(tparms);
    M_GE_to_SE = mmil_jpdfreg_GESE(fname_GE,fname_SE,args{:});
    save(fname_regref_B0uw,'M_GE_to_SE','fname_SE','fname_GE');
  else
    tmp = load(fname_regref_B0uw);
    M_GE_to_SE = tmp.M_GE_to_SE;
  end;
  parms.M_GE_to_SE = M_GE_to_SE;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname,fstem] = set_scan_fname(parms,s,revflag,infix)
  if ~exist('infix','var'), infix = []; end;
  fname = []; fstem = [];
  if ~isempty(s)
    fstem = sprintf('%s/%s%d',parms.ContainerPath,parms.fnamestem,s);
    if revflag
      fstem = [fstem '_rev'];    
    else
      fstem = [fstem '_for'];
    end;
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
  sAreg = parms.ScanInfo(s).B0reg_ref_for;
  sBreg = parms.ScanInfo(s).B0reg_ref_rev;
  if parms.B0unwarp_flag && ~isempty(sA) && ~isempty(sB)
    parms.fname_for = sprintf('%s/%s%d_for%s',...
      parms.ContainerPath,parms.fnamestem,sA,parms.ext);
    parms.fname_rev = sprintf('%s/%s%d_rev%s',...
      parms.ContainerPath,parms.fnamestem,sB,parms.ext);
    parms.EchoSpacing_for = single(parms.ScanInfo(sA).EchoSpacing);
    parms.AcquisitionColumns_for = single(parms.ScanInfo(sA).AcquisitionColumns);
    parms.AcquisitionRows_for = single(parms.ScanInfo(sA).AcquisitionRows);
    parms.EchoSpacing_rev = single(parms.ScanInfo(sB).EchoSpacing);
    parms.AcquisitionColumns_rev = single(parms.ScanInfo(sB).AcquisitionColumns);
    parms.AcquisitionRows_rev = single(parms.ScanInfo(sB).AcquisitionRows);
    parms.EchoSpacing = single(parms.ScanInfo(s).EchoSpacing);
    parms.AcquisitionColumns = single(parms.ScanInfo(s).AcquisitionColumns);
    parms.AcquisitionRows = single(parms.ScanInfo(s).AcquisitionRows);
    parms.PhaseDir = parms.ScanInfo(s).PhaseDir;
  else
    parms.fname_for = [];
    parms.fname_rev = [];
    parms.EchoSpacing_for = [];
    parms.AcquisitionColumns_for = [];
    parms.AcquisitionRows_for = [];
    parms.EchoSpacing_rev = [];
    parms.AcquisitionColumns_rev = [];
    parms.AcquisitionRows_rev = [];
    parms.EchoSpacing = [];
    parms.AcquisitionColumns = [];
    parms.AcquisitionRows = [];
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
    % no slice timing or motion correction for single-frame reference
    parms.tshift_flag = 0;
    parms.mc_flag = 0;
    parms.fname_B0uw_ref = fname;    
    parms.motion_B0uw_iters = 1;
  else
    % set input file name for scan
    [fname,fstem] = set_scan_fname(parms,s,revflag);
    % slice timing pattern
    tpattern = mmil_getfield(parms.ScanInfo(s),'slice_tpattern');
    if isempty(tpattern)
      if parms.verbose && parms.tshift_flag
        fprintf('%s: WARNING: using default slice tpattern %s\n',...
          mfilename,tpattern);
      end;
      tpattern = parms.tpattern;
    end;
    parms.tpattern = tpattern;
    % set parameters for resampling to reference
    if parms.resample_flag
      % set outfix to include '_regT1' or '_resBOLD' depending on regT1flag
      if parms.regT1flag==2
        parms.outfix = [parms.outfix '_regT1'];
      else
        parms.outfix = [parms.outfix '_resBOLD'];
      end;
      % set file names for between-scan motion correction
      if ismember(s,parms.SessInfo.snums_SEpep)
        if ~revflag
          s_ref = parms.SessInfo.regT1_ref_for;
          parms.fname_ref = parms.fname_ref_regT1_for;
        else
          s_ref = parms.SessInfo.regT1_ref_rev;
          parms.fname_ref = parms.fname_ref_regT1_rev;
        end;
      else
        if ~revflag
          s_ref = parms.SessInfo.reg_ref_for;
          parms.fname_ref = parms.fname_ref_for;
        else
          s_ref = parms.SessInfo.reg_ref_rev;
          parms.fname_ref = parms.fname_ref_rev;
        end;
      end;
      if ~isempty(parms.fname_ref)
        [fpath,fstem_ref] = fileparts(parms.fname_ref);
        fstem_ref = [fpath '/' fstem_ref];
        parms.fname_mask = sprintf('%s_mask%s',fstem_ref,parms.ext);
        parms.fname_reg = sprintf('%s_reg%d.mat',fstem,s_ref);
      else
        parms.fname_mask = [];
        parms.fname_reg = [];
      end;
    end;
    % set reference for registering to scans used for estimating B0 distortion
    if ~revflag
      s_ref = parms.ScanInfo(s).B0reg_ref_for;
    else
      s_ref = parms.ScanInfo(s).B0reg_ref_rev;
    end;
    [parms.fname_B0uw_ref,fstem_B0uw_ref] = ...
      set_scan_fname(parms,s_ref,revflag);
    parms.fname_B0uw_mask = sprintf('%s_mask%s',fstem_B0uw_ref,parms.ext);
  end;
  % if scan is one of a B0uw pair, disable registration to self
  if ~revflag
    B0uw_refs = parms.SessInfo.B0uw_refs_for;
  else
    B0uw_refs = parms.SessInfo.B0uw_refs_rev;
  end;
  if ismember(s,B0uw_refs)
    parms.regref_B0uw_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,errcode] = correct_scan(fname,parms)
  errcode = 0;
  fname_out = [];
  args = mmil_parms2args(parms,parms.correct_tags);
%% todo: uncomment try/catch after debugging is finished
%  try
    [fname_out,fname_B0dx,fname_motion_1D] = bold_correct_data(fname,args{:});
%  catch
%    fprintf('%s: WARNING: bold_correct_data failed for %s:\n%s\n',...
%      mfilename,fname,lasterr);
%    errcode = 1;
%    return;
%  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = register_reference_to_T1(parms)
  errcode = 0;
  fname = parms.fname_ref_regT1;
  [fpath,fstem_ref] = fileparts(fname);
  tparms = parms;
  tparms.T2_type = ['T2_' parms.fnamestem];
  if ~parms.bbregflag && ~isempty(parms.T1ContainerPath)
    tparms.FSContainerPath = [];
  end;
  args = mmil_parms2args(tparms,parms.regT1_tags);
  [parms.M_T1_to_BOLD,parms.RegInfo,errcode] = ...
    MMIL_Register_T2_to_T1(fname,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = resample_reference(parms)
  % resample regT1 forward reference scan
  if ~isempty(parms.fname_ref_regT1_for)
    [tparms,fname] = set_resample_parms(parms,parms.fname_ref_regT1_for,1,0);
    [parms.fname_ref_regT1_for,errcode] = resample_scan(fname,tparms);
    if errcode, return; end;
  end;
  % resample regT1 reverse reference scan
  if ~isempty(parms.fname_ref_regT1_rev)
    [tparms,fname] = set_resample_parms(parms,parms.fname_ref_regT1_rev,0,0);
    [parms.fname_ref_regT1_rev,errcode] = resample_scan(fname,tparms);
    if errcode, return; end;
  end;
  if parms.SessInfo.revflag
    parms.fname_ref_regT1 = parms.fname_ref_regT1_rev;
  else
    parms.fname_ref_regT1 = parms.fname_ref_regT1_for;
  end;
  % resample forward reference scan
  if ~isempty(parms.fname_ref_for)
    [tparms,fname] = set_resample_parms(parms,...
      parms.fname_ref_for,0,parms.regref_B0uw_flag);
    [parms.fname_ref_for,errcode] = resample_scan(fname,tparms);
    if errcode, return; end;
  end;
  % resample reverse reference scan
  if ~isempty(parms.fname_ref_rev)
    [tparms,fname] = set_resample_parms(parms,...
      parms.fname_ref_rev,0,parms.regref_B0uw_flag);
    [parms.fname_ref_rev,errcode] = resample_scan(fname,tparms);
    if errcode, return; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,fname] = set_resample_parms(parms,fname,save_reg_flag,regref_flag)
  if ~exist('save_reg_flag','var') || isempty(save_reg_flag)
    save_reg_flag = 0;
  end;
  if ~exist('regref_flag','var') || isempty(regref_flag)
    regref_flag = 0;
  end;
  [fpath,fstem,fext] = fileparts(fname);
  if parms.regT1flag==2
    outfix = 'regT1';
  else
    outfix = 'resBOLD';
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
  parms.save_reg_flag = save_reg_flag;
  parms.regT1flag = (parms.regT1flag==2);
  parms.M_T1_to_EPI = parms.M_T1_to_BOLD;
  parms.EPI_type = parms.fnamestem;
  if regref_flag
    parms.M_reg = inv(parms.M_GE_to_SE);
  else
    parms.M_reg = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,errcode] = resample_scan(fname,parms)
  errcode = 0;
  fname_out = [];
  args = mmil_parms2args(parms,parms.resample_tags);
  try
    epi_resample_data(fname,args{:});
  catch
    fprintf('%s: WARNING: epi_resample_data failed for %s:\n%s\n',...
      mfilename,fname,lasterr);
    errcode = 1;
    return;
  end;
  fname_out = parms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = process_data(parms)
  errcode = 0;
  for s=1:length(parms.SessInfo.snums_valid)
    snum = parms.SessInfo.snums_valid(s);
    revflag_list = set_revflag_list(parms.ScanInfo(snum));
    for r=1:length(revflag_list)
      revflag = revflag_list(r);
      [tparms,fname,errcode] = set_correct_parms(parms,snum,revflag);
      if errcode, continue; end;
      [fname_out,errcode] = correct_scan(fname,tparms);
    end;
  end;
  parms.outfix = tparms.outfix;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function revflag_list = set_revflag_list(ScanInfo)
  revflag_list = [];
  switch ScanInfo.pepolar
    case {0,1}
      revflag_list = ScanInfo.pepolar;
    case {2,3}
      revflag_list = [0,1];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = recombine_ape_scans(parms)
  errcode = 0;
  for i=1:length(parms.SessInfo.snums_valid)
    s = parms.SessInfo.snums_valid(i);
    if strcmp(parms.ScanInfo(s).ScanType,'BOLD_ape') &&...
              parms.ScanInfo(s).pepolar>1
      fstem = sprintf('%s/%s%d',parms.ContainerPath,parms.fnamestem,s);
      fname_rev = [fstem '_rev_' parms.outfix '.mgz'];
      fname_for = [fstem '_for_' parms.outfix '.mgz'];
      tparms = parms;
      tparms.pepolar = parms.ScanInfo(s).pepolar;
      tparms.infix = parms.outfix;
      args = mmil_parms2args(tparms,parms.merge_tags);
      try
        bold_merge_APE(fname_rev,fname_for,args{:});
      catch
        fprintf('\n%s: WARNING: bold_merge_APE failed for %s:\n%s\n\n',...
          mfilename,fstem,lasterr);
        errcode = 1;
        continue;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = calc_stdev(parms)
  errcode = 0;
  snums_list = {...
    parms.SessInfo.snums_for...
    parms.SessInfo.snums_rev...
    parms.SessInfo.snums_SE_for...
    parms.SessInfo.snums_SE_rev...
    parms.SessInfo.snums_SEpep...
    parms.SessInfo.snums_SEpep...
  };
  outdir_list = {...
    'intra_stdev_for'...
    'intra_stdev_rev'...
    'inter_stdev_SE_for'...
    'inter_stdev_SE_rev'...
    'inter_stdev_SE_for'...
    'inter_stdev_SE_rev'...
  };
  intra_flags = [1 1 0 0 0 0];
  rev_flags = [0 1 0 1 0 1];
  
  tp = [];
  tp.outstem = parms.fnamestem;
  tp.plotflag = parms.stdev_plotflag;
  tp.verbose = parms.verbose;
  tp.forceflag = parms.forceflag;
  
  for i=1:length(snums_list)
    snums = snums_list{i};
    outdir = outdir_list{i};
    if isempty(snums)
      if parms.verbose
        fprintf('%s: WARNING: skipping stdev calculations for %s...\n',...
          mfilename,outdir);
      end;
      continue;
    end;
    revflag = rev_flags(i);
    intra_flag = intra_flags(i);
    fnamelist = {};
    for j=1:length(snums)
      s = snums(j);
      fnamelist{j} = set_scan_fname(parms,s,revflag,parms.outfix);
    end;
    tp.outdir = sprintf('%s/%s',parms.ContainerPath,outdir);
    args = mmil_parms2args(tp);
    if intra_flag
      results = mmil_intrascan_stdev(fnamelist,args{:});
    else
      results = mmil_interscan_stdev(fnamelist,args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
