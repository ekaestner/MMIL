function [fname_out,fname_motion_mat,fname_motion_1D] = bold_motion_corr(fname,varargin)
%function [fname_out,fname_motion_mat,fname_motion_1D] = bold_motion_corr(fname,[options])
%
% Purpose: Correct head motion for BOLD scans using
%   rigid-body volume registration between data and reference
%
% Usage:
%  [fname_out,fname_motion_mat,fname_motion_1D] = bold_motion_corr(fname,'key1', value1,...);
%
% Required Input
%  fname: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Parameters:
%  'fname_out': output file name
%    if empty, will add '_mc' to file stem of input file name
%    {default = []}
%  'fname_motion_mat': output file name of mat file containing registration
%    matrices M_ref_to_orig for each frame
%    if empty, will append '_motion.mat' to file stem of fname
%    {default = []}
%  'fname_motion_1D': output file name of text file containing motion parameters
%    if empty, will append '_motion.1D' to file stem of fname
%    {default = []}
%  'afni_flag': use AFNI's 3dvolreg to estimate and correct for motion
%    otherwise use bold_estimate_motion (slower, minimum step size)
%      and mmil_resample_vol
%    {default = 1}
%  'interpm': interpolation method
%     0 = nearest neighbor, 1 = linear, 2 = cubic
%     3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%     if resolution or slice position differs for fname_in and fname_ref
%       and afni_flag=1, will resample volume before motion correction
%    { default: 2 }
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': overwrite existing output files
%     {default: 0}
%
% Optional Parameters for motion estimation (only used if afni_flag = 0):
%  'intra_flag': perform within-scan motion correction
%     {default = 1}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    if empty, will use fname as reference
%    {default = []}
%  'fname_mask': name of mask volume
%    if empty, will generate one from fname_ref using mmil_quick_brainmask
%    {default = []}
%  'ref_frame': reference frame number (1-based)
%    {default = 1}
%  'min_trans': minimum translation (mm)
%     used by bold_estimate_motion, not 3dvolreg
%     {default = 0.05}
%  'min_rot': minimum rotation (degrees)
%     used by bold_estimate_motion, not 3dvolreg
%     {default = 0.05}
%
% Output:
%   fname_out: output file name of motion corrected 4D volume
%   fname_motion_mat: output file name of mat file with registration matrices
%   fname_motion_1D: output file name of text file containing motion parameters
%
% Created:  04/26/10 by Don Hagler
% Last Mod: 10/31/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_motion_mat = []; fname_motion_1D = [];
if ~mmil_check_nargs(nargin,2), return; end;

parms = check_input(fname,varargin);

if output_exists(parms)
  fname_out = parms.fname_out;
  fname_motion_mat = parms.fname_motion_mat;
  fname_motion_1D = parms.fname_motion_1D;
  return;
end;

if parms.afni_flag
  parms = check_mismatch(parms);
  [fname_out,fname_motion_mat,fname_motion_1D] = correct_motion(parms);
else
  [parms,fname_motion_mat,fname_motion_1D] = estimate_motion(parms);
  fname_out = correct_estimated_motion(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
  ... % optional
    'fname_out',[],[],...
    'fname_motion_mat',[],[],...
    'fname_motion_1D',[],[],...
    'afni_flag',true,[false true],...
    'interpm',2,[0:5],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'intra_flag',true,[false true],...
    'fname_ref',fname,[],...
    'ref_frame',1,[1 Inf],...
    'fname_mask',[],[],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
  ... % hidden
    'tmpdir','tmp_mc',[],...
    'tmpext','.mgh',[],...
    'cleanupflag',true,[false true],...
    'mstep',4,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
    'motion_radius',50,[],... % for calculating distance from angle
  ...
    'estimate_tags',{'intra_flag','fname_ref','ref_frame','fname_mask',...
                     'min_trans','min_rot','verbose','mstep','scales'},[],...
    'resamp_tags',{'fname_out','save_reg_flag','fname_reg','qmat',...
      'fname_qmat','native_flag','nvoxels','resolution','deoblique_flag',...
      'std_orient','M_reg','M_ref','volsz_ref','regT1flag','fname_T1',...
      'M_T1_to_EPI','smooth','rot','trans','interpm','bclamp','EPI_type',...
      'ext','verbose','forceflag'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~exist(fname,'file')
    error('file %s not found',fname);
  end;
  if ~exist(parms.fname_ref,'file')
    error('file %s not found',parms.fname_ref);
  end;

  [fpath,fstem,fext] = fileparts(parms.fname_in);
  if ~isempty(parms.fname_out)
    parms.outdir = fileparts(parms.fname_out);
  else
    parms.outdir = fpath;
  end;
  if mmil_isrelative(parms.tmpdir)
    parms.tmpdir = [parms.outdir '/' parms.tmpdir];
  end;
  if isempty(parms.fname_out)
    parms.fname_out = [parms.outdir '/' fstem '_mc' fext];
  end;
  if isempty(parms.fname_motion_mat)
    parms.fname_motion_mat = [parms.outdir '/' fstem '_motion.mat'];
  end;
  if isempty(parms.fname_motion_1D)
    parms.fname_motion_1D = [parms.outdir '/' fstem '_motion.1D'];
  end;
  if parms.afni_flag
    parms.fname_res = [parms.tmpdir '/' fstem '_res' parms.tmpext];
  end;

  % get info about input volume
  [parms.M,parms.volsz] = ...
    mmil_load_mgh_info(parms.fname_in,parms.forceflag);
  parms.nframes = parms.volsz(4);

  % get info about reference volume
  [parms.M_ref,parms.volsz_ref] = ...
    mmil_load_mgh_info(parms.fname_ref,parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exist_flag = output_exists(parms)
  if parms.forceflag
    exist_flag = 0;
    return;
  else
    exist_flag = 1;
  end;
  if ~exist(parms.fname_out,'file')
    exist_flag = 0;
  end;
  if ~exist(parms.fname_motion_mat,'file')
    exist_flag = 0;
  end;
  if ~exist(parms.fname_motion_1D,'file')
    exist_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_mismatch(parms)
  % compare M and M_ref, volsz and volsz_ref
  if all(parms.M(:)==parms.M_ref(:)) &&...
     any(parms.volsz(1:3)==parms.volsz_ref(1:3))
    return;
  end;
  % resample input volume to match ref
  if ~exist(parms.fname_res,'file') || parms.forceflag
    mmil_mkdir(parms.tmpdir);
    tparms = parms;
    tparms.fname_out = parms.fname_res;
    tparms.nvoxels = parms.volsz_ref(1:3);
    tparms.resolution = fs_voxsize(parms.M_ref);
    tparms.deoblique_flag = 0;
    %% todo: test mcflirt
    %%       maybe it handles between-scan motion better than 3dvolreg
    %% todo: use bold_estimate_motion with intra_flag = 0
    %%       needed if two scans from different sessions
    %%       or possibly if after a new localizer or prescription
    args = mmil_parms2args(tparms,parms.resamp_tags);
    epi_resample_data(parms.fname_in,args{:});
  end;
  parms.fname_in = parms.fname_res;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,fname_motion_mat,fname_motion_1D] = correct_motion(parms)
  fname_out = [];
  fname_motion_mat = [];
  fname_motion_1D = [];
  if ~exist(parms.fname_out,'file') ||...
     ~exist(parms.fname_motion_1D,'file') || parms.forceflag
    mmil_mkdir(parms.tmpdir);
    if parms.volsz(4)==1 && strcmp(parms.fname_in,parms.fname_ref)
      if parms.verbose
        fprintf('%s: copying %s to %s\n',...
          mfilename,parms.fname_in,parms.fname_out);
      end;
      % copy from mgz to mgh
      fs_copy_mgh(parms.fname_in,parms.fname_out);
      % create motion.1D file with 0s
      fid = fopen(parms.fname_motion_1D,'wt');
      if fid<0
        error('failed to open %s for writing',parms.fname_motion_1D);
      end;
      fprintf(fid,'0 0 0 0 0 0 0 1 1\n');
      fclose(fid);
    else
      if parms.verbose
        fprintf('%s: using frame %d from %s as motion reference for %s...\n',...
          mfilename,parms.ref_frame,parms.fname_ref,parms.fname_in);
      end;
      % convert reference scan to nii
      [fpath_ref,fstem_ref,fext_ref] = fileparts(parms.fname_ref);
      fname_ref_tmp = [parms.tmpdir '/' fstem_ref '_ref.nii'];
      conv_options = sprintf('-f %d',parms.ref_frame-1);
      fs_mri_convert(parms.fname_ref,fname_ref_tmp,'options',conv_options,...
        'forceflag',parms.forceflag);
      % convert input image to nii
      [fpath,fstem,fext] = fileparts(parms.fname_in);
      fname_in_tmp = [parms.tmpdir '/' fstem '.nii'];
      fs_mri_convert(parms.fname_in,fname_in_tmp,'forceflag',parms.forceflag);
      % correct for head motion using 3dvolreg
      fname_out_tmp = [parms.tmpdir '/' fstem '_mc.nii'];
      if ~exist(fname_out_tmp,'file') | parms.forceflag
        if parms.verbose
          fprintf('%s: correcting motion for %s\n',mfilename,parms.fname_in);
        end;
        cmd = sprintf('3dvolreg -prefix %s',fname_out_tmp);
        cmd = sprintf('%s -base ''%s''',cmd,fname_ref_tmp);
        cmd = sprintf('%s -dfile %s',cmd,parms.fname_motion_1D);
        cmd = sprintf('%s %s',cmd,fname_in_tmp);
        if parms.verbose
          fprintf('%s: cmd = %s\n',mfilename,cmd);
        end;
        [status,result] = unix(cmd);
        if status
          error('3dvolreg failed:\n%s',result);
        end;
      end;
      % convert nifti to mgh
      fs_mri_convert(fname_out_tmp,parms.fname_out);
    end;
  end;

  % create M_motion from fname_motion_1D
  if ~exist(parms.fname_motion_mat,'file') || parms.forceflag
    motion_tseries_afni = mmil_load_motion_1D(parms.fname_motion_1D);
    nframes = size(motion_tseries_afni,1);
    M_motion = cell(nframes,1);
    motion_tseries = zeros(nframes,6);
    for f=1:size(motion_tseries_afni,1)
      tmp_motion = motion_tseries_afni(f,:);
      trans = tmp_motion(1:3);
      rot   = tmp_motion(4:6);
      M_reg = mmil_construct_M('trans',trans,'rot',rot,...
        'order',{'rot','trans'},'rot_order',{'z','x','y'});
      M_motion{f} = inv(M_reg);
      % recalculate motion_tseries
      tmp_motion = mmil_M_mat2vec(M_reg);
      motion_tseries(f,:) = tmp_motion;
    end;
    % calculate absolute difference in head position or rotation
    [mean_motion,mean_trans,mean_rot] = ...
      mmil_mean_motion(motion_tseries,parms.motion_radius);
    save(parms.fname_motion_mat,'M_motion',...
      'motion_tseries_afni','motion_tseries',...
      'mean_motion','mean_trans','mean_rot');
  end;
  
  % delete temporary files
  if parms.cleanupflag & exist(parms.tmpdir,'dir')
    cmd = sprintf('rm -r %s\n',parms.tmpdir);
    [status,result] = unix(cmd);
    if status
      warning('failed to remove tmp dir %s:\n%s',parms.tmpdir,result);
    end;
  end;
  
  fname_out = parms.fname_out;
  fname_motion_1D = parms.fname_motion_1D;
  fname_motion_mat = parms.fname_motion_mat;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,fname_motion_mat,fname_motion_1D] = estimate_motion(parms)
  fname_motion_mat = [];
  fname_motion_1D = [];
  
  % estimate motion for each frame
  if ~exist(parms.fname_motion_mat,'file') || parms.forceflag
    if ~isempty(parms.fname_mask)
      if ~exist(parms.fname_mask) || parms.forceflag
        mmil_quick_brainmask([],'fname_in',parms.fname_ref,...
          'fname_out',parms.fname_mask,...
          'forceflag',parms.forceflag);
      end;
    end;
    args = mmil_parms2args(parms,parms.estimate_tags);
    [M_motion,motion_tseries] = ...
      bold_estimate_motion(parms.fname_in,args{:});
    % calculate absolute difference in head position or rotation
    [mean_motion,mean_trans,mean_rot] = ...
      mmil_mean_motion(motion_tseries,parms.motion_radius);
    save(parms.fname_motion_mat,'M_motion','motion_tseries',...
      'mean_motion','mean_trans','mean_rot');
  else
    load(parms.fname_motion_mat);
  end;
  parms.M_motion = M_motion;
  fname_motion_mat = parms.fname_motion_mat;

  % save motion_tseries to motion_1D file
  if ~exist(parms.fname_motion_1D,'file') || parms.forceflag
    fid = fopen(parms.fname_motion_1D,'wt');
    if fid==-1
      error('failed to open %s for writing',parms.fname_motion_1D);
    end;
    for f=1:size(motion_tseries,1)
      fprintf(fid,'%4d%s\n',f-1,sprintf('%10.4f',motion_tseries(f,:)));
    end;
    fclose(fid);
  end;
  fname_motion_1D = parms.fname_motion_1D;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = correct_estimated_motion(parms)
  if ~exist(parms.fname_out,'file') || parms.forceflag
    [vol_orig,M] = fs_load_mgh(parms.fname_in);
    % initialize output
    vol = zeros([parms.volsz_ref(1:3),parms.nframes],'single');
    if parms.verbose
      if parms.nframes>1
        plurstr = 's';
      else
        plurstr = [];
      end;      
      fprintf('%s: resampling %d frame%s...\n',...
        mfilename,parms.nframes,plurstr);
    end;
    for f = 1:parms.nframes
      vol_tmp = squeeze(vol_orig(:,:,:,f));
      vol_tmp = mmil_resample_vol(vol_tmp,M,...
        'M_ref',parms.M_ref,'nvox_ref',parms.volsz_ref,...
        'M_reg',parms.M_motion{f},...
        'interpm',parms.interpm);
      vol(:,:,:,f) = vol_tmp;
    end;
    fs_save_mgh(vol,parms.fname_out,parms.M_ref);
  end;
  fname_out = parms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

