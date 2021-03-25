function fname_out = mmil_reg2mni(fname_head,fname_brain,varargin)
%function fname_out = mmil_reg2mni(fname_head,fname_brain,[options])
%
% Required Input:
%   fname_head: full path of head image
%   fname_brain:  full path of brain image
%
% Optional Parameters:
%   'fname_in': full path of input image to be warped to atlas
%     if empty, will apply warp to fname_head
%     must be in register with fname_head
%     {default = []}
%   'fname_out': full path of output image warped to atlas
%     if empty, will create in outdir with name based on fname_in
%     {default = []}
%   'outdir': output directory
%     {default = [pwd '/reg2mni']}
%   'verbose': [0|1] display status updates
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_out: output file nonlinear warped to atlas
%
% Created:  11/16/17 by Don Hagler
% Last Mod: 11/16/17 by Don Hagler
%

%% todo: use full paths instead of cd'ing into outdir and using local fstems

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];
if ~mmil_check_nargs(nargin,2), return; end;

% check input parameters
parms = check_input(fname_head,fname_brain,varargin);

% check output file
if exist(parms.fname_out,'file') && ~parms.forceflag
  fname_out = parms.fname_out;
  return;
end;

% copy input images
parms = copy_input(parms);

% copy atlas images
parms = copy_atlas(parms);

% affine registration of brain to atlas
parms = affine_reg(parms);

% nonlinear registration to head to atlas
parms = nonlin_reg(parms);

% apply nonlinear registration
parms = apply_warp(parms);

fname_out = parms.fname_out;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_head,fname_brain,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'fname_head',fname_head,[],...
    'fname_brain',fname_brain,[],...
  ...
    'fname_in',[],[],...
    'fname_out',[],[],...
    'outdir',[pwd '/reg2mni'],[],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'next','.nii',[],...
    'orient','RAS',[],...
    'atlas','T1_2_MNI152_2mm',[],...
    'dof',12,[],...
    'warpres',10,[],...
    'interp','trilinear',[],...
   ...
    'fstem_head','head',[],...
    'fstem_brain','brain',[],...
    'fstem_ref','standard',[],...
  });
  if ~exist(parms.fname_head,'file')
    error('file %s not found',parms.fname_head);
  end;
  if ~exist(parms.fname_brain,'file')
    error('file %s not found',parms.fname_brain);
  end;
  if isempty(parms.fname_in)
    parms.fname_in = parms.fname_head;
  end;
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  [fpath,parms.fstem_in] = fileparts(parms.fname_in);
  if isempty(parms.fname_out)
    parms.fname_out = sprintf('%s/%s_reg2mni%s',parms.outdir,parms.fstem_in,parms.next);
  end;
  parms.fstem_head_atlas = sprintf('%s/data/standard/MNI152_T1_2mm',getenv('FSLDIR'));
  parms.fstem_brain_atlas = sprintf('%s_brain',parms.fstem_head_atlas);
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = copy_input(parms)
  parms.fname_head = convert_mgh2nifti(parms.fname_head,...
    sprintf('%s/%s',parms.outdir,parms.fstem_head),parms);
  parms.fname_brain = convert_mgh2nifti(parms.fname_brain,...
    sprintf('%s/%s',parms.outdir,parms.fstem_brain),parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = convert_mgh2nifti(fname,outstem,parms)
  [tpath,tstem] = fileparts(outstem);
  fname_out = sprintf('%s/%s%s',parms.outdir,tstem,parms.next);
  if ~exist(fname_out,'file') || parms.forceflag
    fs_mri_convert(fname,fname_out,...
      'out_orient',parms.orient,'forceflag',parms.forceflag);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = copy_atlas(parms)
  fstem_brain_atlas = sprintf('%s_%s',parms.fstem_ref,parms.fstem_brain);
  fname_brain_atlas = sprintf('%s/%s%s',parms.outdir,fstem_brain_atlas,parms.next);
  if ~exist(fname_brain_atlas,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: copying brain atlas...\n',mfilename);
    end;
    cmd = sprintf('cd %s',parms.outdir);
    cmd = sprintf('%s; fslmaths %s %s',...
      cmd,parms.fstem_brain_atlas,fstem_brain_atlas);
    run_cmd(cmd,parms);
  end;
  
  fstem_head_atlas = sprintf('%s_%s',parms.fstem_ref,parms.fstem_head);
  fname_head_atlas = sprintf('%s/%s%s',parms.outdir,fstem_head_atlas,parms.next);
  if ~exist(fname_head_atlas,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: copying head atlas...\n',mfilename);
    end;
    cmd = sprintf('cd %s',parms.outdir);
    cmd = sprintf('%s; fslmaths %s %s',...
      cmd,parms.fstem_head_atlas,fstem_head_atlas);
    run_cmd(cmd,parms);
  end;

  fstem_mask_atlas = sprintf('%s_mask',parms.fstem_ref);
  fname_mask_atlas = sprintf('%s/%s%s',parms.outdir,fstem_mask_atlas,parms.next);
  if ~exist(fname_mask_atlas,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: creating atlas mask...\n',mfilename);
    end;
    cmd = sprintf('cd %s',parms.outdir);
    cmd = sprintf('%s; fslmaths %s -bin -dilF -dilF %s -odt char',...
      cmd,parms.fstem_brain_atlas,fstem_mask_atlas);
    run_cmd(cmd,parms);
  end;

  parms.fstem_brain_atlas = fstem_brain_atlas;
  parms.fstem_head_atlas = fstem_head_atlas;
  parms.fstem_mask_atlas = fstem_mask_atlas;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_cmd(cmd,parms)
  if parms.verbose
    fprintf('%s\n',cmd);
  end;
  [s,r] = unix(cmd);
  if s
    error('cmd %s failed:\n%s',cmd,r);
  elseif parms.verbose
    fprintf('%s\n',r);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = affine_reg(parms)
  parms.fstem_brain2ref = sprintf('%s2%s_linear',...
    parms.fstem_brain,parms.fstem_ref);
  parms.fstem_ref2brain = sprintf('%s2%s_linear',...
    parms.fstem_ref,parms.fstem_brain);
  parms.fname_affine_brain2ref = sprintf('%s/%s.mat',parms.outdir,parms.fstem_brain2ref);
  parms.fname_affine_ref2brain = sprintf('%s/%s.mat',parms.outdir,parms.fstem_ref2brain);
  if ~exist(parms.fname_affine_brain2ref,'file') || ...
      ~exist(parms.fname_affine_ref2brain,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: affine registration of brain to atlas...\n',mfilename);
    end;
    cmd = sprintf('cd %s',parms.outdir);
    cmd = sprintf('%s; flirt -ref %s -in %s -out %s',...
      cmd,parms.fstem_brain_atlas,parms.fstem_brain,parms.fstem_brain2ref);
    cmd = sprintf('%s -omat %s.mat -interp %s -cost corratio -dof %d',...
      cmd,parms.fstem_brain2ref,parms.interp,parms.dof);
    cmd = sprintf('%s; convert_xfm -inverse -omat %s %s',...
      cmd,parms.fname_affine_ref2brain,parms.fname_affine_brain2ref);
    run_cmd(cmd,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = nonlin_reg(parms)
  parms.fstem_head2ref = sprintf('%s2%s_nonlinear',parms.fstem_head,parms.fstem_ref);
  parms.fname_nonlin_head2ref = sprintf('%s/%s_warp%s',...
    parms.outdir,parms.fstem_head2ref,parms.next);
  if ~exist(parms.fname_nonlin_head2ref,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: nonlinear registration of head to atlas...\n',mfilename);
    end;
    cmd = sprintf('cd %s',parms.outdir);
    cmd = sprintf('%s; fnirt --in=%s --aff=%s.mat --cout=%s_warp',...
      cmd,parms.fstem_head,parms.fstem_brain2ref,parms.fstem_head2ref);
    cmd = sprintf('%s --iout=%s --jout=%s_jac --config=%s',...
      cmd,parms.fstem_head2ref,parms.fstem_head2ref,parms.atlas);
    cmd = sprintf('%s --ref=%s --refmask=%s --warpres=%d,%d,%d',...
      cmd,parms.fstem_head_atlas,parms.fstem_mask_atlas,...
      parms.warpres,parms.warpres,parms.warpres);
    run_cmd(cmd,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = apply_warp(parms)
  if ~exist(parms.fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: applying warp to mni...\n',mfilename);
    end;
    fname_tmp = convert_mgh2nifti(parms.fname_in,...
      mmil_tempfname(parms.fstem_in,parms.outdir),parms);
    [tmp,fstem_tmp] = fileparts(fname_tmp);
    cmd = sprintf('cd %s',parms.outdir);
    cmd = sprintf('%s; applywarp --ref=%s --in=%s --out=%s_reg2mni --warp=%s_warp',...
      cmd,parms.fstem_head_atlas,fstem_tmp,parms.fstem_in,parms.fstem_head2ref);
    fstem_in = sprintf('%s_reg2mni',parms.fstem_in);
    [tpath,tstem] = fileparts(parms.fname_out);
    if ~strcmp(tpath,parms.outdir) || ~strcmp(tstem,fstem_in)
      cmd = sprintf('%s; immv %s %s',cmd,fstem_in,parms.fname_out);
    end;
    run_cmd(cmd,parms);
    delete(fname_tmp);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

