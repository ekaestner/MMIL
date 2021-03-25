function fname_out = mmil_headmask(fname_in,varargin)
%function fname_out = mmil_headmask(fname_in,varargin)
%
% Purpose: create head mask volume from T1-weighted image
%   using scalp surface from FreeSurfer's mri_watershed
%
% Usage: mmil_headmask(fname_in,'key1',value1,...)
%
% Required Parameters:
%  fname_in: input T1-weighted file name
%
% Optional Parameters:
%  'fname_out': output file name
%    {default = 'headmask.mgz'}
%  'tmpdir': temporary directory containing intermediate output
%    {default = 'tmp_headmask'}
%  'cleanup_flag': [0|1] remove temporary files before quitting
%    {default = 1}
%  'verbose': [0|1] display status messages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  03/22/16 by Don Hagler
% Last Mod: 03/22/16 by Don Hagler
%

%% todo: option to also save brainmask?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname_in,varargin);

% check output file
if exist(parms.fname_out) && ~parms.forceflag
  fname_out = parms.fname_out;
  return;
end;

% create temporary output directory
mmil_mkdir(parms.tmpdir);

% create brain mask
parms = create_brainmask(parms);

% create head mask
parms = create_headmask(parms);

% resample to original slicing and resolution
parms = resamp_to_raw(parms);

% copy to final output
fs_copy_mgh(parms.fname_headmask_raw,parms.fname_out);
fname_out = parms.fname_out;

% remove temporary files
if parms.cleanup_flag
  parms = cleanup_tmpdir(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_out','headmask.mgz',[],...
    'tmpdir','tmp_headmask',[],...
    'cleanup_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ...
    'brainmask_tags',{'fname_out','conform_flag','nu_flag','T1_flag',...
                      'surfs_flag','nu3T_flag','tal_flag','res2raw_flag',...
                      'conform_options','fname_talxfm','surfs_outdir',...
                      'surfs_nverts','watershed_options','tmpdir',...
                      'cleanup_flag','verbose','forceflag'},[],...
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  [fpath,parms.fstem] = fileparts(parms.fname_in);
  if mmil_isrelative(parms.tmpdir)
    parms.tmpdir = [pwd '/' parms.tmpdir];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_brainmask(parms)
  tparms = parms;
  tparms.fname_out = sprintf('%s/%s_brainmask.mgz',parms.tmpdir,parms.fstem);
  tparms.surfs_outdir = [parms.tmpdir '/brainmask_surfs'];
  tparms.res2raw_flag = false;
  tparms.cleanup_flag = false;
  args = mmil_parms2args(tparms,parms.brainmask_tags);
  fname_out = fs_brainmask(parms.fname_in,args{:});
  parms.fname_brainmask = fname_out;
  parms.fname_scalpsurf = sprintf('%s/outer_scalp.tri',tparms.surfs_outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_headmask(parms)
  fname_out = sprintf('%s/%s_headmask.mgz',parms.tmpdir,parms.fstem);
  if ~exist(fname_out,'file') || forceflag
    if parms.verbose
      fprintf('%s: creating head mask...\n',mfilename);
    end;
    vol = ctx_load_mgh(parms.fname_brainmask);
    surf = fs_read_trisurf(parms.fname_scalpsurf);
    surf.vertices(:,1) = -surf.vertices(:,1) + vol.lphcent(1);
    surf.vertices(:,2) = -surf.vertices(:,2) + vol.lphcent(2);
    surf.vertices(:,3) = surf.vertices(:,3) + vol.lphcent(3);
    msurf = preprocessQ(surf);
    [tmp,vol_mask] = getmaskvol(vol,msurf,eye(4));
    ctx_save_mgh(vol_mask,fname_out);
  end;
  parms.fname_headmask = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resamp_to_raw(parms)
  fname_out = sprintf('%s/%s_headmask_raw.mgz',parms.tmpdir,parms.fstem);
  fs_mri_convert(parms.fname_headmask,fname_out,...
    'options',sprintf('-rl %s',parms.fname_in),...
    'forceflag',parms.forceflag);
  parms.fname_headmask_raw = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = cleanup_tmpdir(parms)
  if exist(parms.tmpdir,'dir')
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,msg] = unix(cmd);
    if status
      error('tmpdir cleanup failed:\n%s\n%s\n',cmd,msg);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

