function fname_out = fs_nu_corr(fname_in,varargin)
%function fname_out = fs_nu_corr(fname_in,varargin)
%
% Purpose: perform N3 non-uniformity intensity correction on T1-weighted image
%
% Usage: fs_nu_corr(fname_in,'key1',value1,...)
%
% Required Parameters:
%  fname_in: input T1-weighted file name
%
% Optional Parameters:
%  'fname_out': output file name
%    {default = 'nu.mgz'}
%  'nu3T_flag': [0|1] whether to use nu intensity correction
%    with settings optimized for 3T
%    {default = 0}
%  'tal_flag': [0|1] register to Talairach space
%    for standardization of white matter intensity values
%    {default = 1}
%  'fname_talxfm': [0|1] existing Talairach transform file
%    if supplied, will skip Talairach registration
%    {default = []}
%  'tmpdir': temporary directory containing intermediate output
%    {default = 'tmp_nu_corr'}
%  'cleanup_flag': [0|1] remove temporary files before quitting
%    {default = 1}
%  'verbose': [0|1] display status messages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  12/29/14 by Don Hagler
% Last Mod: 11/27/15 by Don Hagler
%

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

% N3 nonuniformity correction
% registration to Talairach
if parms.tal_flag && isempty(parms.fname_talxfm)
  parms.fname_talxfm = tal_reg(parms);
else
  parms.fname_talxfm = [];
end;
parms.fname_in = nu_correction(parms,'nu');

% copy to final destination
fname_out = copy_output(parms);

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
    'fname_out','nu.mgz',[],...
    'nu3T_flag',false,[false true],...
    'tal_flag',true,[false true],...
    'fname_talxfm',[],[],...
    'tmpdir','tmp_nu_corr',[],...
    'cleanup_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  if parms.tal_flag &&...
     ~isempty(parms.fname_talxfm) && ~exist(parms.fname_talxfm,'file')
    error('file %s not found',parms.fname_talxfm);
  end;
  if mmil_isrelative(parms.tmpdir)
    parms.tmpdir = [pwd '/' parms.tmpdir];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = tal_reg(parms)
  % initial nu correction
  fname_tal_input = nu_correction(parms,'orig_nu');
  % talairach registration  
  fname_out = sprintf('%s/talairach.xfm',parms.tmpdir);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: registering to MNI Talairach...\n',mfilename);
      tic;
    end;
    cmd = sprintf('cd %s',parms.tmpdir);
    cmd = sprintf('%s; talairach_avi',cmd);
    cmd = sprintf('%s --i %s --xfm %s',...
      cmd,fname_tal_input,fname_out);
    [status,msg] = unix(cmd);
    if status || ~exist(fname_out,'file')
      error('MNI Talairach registration failed:\n%s\n%s\n',cmd,msg);
    end;
    if parms.verbose
      toc;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = nu_correction(parms,outstem)
  fname_out = sprintf('%s/%s.mgz',parms.tmpdir,outstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: correcting intensity non-uniformity...\n',mfilename);
      tic;
    end;
    cmd = 'mri_nu_correct.mni';
    cmd = sprintf('%s --i %s --o %s',cmd,parms.fname_in,fname_out);
    if strcmp(outstem,'orig_nu')
      cmd = sprintf('%s --proto-iters 1000 --distance 50 --n 1 --no-rescale',cmd);
    elseif parms.nu3T_flag  
      cmd = sprintf('%s --proto-iters 1000 --distance 50 --n 1',cmd);
    else
      cmd = sprintf('%s --n 2',cmd);
    end;
    if ~isempty(parms.fname_talxfm)
      cmd = sprintf('%s --uchar %s',cmd,parms.fname_talxfm);
    end;
    [status,msg] = unix(cmd);
    if status || ~exist(fname_out,'file')
      error('non-uniformity correction failed:\n%s\n%s\n',cmd,msg);
    end;
    if parms.verbose
      toc;
    end;
    % add transform file to header
    if ~isempty(parms.fname_talxfm)
      if parms.verbose
        fprintf('%s: adding transform to header...\n',mfilename);
        tic;
      end;
      cmd = sprintf('mri_add_xform_to_header -c %s %s %s\n',...
        parms.fname_talxfm,fname_out,fname_out);
      [status,msg] = unix(cmd);
      if status || ~exist(fname_out,'file')
        error('failed to add transform to header:\n%s\n%s\n',cmd,msg);
      end;
      if parms.verbose
        toc;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = copy_output(parms)
  fname_out = parms.fname_out;
  if ~strcmp(fname_out,parms.fname_in)
    cmd = sprintf('cp %s %s',parms.fname_in,fname_out);
    [status,msg] = unix(cmd);
    if status
      error('failed to copy output file:\n%s\n%s\n',cmd,msg);
    end;
  end;
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

