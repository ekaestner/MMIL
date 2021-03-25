function mmil_mriApplyTransforms(fname_in,varargin)
%function mmil_mriApplyTransforms(fname_in,[options])
%
% Purpose: apply affine transformation matrix to a volume
%   This is a wrapper for Dominic Holland's mriApplyTransforms
%
% Required Input:
%   fname_in:  file name of input volume
%
% Optional Input:
%   'fname_out': file name of output volume
%     {default = []}
%   'fname_dx':  file name of x displacement volume
%     {default = []}
%   'fname_dy':  file name of y displacement volume
%     {default = []}
%   'fname_dz':  file name of z displacement volume
%     {default = []}
%   'fname_mat':  file name of affine registration matrix text file
%     {default = []}
%   'invflag': [0|1] whether to invert the affine reg matrix
%     {default = 0}
%   'fname_targ':  file name of target volume
%     {default = []}
%   'binfile' - full or relative path of binary file
%     {default = 'mriApplyTransforms'}
%   'sincflag' - [0|1] whether to use sinc interpolation (otherwise linear)
%     {default = 1}
%   'forceflag' - [0|1] overwrite existing output files
%     {default = 0}
% 
% Created:  04/22/10 by Don Hagler
% Last Mod: 04/22/10 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,5)) return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_in',fname_in,[],...
  'fname_out',[],[],...
  'fname_dx',[],[],...
  'fname_dy',[],[],...
  'fname_dz',[],[],...
  'fname_mat',[],[],...
  'invflag',false,[false true],...
  'fname_targ',[],[],...
  'bindir',[],[],...
  'binfile','mriApplyTransforms',[],...
  'sincflag',false,[false true],...
  'forceflag',false,[false true],...
});

if isempty(parms.bindir)
  parms.bindir = [getenv('MMPS_DIR') '/bin'];
end;
if isempty(regexp(parms.binfile,'^/')) % relative path
  parms.binfile = [parms.bindir '/' parms.binfile];
end;

% check binary exists
if ~exist(parms.binfile,'file')
  error('binary file %s not found',parms.binfile);
end;

% check if output already exists
if ~isempty(parms.fname_out) & exist(parms.fname_out,'file') & ~parms.forceflag
  return;
end;

% check input files exist
if ~exist(parms.fname_in,'file'), error('file %s not found',parms.fname_in); end;

if ~isempty(parms.fname_dx) & ~exist(parms.fname_dx,'file')
  error('file %s not found',parms.fname_dx);
end;
if ~isempty(parms.fname_dy) & ~exist(parms.fname_dy,'file')
  error('file %s not found',parms.fname_dy);
end;
if ~isempty(parms.fname_dz) & ~exist(parms.fname_dz,'file')
  error('file %s not found',parms.fname_dz);
end;
if ~isempty(parms.fname_mat) & ~exist(parms.fname_mat,'file')
  error('file %s not found',parms.fname_mat);
end;
if ~isempty(parms.fname_targ) & ~exist(parms.fname_targ,'file')
  error('file %s not found',parms.fname_targ);
end;
  
if ~exist(parms.fname_out,'file') || parms.forceflag
  % run mriApplyTransforms
  cmd = sprintf('%s -f %s',parms.binfile,parms.fname_in);
  if ~isempty(parms.fname_dx)
    cmd = sprintf('%s -dx %s',cmd,parms.fname_dx);
  end;
  if ~isempty(parms.fname_dy)
    cmd = sprintf('%s -dy %s',cmd,parms.fname_dy);
  end;
  if ~isempty(parms.fname_dz)
    cmd = sprintf('%s -dz %s',cmd,parms.fname_dz);
  end;
  if ~isempty(parms.fname_mat)
    cmd = sprintf('%s -m %s',cmd,parms.fname_mat);
    if parms.invflag
      cmd = sprintf('%s -inv',cmd);
    end;
  end;
  if ~isempty(parms.fname_targ)
    cmd = sprintf('%s -t %s',cmd,parms.fname_targ);
  end;
  if parms.sincflag
    cmd = sprintf('%s -sinc',cmd);
  end;
  if ~isempty(parms.fname_out)
    cmd = sprintf('%s -o %s',cmd,parms.fname_out);
  end;
  % execute command
  [status,result] = unix(cmd);
  if status
    error('cmd %s failed:\n%s',cmd,result);
  else
    fprintf('%s: cmd %s succeeded:\n%s',mfilename,cmd,result);
  end;
end;

