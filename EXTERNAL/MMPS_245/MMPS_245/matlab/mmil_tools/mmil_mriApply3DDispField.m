function mmil_mriApply3DDistField(fname_in,fname_out,fname_dx,fname_dy,fname_dz,varargin)
%function mmil_mriApply3DDistField(fname_in,fname_out,...
%                                  fname_dx,fname_dy,fname_dz,[options])
%
% Purpose: apply 3D displacement fields to a volume
%   This is a wrapper for Dominic Holland's mriApply3DDistField
%
% Required Input:
%   fname_in:  file name of input volume
%   fname_out: file name of output volume
%   fname_dx:  file name of x displacement volume
%   fname_dy:  file name of y displacement volume
%   fname_dz:  file name of z displacement volume
%
% Optional Input:
%   'binfile' - full or relative path of binary file
%     {default = 'mriApply3DDispField'}
%   'sincflag' - [0|1] whether to use sinc interpolation (otherwise linear)
%     {default = 0}
%   'forceflag' - [0|1] overwrite existing output files
%     {default = 0}
% 
% Created:  09/08/09 by Don Hagler
% Last Mod: 04/02/10 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,5)) return; end;
parms = mmil_args2parms(varargin, { ...
  'bindir',[],[],...
  'binfile','mriApply3DDispField',[],...
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

% check input files exist
if ~exist(fname_in,'file'), error('file %s not found',fname_in); end;
if ~exist(fname_dx,'file'), error('file %s not found',fname_dx); end;
if ~exist(fname_dy,'file'), error('file %s not found',fname_dy); end;
if ~exist(fname_dz,'file'), error('file %s not found',fname_dz); end;

if ~exist(fname_out,'file') || parms.forceflag
  % run mriApply3DDispField
  cmd = sprintf('%s -i %s -o %s -dx %s -dy %s -dz %s %s',...
    parms.binfile,...
    fname_in,fname_out,...
    fname_dx,fname_dy,fname_dz);
  if parms.sincflag
    cmd = sprintf('%s -sinc',cmd);
  end;
  % execute command
  [status,result] = unix(cmd);
  if status
    error('cmd %s failed:\n%s',cmd,result);
  else
    fprintf('%s: cmd %s succeeded:\n%s',mfilename,cmd,result);
  end;
end;

