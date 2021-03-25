function fname_out = fs_T1_corr(fname_in,varargin)
%function fname_out = fs_T1_corr(fname_in,varargin)
%
% Purpose: perform T1 intensity normalization using mri_normalize
%
% Usage: fs_T1_corr(fname_in,'key1',value1,...)
%
% Required Parameters:
%  fname_in: input T1-weighted file name
%
% Optional Parameters:
%  'fname_out': output file name
%    {default = 'wmseg.mgz'}
%  'fname_mask': input brain mask file name
%    {default = []}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  12/29/14 by Don Hagler
% Last Mod: 02/19/15 by Don Hagler
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

% T1 intensity normalization
fname_out = T1_correction(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_out','T1.mgz',[],...
    'fname_mask',[],[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  if ~isempty(parms.fname_mask) && ~exist(parms.fname_mask,'file')
    error('file %s not found',parms.fname_mask);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = T1_correction(parms)
  fname_out = parms.fname_out;
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: normalizing T1 intensities...\n',mfilename);
      tic;
    end;
    cmd = 'mri_normalize';
    if ~isempty(parms.fname_mask)
      cmd = sprintf('%s -mask %s',cmd,parms.fname_mask);
    end;
    cmd = sprintf('%s %s %s',cmd,parms.fname_in,fname_out);
    [status,msg] = unix(cmd);
    if status || ~exist(fname_out,'file')
      error('T1 normalization failed:\n%s\n%s\n',cmd,msg);
    end;
    if parms.verbose
      toc;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

