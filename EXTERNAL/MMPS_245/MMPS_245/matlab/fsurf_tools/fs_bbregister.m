function fs_bbregister(fname,subj,varargin)
%function fs_bbregister(fname,subj,[options])
%
% Purpose: a wrapper around FreeSurfer's bbregister (v4.3 and later)
%   for registering an MRI volume to a FreeSurfer recon
%
% Usage:
%  fs_bbregister(fname,subj,'key1', value1,...); 
%
% Required Parameters:
%   fname: full path name of volume to be registered
%   subj: FreeSurfer recon subject name
%
% Optional parameters:
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'fname_regdat': output name of registration file (e.g. register.dat)
%    If empty, will append '_register.dat' to file stem of fname_mov
%    {default: []}
%  'options' - bbregister options (see bbregister --help)
%    {default: '--init-fsl --t2'}
%  'forceflag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% Created:  06/11/09 by Josh Kuperman
% Last Mod: 08/25/09 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin, 2)), return; end;
parms = mmil_args2parms(varargin, { ...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'fname_regdat',[],[],...
  'options','--init-fsl --t2',[],...
  'binfile','bbregister',[],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subjpath = [parms.subjdir '/' subj];
if ~exist(subjpath,'dir')
  error('directory %s not found',subjpath);
end;
if ~exist(fname,'file')
  error('file %s not found',fname);
end;

% set fname_regdat
if isempty(parms.fname_regdat)
  [tmp_path,tmp_fstem] = fileparts(fname);
  parms.fname_regdat = [tmp_path '/' tmp_fstem '_register.dat'];
end;

% check if output exists already
if exist(parms.fname_regdat,'file') & ~parms.forceflag
  fprintf('%s: WARNING: not overwriting existing reg file %s\n',...
    mfilename,parms.fname_regdat);
  return;
end;

% set subjects dir
orig_subjdir = getenv('SUBJECTS_DIR');
setenv('SUBJECTS_DIR', parms.subjdir);

% create command
cmd = sprintf('bbregister --s %s --mov %s --reg %s %s\n',...
  subj,fname,parms.fname_regdat,parms.options);

% run registration
[status,result] = mmil_unix(cmd)
if ~status
  error('cmd: %s failed:\n%s',cmd,result);
else
  fprintf('%s\n%s\n',cmd,result);
end;

% restore the original subjects_dir
setenv('SUBJECTS_DIR',orig_subjdir);

