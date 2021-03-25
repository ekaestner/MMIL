function mmil_bbregister(fname,FSpath,varargin)
%function mmil_bbregister(fname,FSpath,[options])
%
% Purpose: a wrapper around FreeSurfer's bbregister (v4.3 and higher)
%   for registering an MRI volume to a FreeSurfer recon
%
% Usage:
%  mmil_bbregister(FSpath, fname,'key1', value1,...); 
%
% Required Parameters:
%   fname: full path name of volume to be registered
%   FSpath: full path name of FreeSurfer recon
%
% Optional parameters:
%  'fname_regdat': output name of registration file (e.g. register.dat)
%    If empty, will append '_register.dat' to file stem of fname_mov
%    {default: []}
%  'options' - bbregister options (see bbregister --help)
%    {default: '--init-fsl --t2'}
%  'FS_version': FreeSurfer version number (e.g. 430, 450, etc.)
%    {default: 450}
%  'forceflag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% Created:  08/19/09 by Don Hagler
% Rcnt Mod: 04/21/10 by Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin, 2)), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_regdat',[],[],...
  'options','--init-fsl --t2',[],...
  'FS_version',450,[],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get subjdir and subj name
[subjdir,tmp_stem,tmp_ext] = fileparts(FSpath);
subj = [tmp_stem,tmp_ext];

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
setenv('SUBJECTS_DIR', subjdir);

% create command
cmd = sprintf('setenv SUBJECTS_DIR %s\n',subjdir);
cmd = sprintf('%ssource $PUBSH/bin/SetUpFreeSurfer.csh %d\n',...
  cmd,parms.FS_version);
cmd = sprintf('%sbbregister --s %s --mov %s --reg %s %s\n',...
  cmd,subj,fname,parms.fname_regdat,parms.options);

% run registration
[status,result] = mmil_unix(cmd);
if status
  error('cmd: %s failed:\n%s',cmd,result);
else
  fprintf('%s\n%s\n',cmd,result);
end;

% restore the original subjects_dir
setenv('SUBJECTS_DIR',orig_subjdir);

