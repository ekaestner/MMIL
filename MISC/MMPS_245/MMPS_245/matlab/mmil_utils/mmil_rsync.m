function mmil_rsync(source,dest,varargin)
%function mmil_rsync(source,dest,[options])
%
% Purpose: use rsync to backup files
%
% Required Parameters:
%   source: source directory or cell array of multiple source directories
%     relative to rootdir
%   dest: full path of destination directory
%
% Optional Parameters:
%   'name': identifying name for set of files to be backed up
%     will be used to create subdirectory in dest and for log file
%     {default: []}
%   'rootdir': root directory containing source
%     if empty, source is assumed to be full path
%     {default: []}
%   'options': rsync command line options
%     {default: '-a --delete'}
%
% Created:  07/23/12 by Don Hagler
% Last Mod: 07/23/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'name',[],[],...
  'rootdir',[],[],...
  'options','-a --delete',[],...
});

if ~iscell(source), source = {source}; end;
ndirs = length(source);

if ~isempty(parms.name)
  fname_log = sprintf('%s/%s_bkup.log',dest,parms.name);
  dest = sprintf('%s/%s',dest,parms.name);
else
  fname_log = sprintf('%s/bkup.log',dest);
end;

mmil_mkdir(dest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:ndirs
  if ~isempty(parms.rootdir)
    sdir = sprintf('%s/%s',parms.rootdir,source{s});
  else
    sdir = source{s};
  end;
  ddir = dest;
  if ndirs==1
    [spath,sstem,sext] = fileparts(sdir);
    [dpath,dstem,dext] = fileparts(ddir);
    tmp_sdir = [sstem sext];
    tmp_ddir = [dstem dext];
    if strcmp(tmp_sdir,tmp_ddir)
      ddir = dpath;
    end;
  end;
  if ~exist(sdir,'file')
    fprintf('%s: WARNING: %s not found\n',mfilename,sdir);
  else
    cmd = sprintf('rsync %s %s %s',parms.options,sdir,ddir);
    [s,r] = mmil_unix(cmd);
    if s, error('cmd "%s" failed:\n%s',cmd,r); end;  
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fname_log,'wt');
if fid==-1
  error('failed to create log file %s',fname_log);
end;
if ~isempty(parms.name)
  fprintf(fid,'%s : ',parms.name);
end;
fprintf(fid,'last backup on %s\n',datestr(now,'mm/dd/yyyy at HH:MM:SS'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
