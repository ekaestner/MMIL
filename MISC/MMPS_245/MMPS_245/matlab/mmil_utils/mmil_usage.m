function usage = mmil_usage(rootdir,fname)
%function usage = mmil_usage([rootdir],[fname])
%
% Purpose: summarize the size of directories 
%
% Optional Parameters:
%   rootdir: root directory from which to descend two levels
%     {default: $HOME}
%   fname: output file name
%     {default: []}
%
% Created:  07/23/12 by Don Hagler
% Last Mod: 07/26/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rootdir','var') || isempty(rootdir)
  rootdir = getenv('HOME');
end;
if ~exist('fname','var')
  fname = [];
end;
usage = {};
if ~isempty(fname)
  fid = fopen(fname,'wt');
else
  fid = 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1;
dlist = dir(rootdir);
for d=1:length(dlist)
  dname = dlist(d).name;
  full_dname = sprintf('%s/%s',rootdir,dname);
  if strncmp(dname,'.',1) || ~dlist(d).isdir || strcmp(dname,'lost+found')
    continue;
  end;
  [status,result] = mmil_unix(sprintf('du -sh %s',full_dname));
  result = deblank(result);
  [usage,n] = display_usage(result,fid,usage,n,0);
  slist = dir(full_dname);
  for s=1:length(slist)
    sname = slist(s).name;
    full_sname = sprintf('%s/%s',full_dname,sname);
    if strncmp(sname,'.',1) || ~slist(s).isdir || strcmp(dname,'lost+found')
      continue;
    end;
    [status,result] = mmil_unix(sprintf('du -sh %s',full_sname));
    [usage,n] = display_usage(result,fid,usage,n,1);
  end;
end;

if ~isempty(fname), fclose(fid); end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [usage,n] = display_usage(result,fid,usage,n,subflag)
  tmp = regexp(result,...
    '(?<size>[\d\w\.]+)\s+(?<location>[\w\/\.]+)','names','once');
  if isempty(tmp), return; end;
  usage{1,n} = tmp.size;
  usage{2,n} = tmp.location;
  n = n + 1;
  if ~subflag
    fprintf(fid,'%-6s %s\n',tmp.size,tmp.location);
  else
    fprintf(fid,'%-12s %s\n',tmp.size,tmp.location);
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

