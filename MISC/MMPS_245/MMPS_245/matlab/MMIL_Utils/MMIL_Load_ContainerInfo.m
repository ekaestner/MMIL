function [ContainerInfo, errcode] = MMIL_Load_ContainerInfo(ContainerPath)
%function [ContainerInfo, errcode] = MMIL_Load_ContainerInfo(ContainerPath)
%
% Early Mod:  10/16/09 by Don Hagler
% Last Mod:   02/26/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
ContainerInfo = [];
errcode = 0;

fname = sprintf('%s/ContainerInfo.mat',ContainerPath);
fname_dir = dir(fname);
if ~exist(fname,'file')
  fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
  errcode = 1;
  return;
elseif fname_dir.bytes == 0
  fprintf('%s: ERROR: file %s is empty (0 bytes)\n',mfilename,fname);
  return;
end
try
  load(fname);
catch
  fprintf('%s: ERROR: cannot read file %s:\n%s\n',mfilename,fname,lasterr);
  errcode = 1;
end

