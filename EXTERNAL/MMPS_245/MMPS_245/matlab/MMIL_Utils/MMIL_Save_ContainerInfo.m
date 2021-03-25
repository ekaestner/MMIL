function errcode = MMIL_Save_ContainerInfo(ContainerPath,ContainerInfo)
%function errcode = MMIL_Save_ContainerInfo(ContainerPath,ContainerInfo)
%
% Required Input:
%   ContainerPath: full path of output container
%   ContainerInfo: container info struct
%
% Created:   12/26/10 by Don Hagler
% Last Mod:  09/12/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
errcode = 0;
fname = sprintf('%s/ContainerInfo.mat',ContainerPath);

try
  save(fname,'ContainerInfo');
catch
  fprintf('%s: ERROR: cannot save file %s:\n%s\n',mfilename,fname,lasterr);
  errcode = 1;
end;

