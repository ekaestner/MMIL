function errcode = MMIL_Save_ContainerInfo_ScanInfo(ContainerPath,ContainerInfo,ScanInfo)
%function errcode = MMIL_Save_ContainerInfo_ScanInfo(ContainerPath,ContainerInfo,ScanInfo)
%
% Required Input:
%   ContainerPath: full path of output container
%   ContainerInfo: container info struct
%   ScanInfo: scan info struct containing struct arrays for multiple scan types
%
% Created:  09/12/12 by Don Hagler
% Last Mod: 09/12/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,3) return; end;
errcode = 0;

ScanTypes = fieldnames(ScanInfo);
for s=1:length(ScanTypes)
  tag = ScanTypes{s};
  ContainerInfo.([tag '_cntr']) = length(ScanInfo.(tag));
end;
ContainerInfo.ScanInfo = ScanInfo;
ContainerInfo.Updated = datestr(now);
errcode = MMIL_Save_ContainerInfo(ContainerPath,ContainerInfo);

