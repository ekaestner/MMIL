function [FS_version,errcode] = MMIL_Get_FS_version(ContainerPath,default_version)
%function [FS_version,errcode] = MMIL_Get_FS_version(ContainerPath,[default_version])
%
% Purpose: loads ContainerInfo and gets ContainerInfo.FS_version
%
% Required input:
%   ContainerPath: full path of freesurfer recon directory (for one subject)
%
% Optional input:
%   default_version: value to return if ContainerInfo
%     or FS_version field not found
%     {default = []}
%
% Output:
%  FS_version: freesurfer version number (e.g. 305, 450, etc.)
%
% Created:  11/04/11 Don Hagler
% Last Mod: 11/04/11 Don Hagler
%

% check parameters
if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('default_version','var'), default_version = []; end;
FS_version = default_version;
errcode = 0;

[ContainerInfo, errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode, return; end;

FS_version = mmil_getfield(ContainerInfo,'FS_version',default_version);


