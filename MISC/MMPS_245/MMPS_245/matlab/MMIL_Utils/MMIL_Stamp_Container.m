function [errcode,warncode] = MMIL_Stamp_Container(ContainerPath,ProjID)
%function [errcode,warncode] = MMIL_Stamp_Container(ContainerPath,ProjID)
%
% Purpose: add ProjID and MMPSVER to ContainerInfo
%
% Required Input:
%   ContainerPath: full path of data Container
%
% Optional Input:
%   ProjID: character string specifying Project ID
%
% Output:
%   errcode: returns 1 if error reading/writing ContainerInfo
%   warncode: returns 1 if warning about mismatched ProjID or MMPSVER
%
% Created:  12/26/10 by Don Hagler
% Last Mod: 04/12/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
errcode = 0;
warncode = 0;
if ~exist('ProjID','var'), ProjID = []; end;
MMPSVER = getenv('MMPSVER');

% load ContainerInfo
[ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode ~= 0, return; end;

% check ProjID
if isfield(ContainerInfo,'ProjID')
  if isempty(ContainerInfo.ProjID) && isempty(ProjID)
  elseif ~strcmp(ContainerInfo.ProjID,ProjID)
    fprintf('%s: WARNING: ProjID %s does not match ContainerInfo.ProjID %s\n',...
      mfilename,ProjID,ContainerInfo.ProjID);
    warncode = 1;
  end;
else
  ContainerInfo.ProjID = ProjID;
end;

% check MMPSVER
if isfield(ContainerInfo,'MMPSVER')
  if isempty(ContainerInfo.MMPSVER) && isempty(MMPSVER)
  elseif ~strcmp(ContainerInfo.MMPSVER,MMPSVER)
    fprintf('%s: WARNING: MMPSVER %s does not match ContainerInfo.MMPSVER %s\n',...
      mfilename,MMPSVER,ContainerInfo.MMPSVER);
    warncode = 1;
  end;
else
  ContainerInfo.MMPSVER = MMPSVER;
end;

% check for MMPS_nonsvn dir
user = getenv('USER');
nonsvn = ['/home/' user '/matlab/MMPS_nonsvn'];
if exist(nonsvn,'dir')
  dlist = dir([nonsvn '/*.m']);
  if ~isempty(dlist)
    ContainerInfo.MMPS_nonsvn = dlist;
  end;
else
  ContainerInfo.MMPS_nonsvn = [];
end;

% save modified ContainerInfo
errcode = MMIL_Save_ContainerInfo(ContainerPath,ContainerInfo);

