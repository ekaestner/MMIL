function StudyInfo = MMIL_Get_ContInfo(StudyInfo,RootDirs)
%function StudyInfo = MMIL_Get_ContInfo(StudyInfo,RootDirs)
%
% Purpose: load device and processing information from ContainerInfo
%   of each visit included in StudyInfo, located in RootDirs
%
% Required Parameters:
%  StudyInfo: struct array of study information (see MMIL_Get_StudyInfo)
%  RootDirs: struct specifying root directories
%
% Output:
%   StudyInfo: struct array with info for each "study" (scan session)
%
% Created:  10/25/12 by Don Hagler
% Prev Mod: 01/02/13 by Don Hagler
% Last Mod: 05/16/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

tags = {'Manufacturer','ManufacturersModelName',...
        'DeviceSerialNumber','MagneticFieldStrength'};

for s=1:length(StudyInfo)
  if isempty(mmil_getfield(StudyInfo(s),'raw'))
    [ContainerPath,StudyInfo(s).raw] = ...
      MMIL_Get_Container(RootDirs,StudyInfo(s).VisitID,'raw');
  else
    ContainerPath = [RootDirs.raw '/' StudyInfo(s).raw];
  end;

  if isempty(StudyInfo(s).raw) || ~exist(ContainerPath,'dir'), continue; end;
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
  if errcode, continue; end;
  for t=1:length(tags)
    tag = tags{t};
    StudyInfo(s).(tag) = mmil_getfield(ContainerInfo,tag);
  end;
  if isfield(ContainerInfo,'MMPSVER'),
    StudyInfo(s).MMPS_version = ContainerInfo.MMPSVER;
  end
  if isfield(ContainerInfo,'ClassificationCompleted'),
    StudyInfo(s).ProcDate = ...
      str2double(datestr(ContainerInfo.ClassificationCompleted,'yyyymmdd'));
  end
end;

