function StudyInfo = MMIL_Get_SubjInfo(StudyInfo,RootDirs,ProjID)
%function StudyInfo = MMIL_Get_SubjInfo(StudyInfo,RootDirs,ProjID)
%
% Purpose: load subject specific information and attach it to
%   each matching visit included in StudyInfo
%
% Required Parameters:
%  StudyInfo: struct array of study information (see MMIL_Get_StudyInfo)
%  RootDirs: struct specifying root directories
%  ProjID: project ID string
%
% Notes:
%   for this function to work, this file must exist:
%     {RootDirs.home}/ProjInfo/{ProjID}/{ProjID}_SubjInfo.csv
%   one column must have SubjID as the column header
%     and contain subject ID strings that match those in StudyInfo
%   all other column headers will be used to create fields in StudyInfo struct
%     e.g. Sex, Site, DOB, Group
%   if DOB is supplied, 'Age' will be calculated based on StudyDate
%     should have 'yyyymmdd' or 'yyyymm' format
%     if 'yyyymm' format is used, the day will be set to 15
%
% Output:
%   StudyInfo: struct array with info for each "study" (scan session)
%
% Created:  10/25/12 by Don Hagler
% Last Mod: 04/23/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

fname_info = [RootDirs.home '/ProjInfo/' ProjID '/' ProjID '_SubjInfo.csv'];
if ~exist(fname_info,'file')
  fprintf('%s: WARNING: SubjInfo file %s not found\n',mfilename,fname_info);
  return;
end;

SubjInfo = MMIL_Read_StudyInfo(fname_info);

tags = fieldnames(SubjInfo);
if ~ismember('SubjID',tags)
  fprintf('%s: WARNING: SubjID column not found in %s\n',mfilename,fname_info);
  return;
end;
SubjIDs = {SubjInfo.SubjID};

for s=1:length(StudyInfo)
  SubjID = StudyInfo(s).SubjID;
  ind = find(strcmp(SubjID,SubjIDs));
  if isempty(ind)
    fprintf('%s: WARNING: SubjID %s not found in SubjInfo\n',...
      mfilename,SubjID);
    continue;
  end;
  for t=1:length(tags)
    tag = tags{t};
    StudyInfo(s).(tag) = mmil_getfield(SubjInfo(ind),tag);
  end;

  % calculate Age from DOB and StudyDate
  if ~isfield(StudyInfo,'Age')
    StudyInfo(s).Age = [];
  end;
  DOB = mmil_getfield(StudyInfo(s),'DOB');
  StudyDate = mmil_getfield(StudyInfo(s),'StudyDate');
  if ~isempty(DOB) && ~isempty(StudyDate)
    if isnumeric(DOB), DOB = num2str(DOB); end;
    if isnumeric(StudyDate), StudyDate = num2str(StudyDate); end;
    switch length(DOB)
      case 4, DOB = [DOB '0630'];
      case 6, DOB = [DOB '15'];
      case 8
      otherwise
        fprintf('%s: WARNING: invalid DOB "%s" for %s\n',...
          mfilename,DOB,SubjID);
    end
    StudyInfo(s).DOB = str2double(DOB);
    try 
      StudyInfo(s).Age = ...
        (datenum(StudyDate,'yyyymmdd') - datenum(DOB,'yyyymmdd'))/365.25;
    catch 
      fprintf('%s: WARNING: unable to calculate age for %s\n',...
        mfilename,StudyInfo(s).VisitID);
    end
  end; 
end;

