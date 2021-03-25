function StudyInfo = MMIL_Read_StudyInfo(fname)
%function StudyInfo = MMIL_Read_StudyInfo(fname)
%
%  fname: file name of csv file containing information about each study
%           each column will become a field in the StudyInfo output struct
%           each row is from a different subject or exam
%
% Created:  03/31/09 by Don Hagler
% Last Mod: 05/02/12 by Don Hagler
%

StudyInfo = [];
if (~mmil_check_nargs(nargin,1)) return; end;
if ~exist(fname,'file'), error('file %s not found',fname); end;

% create struct array with field names from the column headers of csv file
StudyInfo = mmil_csv2struct(fname);

% ignore lines that start with % or empty SubjID/VisitID
ID = [];
if isfield(StudyInfo,'SubjID')
  ID = 'SubjID';
elseif isfield(StudyInfo,'VisitID')
  ID = 'VisitID';
end;
if ~isempty(ID)
  ind_keep = [];
  for i=1:length(StudyInfo)
    if ~isempty(StudyInfo(i).(ID)) & ~strncmp(StudyInfo(i).(ID),'%',1)
      ind_keep = [ind_keep, i];
    end;
  end;
  StudyInfo = StudyInfo(ind_keep);
end;

