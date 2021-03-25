function StudyInfo = MMIL_Import_Meta(StudyInfo,fname)
%function StudyInfo = MMIL_Import_Meta(StudyInfo,fname)
%
% Created:  02/10/09 by Alain Koyama
% Last Mod: 03/24/11 by Don Hagler
%

if ~exist(fname,'file')
  fprintf('%s: file %s not found.. no Subject MetaData info will be imported\n',...
    mfilename,fname);
  return;
end;

Meta = mmil_readtext(fname,',','','"');

SubjIDs = {StudyInfo.SubjID};
StudyDates = [StudyInfo.StudyDate];

for i=2:size(Meta,1) % loop thru qc spreadsheet and update StudyInfo
  % pull directory name from first column and extract Subject ID, Study Date
  SubjID = Meta{i,1};
  if isnumeric(Meta{i,2})
    StudyDate = Meta{i,2};
  else 
    StudyDate = str2num(Meta{i,2});
  end
  ncols = size(Meta,2);

  % find subject and visit in StudyInfo
  S_ind=find(strcmp(SubjIDs,SubjID) & StudyDates==StudyDate);

  if isempty(S_ind)
    fprintf('%s: WARNING: skipping SubjID %s StudyDate %d (not found in StudyInfo)\n',...
      mfilename,SubjID,StudyDate);
    continue;
  end;
  if length(S_ind) > 1
    fprintf('%s: WARNING: skipping SubjID %s StudyDate %d (duplicates found in StudyInfo)\n',...
      mfilename,SubjID,StudyDate);
    continue;
  end;

  if ncols < 3
    fprintf('%s: WARNING: No MetaData columns found in file %s\n',...
      mfilename,fname);
    return;
  end

  for j=3:ncols
    MetaCol = Meta{1,j}; % grab header name for MetaData
    % MetaCol = regexprep(MetaCol,'^[^a-zA-Z]','MetaData_'); % first character of field must be a letter so it can be a struct field
    MetaCol = regexprep(MetaCol,'\W','_'); % remove invalid characters so it can be a struct field
    StudyInfo(S_ind).(MetaCol) = Meta{i,j};
  end
end

