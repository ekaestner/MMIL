function StudyInfo = MMIL_Import_PET_QC(StudyInfo,fname)
%function StudyInfo = MMIL_Import_PET_QC(StudyInfo,fname)
%
% Created:  11/24/08 by Alain Koyama
% Last Mod: 03/23/11 by Don Hagler
%

if ~exist(fname,'file')
  fprintf('mfilename: file %s not found.. no PET QC info will be imported\n',...
    mfilename,fname);
  return;
end;

PET_QC = mmil_readtext(fname,',','','"');

% Find relevant columns in spreadsheet
headers = PET_QC(1,:);
n = regexpi(headers,'^SubjID$'); n = ~cellfun('isempty', n);
SubjID_col = find(n==1);
n = regexpi(headers,'^StudyDate_PET$'); n = ~cellfun('isempty', n);
StudyDate_col = find(n==1);
n = regexpi(headers,'^PETRegQC$'); n = ~cellfun('isempty', n);
PETRegQC_col = find(n==1);
n = regexpi(headers,'^PETRegQCNotes$'); n = ~cellfun('isempty', n);
PETRegQCNotes_col = find(n==1);

% Find SubjectID and StudyDate in StudyInfo struct
SubjIDs = {StudyInfo.SubjID};
StudyDates = cell2mat({StudyInfo.StudyDate_PET});

% remove those subjects w/o a PET study date
SubjIDs = SubjIDs(~cellfun('isempty', {StudyInfo.StudyDate_PET})==1);

for i=2:size(PET_QC,1)
  % pull relevant values from spreadsheet
  SubjID = PET_QC{i,SubjID_col};
  StudyDate = PET_QC{i,StudyDate_col};
  PETRegQC = PET_QC{i,PETRegQC_col};
  PETRegQCNotes = PET_QC{i,PETRegQCNotes_col};

  if isempty(StudyDate), continue; end; % skip if no PET for this subject
  
  S_ind=find(strcmp(SubjIDs,SubjID) & StudyDates==StudyDate);
  
  if isempty(S_ind)
    fprintf('%s: WARNING: skipping SubjID %s StudyDate %d (not found in StudyInfo)\n',...
      mfilename,SubjID,StudyDate);
    continue;
  end;
  if length(S_ind) > 1
    fprintf('%s: WARNING: multiple studies with SubjID %s\n',...
      mfilename,SubjID);
    S_ind = S_ind(1);
  end;

  StudyInfo(S_ind).PETRegQC = PETRegQC;
  StudyInfo(S_ind).PETRegQCNotes = PETRegQCNotes;
end

