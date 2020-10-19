function StudyInfo = MMIL_Import_ViewerQC(StudyInfo,fname,projname)
%function StudyInfo = MMIL_Import_ViewerQC(StudyInfo,fname,projname)
%
% Created:  02/05/09 by Alain Koyama
% Last Mod: 03/24/11 by Don Hagler
%

% NOTE: MMIL_Import_ViewerQC now expects that subjectID is saved in
%       qcinfo.SubjID, not qcinfo.name. (as of 11/22/10 JCR)
%
% Note: a QC value of 1 or 2 in MMIL_Viewer becomes 1 (good), a value of 3 or 4 
% becomes -1 (bad), and 0 is still 0 (not done).. so basically it's the same as 
% how ADNIdvQC values are 'converted' when imported into ADNI StudyInfo


% determine what type of qc this is
[qcpath,qcfile] = fileparts(fname);
qctype = char(regexp(qcfile,['(?<=' projname '_).+(?=QC)'],'match'));

if ~exist(fname,'file')
  fprintf('%s: file %s not found... no %s QC info will be imported\n',...
    mfilename,fname,qctype);
  return;
end;

load(fname);

SubjIDs = {qcinfo.SubjID};
StudyDates = cell2mat({qcinfo.studydate});

if ~isfield(qcinfo,'input'), qcinfo(1).input = []; end;

for i=1:length(StudyInfo)
  SubjID = StudyInfo(i).SubjID;
  if isempty(SubjID) || strcmp('Unknown',SubjID), 
     fprintf('%s: WARNING: skipping subject with unknown or missing ID.\n', mfilename);
     continue; 
  end
  StudyDate = StudyInfo(i).StudyDate;
  if isempty(StudyDate), fprintf('%s: WARNING: skipping subject %s due to missing StudyDate\n', ...
     mfilename,SubjID); continue; 
  end
  
  % find subject and visit in qcinfo
  S_ind=find(strcmp(SubjIDs,SubjID) & StudyDates==StudyDate);
  if isempty(S_ind)
    fprintf('%s: WARNING: skipping SubjID %s StudyDate %d (not found in qcinfo)\n',...
      mfilename,SubjID,StudyDate);
    continue;
  end
  
  if length(S_ind) == 1 % single series QC
    if qcinfo(S_ind).qc == 1 || qcinfo(S_ind).qc == 2
      QCvalue = 1;
    elseif qcinfo(S_ind).qc == 3 || qcinfo(S_ind).qc == 4
      QCvalue = -1;
    else
      QCvalue = 0;
    end
    StudyInfo(i).([qctype 'QC']) = QCvalue;
    StudyInfo(i).([qctype 'QCNotes']) = qcinfo(S_ind).notes;
  else % multiple series QC
    QCvalues = cell2mat({qcinfo(S_ind).qc});
    n12 = find(QCvalues==1 | QCvalues==2);
    n34 = find(QCvalues==3 | QCvalues==4);
    if ~isempty(n12)
      StudyInfo(i).([qctype 'QC']) = 1;
    elseif ~isempty(n34)
      StudyInfo(i).([qctype 'QC']) = -1;
    else
      StudyInfo(i).([qctype 'QC']) = 0;
    end
    for j=1:length(S_ind) % concatenate comments from multiple series
      if j==1
        QCnotes = sprintf('%s: %s',char(qcinfo(S_ind(j)).input),qcinfo(S_ind(j)).notes);
      else
        QCnotes = sprintf('%s; %s: %s',QCnotes,char(qcinfo(S_ind(j)).input),qcinfo(S_ind(j)).notes);
      end
    end
    StudyInfo(i).([qctype 'QCNotes']) = QCnotes;
  end
end



