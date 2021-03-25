function StudyInfo = MMIL_Import_DTI_QC(StudyInfo,fname)
%function StudyInfo = MMIL_Import_DTI_QC(StudyInfo,fname)
%
% Created:  11/07/12 by Don Hagler
% Prev Mod: 01/29/16 by Don Hagler
% Last Mod: 05/30/17 by Don Hagler
%

if ~exist(fname,'file')
  fprintf('%s: file %s not found. no DTI QC info will be imported\n',...
    mfilename,fname);
  return;
end;

QCInfo = mmil_csv2struct(fname);

col_labels = fieldnames(QCInfo);

if ~ismember('VisitID',col_labels)
  fprintf('%s: WARNING: %s missing VisitID column. skipping import\n',...
    mfilename,fname);
  return;
end;
col_labels = setdiff(col_labels,'VisitID');

if ismember('QC',col_labels)
  qcflag = 1;
  col_labels = setdiff(col_labels,'QC');
else
  qcflag = 0;
  if any(~ismember({'DTIQC','DTI_regQC'},col_labels))
    fprintf('%s: WARNING: %s missing QC info. mmust have QC column or DTIQC and DTI_regQC columns. skipping import\n',...
      mfilename);
    return;
  end;
  col_labels = setdiff(col_labels,{'DTIQC','DTI_regQC'});
end;

VisitIDs = {StudyInfo.VisitID};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StudyInfo(1).DTIQC = [];

for i=1:length(QCInfo)
  VisitID = QCInfo(i).VisitID;
  if isnumeric(VisitID),
     VisitID = num2str(VisitID);
  end
  % find subject and visit in StudyInfo
  s=strmatch(VisitID,VisitIDs,'exact');
  if isempty(s), continue; end;
  % copy QCInfo values
  for j=1:length(col_labels)
    StudyInfo(s).(['DTIQC_' col_labels{j}]) = QCInfo(i).(col_labels{j});
  end;  
  % set DTI QC values
  if qcflag
    StudyInfo(s).DTIQC = QCInfo(i).QC;
  else
    StudyInfo(s).DTIQC = QCInfo(i).DTIQC;
    StudyInfo(s).DTI_regQC = QCInfo(i).DTI_regQC;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

