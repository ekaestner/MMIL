function StudyInfo = MMIL_Import_FSRecon_QC(StudyInfo,fname)
%function StudyInfo = MMIL_Import_FSRecon_QC(StudyInfo,fname)
%
% Created:  02/05/09 by Alain Koyama
% Last Mod: 01/31/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname,'file')
  fprintf('%s: file %s not found. no FSRecon QC info will be imported\n',...
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
col_labels = setdiff(col_labels,{'VisitID'});

if ismember('QC',col_labels)
  qcflag = 1;
  col_labels = setdiff(col_labels,{'QC'});
else
  qcflag = 0;
  if any(~ismember({'ASEGQC','SurfaceQC'},col_labels))
    fprintf('%s: WARNING: %s missing QC info. mmust have QC column or ASEGQC and SurfaceQC columns. skipping import\n',...
      mfilename);
    return;
  end;
  col_labels = setdiff(col_labels,{'ASEGQC','SurfaceQC'});
end;

VisitIDs = {QCInfo.VisitID};
% make sure VisitIDs are strings
for s=1:length(VisitIDs)
  if isnumeric(VisitIDs{s})
    VisitIDs{s} = num2str(VisitIDs{s});
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StudyInfo(1).FSQC = [];

for s=1:length(StudyInfo)
  if isfield(StudyInfo,'STRUCT_VisitID') && ~isempty(StudyInfo(s).STRUCT_VisitID)
    VisitID = StudyInfo(s).STRUCT_VisitID;
  else
    VisitID = StudyInfo(s).VisitID;
  end;
  % find subject and visit in StudyInfo
  i=strmatch(VisitID,VisitIDs,'exact');
  if isempty(i), continue; end;
  % copy QCInfo values
  for j=1:length(col_labels)
    StudyInfo(s).(['FSQC_' col_labels{j}]) = QCInfo(i).(col_labels{j});
  end;  
  % set FSQC values
  if qcflag
    StudyInfo(s).FSQC = QCInfo(i).QC;
  else
    StudyInfo(s).AsegQC = QCInfo(i).ASEGQC;
    StudyInfo(s).SurfQC = QCInfo(i).SurfaceQC;
    if ~strcmp(StudyInfo(s).AsegQC,'bad') && ~strcmp(StudyInfo(s).SurfQC,'bad')
      StudyInfo(s).FSQC = 1;
    else
      StudyInfo(s).FSQC = 0;
    end;
  end;
end;

return;

