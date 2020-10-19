function StudyInfo = MMIL_Import_BOLD_QC(StudyInfo,fname)
%function StudyInfo = MMIL_Import_BOLD_QC(StudyInfo,fname)
%
% Created:  12/31/12 by Don Hagler
% Last Mod: 12/31/12 by Don Hagler
%

%% todo: make BOLD_regQC optional (set to good if not included?)

if ~exist(fname,'file')
  fprintf('%s: file %s not found... no BOLD QC info will be imported\n',...
    mfilename,fname);
  return;
end;

nsubs = length(StudyInfo);

VisitIDs = {StudyInfo.VisitID};
BOLD_QC = repmat({''},nsubs,1);

%% todo: replace use of readtext_jcr with mmil_readtext
[boldqc_dat,res] = readtext_jcr(fname,',','','"');
%% todo: replace rmquotes with mmil_rmquotes
boldqc_dat = rmquotes(boldqc_dat);
boldqc_colnames = boldqc_dat(1,:);
boldqc_visID = boldqc_dat(2:end,strmatch('VisitID',boldqc_colnames,'exact'));
bold_qc =  boldqc_dat(2:end,strmatch('BOLDQC',boldqc_colnames,'exact'));
bold_regqc =  boldqc_dat(2:end,strmatch('BOLD_regQC',boldqc_colnames,'exact'));

for s=1:nsubs
  VisitID = StudyInfo(s).VisitID;
  ind = strmatch(VisitID,boldqc_visID,'exact');
  if isempty(ind), 
    fprintf('%s: WARNING: VisitID %s not found in %s\n',...
      mfilename,VisitID,fname);
    continue;
  end
  StudyInfo(s).BOLDQC = lower(bold_qc{ind(1)});
  StudyInfo(s).BOLD_regQC = lower(bold_regqc{ind(1)});
end

