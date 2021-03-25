function [StudyInfo,RootDirs,ProjInfo]=REC_MMIL_Get_StudyInfo(ProjID,varargin)
%function [StudyInfo,RootDirs,ProjInfo]=REC_MMIL_Get_StudyInfo(ProjID,[options])
%
%  Usage: [StudyInfo,RootDirs,ProjInfo]=...
%               REC_MMIL_Get_StudyInfo(ProjID,'key',value,...)
%
% Required Parameters:
%   ProjID: project ID string
%
% Optional Parameters that determine which studies have overall QC=1:
%   'QC_raw' : [0|1] use good raw QC if it exists for this project
%     {default = 1}
%   'QC_recon' : [0|1] use recon QC if it exists for this project
%     {default = 1}
%   'QC_dv' : [0|1] use dv QC if it exists for this project
%     {default = 1}
%   'QC_DTI' : [0|1] use DTI QC if it exists for this project
%     {default = 1}
%   'QC_PET' : [0|1] use PET QC if it exists for this project
%     {default = 1}
%
% Output:
%   StudyInfo: struct array with info for each "study" (scan session)
%   RootDirs: struct containing locations of data directories
%   ProjInfo: struct containing parameters for this ProjID
%
% Created:  07/21/09 by Don Hagler
% Last Mod: 11/08/12 by Don Hagler
%

%% todo: remove this function and transition completely to use of MMIL_Check_ProjID

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin,{...
  'QC_raw',true,[false true],...
  'QC_recon',true,[false true],...
  'QC_dv',true,[false true],...
  'QC_DTI',true,[false true],...
  'QC_PET',true,[false true],...
  'metaflag',true,[false true],...
...
  'PETflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StudyInfo = [];
RootDirs = [];

if (~mmil_check_nargs(nargin,1)) return; end;
if isempty(ProjID)
  error('must specify a ProjID');
end;
if iscell(ProjID) & length(ProjID)>1
  error('may only specify one ProjID');
elseif iscell(ProjID)
  ProjID = ProjID{1};
end;

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);
if isfield(RootDirs,'raw_pet') & ~isempty(RootDirs.raw_pet), parms.PETflag = 1; end;

ProjInfoDir = sprintf('%s/ProjInfo/%s',RootDirs.home,ProjID);

% set input file names
fname_visit = sprintf('%s/%s_VisitInfo.csv',ProjInfoDir,ProjID); % file with visit numbers (for longitudinal)
fname_meta = sprintf('%s/REC_%s_Meta.csv',ProjInfoDir,ProjID); % file with demographic info, group info, neuropsych scores.. etc
fname_pet_qc = sprintf('%s/REC_%s_PET_QC.csv',ProjInfoDir,ProjID); % PET QC, manually updated
fname_reconqc = sprintf('%s/%s_FSReconQC.csv',ProjInfoDir,ProjID); % MRI QC, manually updated
fname_dvqc = sprintf('%s/%s_dvQC.mat',ProjInfoDir,ProjID); % dv QC, updated using MMIL_Viewer
fname_rawqc = sprintf('%s/%s_RawQC.mat',ProjInfoDir,ProjID); % Raw QC, updated using MMIL_Viewer
fname_dtiqc = sprintf('%s/%s_DTIQC.csv',ProjInfoDir,ProjID); % DTI Reg QC, updated manually

% check for existence of Subject MetaData, PET QC, Recon QC, dv QC, and Raw QC
if ~exist(fname_meta,'file')
  parms.metaflag = 0;
end;
if ~parms.PETflag
  parms.QC_PET = 0;
end;
if ~exist(fname_reconqc,'file')
  parms.QC_recon = 0;
end;
if ~exist(fname_dvqc,'file')
  parms.QC_dv = 0;
end;
if ~isfield(ProjInfo,'STRUCT_rawQCflag') | isempty(ProjInfo.STRUCT_rawQCflag) |...
 ~ProjInfo.STRUCT_rawQCflag | ~exist(fname_rawqc,'file')
  parms.QC_raw = 0;
end;
if ~exist(fname_dtiqc,'file')
  parms.QC_DTI = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load visit info if exists
if exist(fname_visit,'file')
  StudyInfo = MMIL_Read_StudyInfo(fname_visit);
end
% find containers for each visit, or create studyinfo based on raw containers
[StudyInfo,RootDirs] = MMIL_Check_StudyInfo(StudyInfo,RootDirs);

% if empty, project has not been unpacked
if isempty(StudyInfo), return; end;

if ~isfield(StudyInfo,'VisitNumber')
  for i=1:length(StudyInfo)
    StudyInfo(i).VisitNumber = 1;
  end;
end;

% get list of PET dirs
if parms.PETflag
  dirlist_PET = dir(sprintf('%s/PETRAW_*',RootDirs.raw_pet));
end

% initialize extra fields
fieldnames = {'StudyDate_PET',...
  'RawQC','AsegQC','AsegQCNotes',...
  'SurfQC','SurfQCNotes','PETRegQC','PETRegQCNotes'};
for f=1:length(fieldnames)
  if ~isfield(StudyInfo,fieldnames{f})
    StudyInfo = setfield(StudyInfo,{1},fieldnames{f},[]);
  end;
end;

%% todo: is this necessary?  should be in MMIL_Check_StudyInfo
% add extra info about PET (needed when MRI and PET dates do not match)
for i=1:length(StudyInfo)
  StudyInfo(i).StudyDate_PET = [];
  if parms.PETflag
    if isempty(StudyInfo(i).raw_pet)
      n = regexp({dirlist_PET.name},['PETRAW_' StudyInfo(i).VisitID '_\d{8}.+'],'match');
      n(cellfun(@isempty,n)) = [];
      if ~isempty(n) % if this particular subject has PET
        StudyInfo(i).StudyDate_PET = str2num(char(regexp(char(n{1}),'(?<=_)\d{8}(?=\.)','match')));
        StudyInfo(i).raw_pet = char(n{1});
      end
    end;
    StudyInfo(i).StudyDate_PET = str2num(char(regexp(char(StudyInfo(i).raw_pet),'(?<=_)\d{8}(?=\.)','match')));
    if isempty(StudyInfo(i).proc_pet) & ~isempty(RootDirs.proc_pet)
      [tmp,StudyInfo(i).proc_pet] = MMIL_Get_Container(RootDirs,...
        StudyInfo(i).VisitID,'proc_pet');
    end;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import QC and meta info
if parms.QC_PET % import existing PET QC info from spreadsheet
  fprintf('%s: importing PET QC info from %s...\n',mfilename,fname_pet_qc);
  StudyInfo = MMIL_Import_PET_QC(StudyInfo,fname_pet_qc);
end
if parms.QC_recon % import existing Recon QC info from spreadsheet
  fprintf('%s: importing Recon QC info from %s...\n',mfilename,fname_reconqc);
  StudyInfo = MMIL_Import_MRI_QC(StudyInfo,fname_reconqc);
end
if parms.QC_dv  % import existing dvQC info from .mat file (do not need to export, as MMIL_Viewer will update)
  fprintf('%s: importing dv QC info from %s...\n',mfilename,fname_dvqc);
  StudyInfo = MMIL_Import_ViewerQC(StudyInfo,fname_dvqc,ProjID);
end
if parms.QC_raw % import existing rawQC info from .mat file (do not need to export, as MMIL_Viewer will update)
  fprintf('%s: importing Raw QC info from %s...\n',mfilename,fname_rawqc);
  StudyInfo = MMIL_Import_ViewerQC(StudyInfo,fname_rawqc,ProjID);
end
if parms.QC_DTI % import existing DTI QC info from csv spreadsheet file
  fprintf('%s: importing DTI Reg QC info from %s...\n',mfilename,fname_dti_qc);
  StudyInfo = MMIL_Import_DTI_QC(StudyInfo,fname_dti_qc,ProjID);
end
if parms.metaflag % import Subject MetaData from spreadsheet (group info, neuropsych scores, etc..)
  fprintf('%s: importing Subject MetaData from %s...\n',mfilename,fname_meta);
  StudyInfo = MMIL_Import_Meta(StudyInfo,fname_meta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% overall QC
% loop thru all QC's, if any are incomplete, then QC = 0; if any are bad, then QC = -1
for i=1:length(StudyInfo)
  StudyInfo(i).QC = 1; % overall QC is good unless we find otherwise below

  if parms.QC_PET && isempty(StudyInfo(i).PETRegQC)
    StudyInfo(i).QC_PETReg = 0;
    StudyInfo(i).QC = 0;
  elseif parms.QC_PET && ~isempty(regexpi(StudyInfo(i).PETRegQC,'bad'))
    StudyInfo(i).QC_PETReg = -1;
    StudyInfo(i).QC = -1;
  else
    StudyInfo(i).QC_PETReg = 1;
  end

  if parms.QC_raw && (isempty(StudyInfo(i).RawQC) || StudyInfo(i).RawQC==0)
    StudyInfo(i).RawQC = 0;
    StudyInfo(i).QC = 0;
  elseif parms.QC_raw && StudyInfo(i).RawQC==-1
    StudyInfo(i).QC = -1;
  end

  if parms.QC_recon && isempty(StudyInfo(i).AsegQC)
    StudyInfo(i).QC_ASEG = 0;
    StudyInfo(i).QC = 0;
  elseif parms.QC_recon && ~isempty(regexpi(StudyInfo(i).AsegQC,'bad'))
    StudyInfo(i).QC = -1;
    StudyInfo(i).QC_ASEG = -1;
  else
    StudyInfo(i).QC_ASEG = 1;
  end

  if parms.QC_recon && isempty(StudyInfo(i).SurfQC)
    StudyInfo(i).QC_surf = 0;
    StudyInfo(i).QC = 0;
  elseif parms.QC_recon && ~isempty(regexpi(StudyInfo(i).SurfQC,'bad'))
    StudyInfo(i).QC = -1;
    StudyInfo(i).QC_surf = -1;
  else
    StudyInfo(i).QC_surf = 1;
  end

  if parms.QC_DTI && isempty(StudyInfo(i).DTIQC)
    StudyInfo(i).QC_DTI = 0;
    StudyInfo(i).QC_DTIreg = 0
    StudyInfo(i).QC = 0;
  elseif parms.QC_DTI
    if ~isempty(regexpi(StudyInfo(i).DTIQC,'bad'))
      StudyInfo(i).QC_DTI = -1;
    end;
    if ~isempty(regexpi(StudyInfo(i).DTI_regQC,'bad'))
      StudyInfo(i).QC_DTIreg = -1;
    end;
    if any([StudyInfo(i).QC_DTI,StudyInfo(i).QC_DTIreg]==-1)
      StudyInfo(i).QC = -1;
    end;
  else
    StudyInfo(i).QC_DTI = 1;
  end
end

