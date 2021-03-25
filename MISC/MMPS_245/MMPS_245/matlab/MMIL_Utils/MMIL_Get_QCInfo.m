function StudyInfo = MMIL_Get_QCInfo(StudyInfo,RootDirs,ProjInfo,varargin)
%function StudyInfo = MMIL_Get_QCInfo(StudyInfo,RootDirs,ProjInfo,[options])
%
% Required Parameters:
%  StudyInfo: struct array of study information
%  RootDirs: struct specifying root directories
%  ProjInfo: struct containing project-specific parameters
%    NOTE: these input parameters are loaded and checked in MMIL_Get_StudyInfo
%
% Optional Parameters that determine which studies have overall QC=1:
%   'QC_raw' : [0|1] use good raw QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_RawQC.mat
%     {default = 1}
%   'QC_dv' : [0|1] use dv QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_dvQC.mat
%     {default = 1}
%   'QC_recon' : [0|1] use recon QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_FSReconQC.csv
%     {default = 1}
%   'QC_DTI' : [0|1] use DTIQC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_DTIQC.csv
%     {default = 1}
%   'QC_RSI' : [0|1] use RSIQC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_RSIQC.csv
%     {default = 1}
%   'QC_BOLD' : [0|1] use BOLDQC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_BOLDQC.csv
%     {default = 1}
%   'QC_PET' : [0|1] use PET QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_PET_QC.csv
%     {default = 1}
%   'metaflag': [0|1] use meta data csv file if found
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_Meta.csv
%     {default = 1}
%
% Output:
%   StudyInfo: struct array with info for each "study" (scan session)
%
% Created:  03/18/11 by Don Hagler
% Prev Mod: 01/27/16 by Don Hagler
% Last Mod: 05/30/17 by Don Hagler
%

% based on REC_MMIL_Get_StudyInfo, created 07/21/09 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: use csv file for all QC (need to change MMIL_Viewer or not use)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms(varargin,{...
  'QC_raw',true,[false true],...
  'QC_recon',true,[false true],...
  'QC_dv',true,[false true],...
  'QC_DTI',true,[false true],...
  'QC_RSI',true,[false true],...
  'QC_BOLD',true,[false true],...
  'QC_PET',true,[false true],...
  'metaflag',true,[false true],...
});

ProjID = ProjInfo.ProjID;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load files with QC information
indir = sprintf('%s/ProjInfo/%s',RootDirs.home,ProjID);
% file with demographic info, group info, neuropsych scores.. etc
fname_meta = sprintf('%s/%s_Meta.csv',indir,ProjID);
% Raw QC, updated using MMIL_Viewer
fname_raw_qc = sprintf('%s/%s_RawQC.mat',indir,ProjID);
% dv QC, updated using MMIL_Viewer
fname_dv_qc = sprintf('%s/%s_dvQC.mat',indir,ProjID);
% PET QC, manually updated
fname_pet_qc = sprintf('%s/%s_PET_QC.csv',indir,ProjID);
% FreeSurfer QC, manually updated
fname_recon_qc = sprintf('%s/%s_FSReconQC.csv',indir,ProjID);
% DTI QC, manually updated
fname_dti_qc = sprintf('%s/%s_DTIQC.csv',indir,ProjID);
% RSI QC, manually updated
fname_rsi_qc = sprintf('%s/%s_RSIQC.csv',indir,ProjID);
% BOLD QC, manually updated
fname_bold_qc = sprintf('%s/%s_BOLDQC.csv',indir,ProjID);

% check for existence of Subject MetaData, PET QC, Recon QC, dv QC, and Raw QC
if ~exist(fname_meta,'file'), parms.metaflag = 0; end;
if ~exist(fname_pet_qc,'file'), parms.QC_PET = 0; end;
if ~exist(fname_recon_qc,'file'),parms.QC_recon = 0; end;
if ~exist(fname_dv_qc,'file'), parms.QC_dv = 0; end;
rawQCflag = mmil_getfield(ProjInfo,'STRUCT_rawQCflag',0);
if ~rawQCflag || ~exist(fname_raw_qc,'file')
  parms.QC_raw = 0;
end;
if ~exist(fname_dti_qc,'file'),parms.QC_DTI = 0; end;
if ~exist(fname_rsi_qc,'file'),parms.QC_RSI = 0; end;
if ~exist(fname_bold_qc,'file'),parms.QC_BOLD = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if empty, project has not been unpacked
if isempty(StudyInfo), return; end;

if ~isfield(StudyInfo,'VisitNumber')
  for i=1:length(StudyInfo)
    StudyInfo(i).VisitNumber = 1;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import QC and meta info
if parms.QC_PET % import existing PET QC info from spreadsheet
  fprintf('%s: importing PET QC info from %s...\n',mfilename,fname_pet_qc);
  StudyInfo = MMIL_Import_PET_QC(StudyInfo,fname_pet_qc);
end
if parms.QC_recon % import existing Recon QC info from csv spreadsheet file
  fprintf('%s: importing Recon QC info from %s...\n',mfilename,fname_recon_qc);
  StudyInfo = MMIL_Import_FSRecon_QC(StudyInfo,fname_recon_qc);
end
if parms.QC_dv  % import existing dvQC info from .mat file (do not need to export, as MMIL_Viewer will update)
  fprintf('%s: importing dv QC info from %s...\n',mfilename,fname_dv_qc);
  StudyInfo = MMIL_Import_ViewerQC(StudyInfo,fname_dv_qc,ProjID);
end
if parms.QC_raw % import existing rawQC info from .mat file (do not need to export, as MMIL_Viewer will update)
  fprintf('%s: importing Raw QC info from %s...\n',mfilename,fname_raw_qc);
  StudyInfo = MMIL_Import_ViewerQC(StudyInfo,fname_raw_qc,ProjID);
end
if parms.QC_DTI % import existing DTI QC info from csv spreadsheet file
  fprintf('%s: importing DTI QC info from %s...\n',mfilename,fname_dti_qc);
  StudyInfo = MMIL_Import_DTI_QC(StudyInfo,fname_dti_qc);
end
if parms.QC_RSI % import existing RSI QC info from csv spreadsheet file
  fprintf('%s: importing RSI QC info from %s...\n',mfilename,fname_rsi_qc);
  StudyInfo = MMIL_Import_RSI_QC(StudyInfo,fname_rsi_qc);
end
if parms.QC_BOLD % import existing BOLD QC info from csv spreadsheet file
  fprintf('%s: importing BOLD QC info from %s...\n',mfilename,fname_bold_qc);
  StudyInfo = MMIL_Import_BOLD_QC(StudyInfo,fname_bold_qc);
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

  if parms.QC_PET
    StudyInfo(i).QC_PETReg = 1;
    if isempty(StudyInfo(i).PETRegQC)
      StudyInfo(i).QC_PETReg = 0;
      StudyInfo(i).QC = 0;
    elseif ~isempty(regexpi(StudyInfo(i).PETRegQC,'bad'))
      StudyInfo(i).QC_PETReg = -1;
      StudyInfo(i).QC = -1;
    end
  end;

  if parms.QC_raw
    StudyInfo(i).QC_raw = 1;
    if (isempty(StudyInfo(i).RawQC) || StudyInfo(i).RawQC==0)
      StudyInfo(i).QC_raw = 0;
      StudyInfo(i).QC = 0;
    elseif StudyInfo(i).RawQC==-1
      StudyInfo(i).QC_raw = -1;
      StudyInfo(i).QC = -1;
    end
  end;

  if parms.QC_recon
    if isfield(StudyInfo,'FSQC')
      StudyInfo(i).QC_FSRecon = 1;
      if isempty(StudyInfo(i).FSQC)
        StudyInfo(i).QC_FSRecon = 0;
        StudyInfo(i).QC = 0;
      elseif ~StudyInfo(i).FSQC
        StudyInfo(i).QC_FSRecon = -1;
        StudyInfo(i).QC = -1;
      end;
    elseif isfield(StudyInfo,'AsegQC') && isfield(StudyInfo,'SurfQC')
      StudyInfo(i).QC_ASEG = 1;
      if isempty(StudyInfo(i).AsegQC)
        StudyInfo(i).QC_ASEG = 0;
        StudyInfo(i).QC = 0;
      elseif ~isempty(regexpi(StudyInfo(i).AsegQC,'bad'))
        StudyInfo(i).QC_ASEG = -1;
        StudyInfo(i).QC = -1;
      end
      StudyInfo(i).QC_surf = 1;
      if isempty(StudyInfo(i).SurfQC)
        StudyInfo(i).QC_surf = 0;
        StudyInfo(i).QC = 0;
      elseif ~isempty(regexpi(StudyInfo(i).SurfQC,'bad'))
        StudyInfo(i).QC_surf = -1;
        StudyInfo(i).QC = -1;
      end
    end;
  end;

  if parms.QC_DTI
    StudyInfo(i).QC_DTI = 1;
    if isempty(StudyInfo(i).DTIQC)
      StudyInfo(i).QC_DTI = 0;
      StudyInfo(i).QC = 0;
    else
      if isnumeric(StudyInfo(i).DTIQC)
        if ~StudyInfo(i).DTIQC
          StudyInfo(i).QC_DTI = -1;
          StudyInfo(i).QC = -1;
        end;
      else
        if isfield(StudyInfo,'DTIQC')
					if ~isempty(regexpi(StudyInfo(i).DTIQC,'bad'))
	          StudyInfo(i).QC_DTI = -1;
					else
	          StudyInfo(i).QC_DTI = 1;
					end;
				else
					StudyInfo(i).QC_DTI = 0;
        end;
				if isfield(StudyInfo,'DTI_regQC')
					if ~isempty(regexpi(StudyInfo(i).DTI_regQC,'bad'))
	 	        StudyInfo(i).QC_DTIreg = -1;
					else
					  StudyInfo(i).QC_DTIreg = 1;
					end;
				else
				  StudyInfo(i).QC_DTIreg = 0;
				end;					
      	if any([StudyInfo(i).QC_DTI,StudyInfo(i).QC_DTIreg]==-1)
	       	StudyInfo(i).QC = -1;
       	end;
      end;
    end;
  end

  if parms.QC_RSI
    StudyInfo(i).QC_RSI = 1;
    if isempty(StudyInfo(i).RSIQC)
      StudyInfo(i).QC_RSI = 0;
      StudyInfo(i).QC = 0;
    else
      if isnumeric(StudyInfo(i).RSIQC)
        if ~StudyInfo(i).RSIQC
          StudyInfo(i).QC_RSI = -1;
          StudyInfo(i).QC = -1;
        end;
      else
        if isfield(StudyInfo,'RSIQC')
          if ~isempty(regexpi(StudyInfo(i).RSIQC,'bad'))
            StudyInfo(i).QC_RSI = -1;
          else
            StudyInfo(i).QC_RSI = 1;
          end;
        else
          StudyInfo(i).QC_RSI = 0;
        end;
        if isfield(StudyInfo,'RSI_regQC')
          if ~isempty(regexpi(StudyInfo(i).RSI_regQC,'bad'))
            StudyInfo(i).QC_RSIreg = -1;
          else
            StudyInfo(i).QC_RSIreg = 1;
          end;
        else
          StudyInfo(i).QC_RSIreg = 0;
        end;          
        if any([StudyInfo(i).QC_RSI,StudyInfo(i).QC_RSIreg]==-1)
          StudyInfo(i).QC = -1;
        end;
      end;
    end;
  end
  
  if parms.QC_BOLD
    StudyInfo(i).QC_BOLD = 1;
    StudyInfo(i).QC_BOLDreg = 1;
    if isempty(StudyInfo(i).BOLDQC)
      StudyInfo(i).QC_BOLD = 0;
      StudyInfo(i).QC_BOLDreg = 0;
      StudyInfo(i).QC = 0;
    else
      if ~isempty(regexpi(StudyInfo(i).BOLDQC,'bad'))
        StudyInfo(i).QC_BOLD = -1;
      end;
      if ~isempty(regexpi(StudyInfo(i).BOLD_regQC,'bad'))
        StudyInfo(i).QC_BOLDreg = -1;
      end;
      if any([StudyInfo(i).QC_BOLD,StudyInfo(i).QC_BOLDreg]==-1)
        StudyInfo(i).QC = -1;
      end;
    end;
  end
end

