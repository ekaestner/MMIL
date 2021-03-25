function [StudyInfo,RootDirs] = MMIL_Check_StudyInfo(StudyInfo,RootDirs,varargin)
%function [StudyInfo,RootDirs] = MMIL_Check_StudyInfo(StudyInfo,RootDirs,[options])
%
% Usage:
%  [StudyInfo,RootDirs] = MMIL_Check_StudyInfo(StudyInfo,RootDirs,'key1', value1,...);
%
% Required Parameters:
%  StudyInfo: struct array of study information
%   (e.g. read from csv file with MMIL_Read_StudyInfo)
%   must contain one of these fields:
%     VisitID: name of directory in orig data dir
%     MEG_VisitID: name of directory in orig_meg data dir
%     PET_VisitID: name of directory in orig_pet data dir
%   may contain these MRI fields (or will be determined from VisitID)
%     raw: name of raw data container
%     proc: name of processed data container
%     proc_dti: name of processed DTI data container
%     proc_bold: name of processed BOLD data container
%     fsurf: name of freesurfer recon container
%   may container these MEG fields (or will be determined from MEG_VisitID)
%     raw_meg: name of raw MEG data container
%     proc_meg: name of processed MEG data container
%   may contain these PET fields (or will be determined from PET_VisitID)
%     raw_pet: name of raw (unpacked) PET data container
%     proc_pet: name of processed PET data container
%   if empty, use all subjects found in RootDirs.raw and other
%  RootDirs: struct that may contain the following fields:
%   home, batch
%   orig, raw, proc, proc_dti, proc_bold
%   fsurf, fsico, long
%   orig_meg, raw_meg, proc_meg
%   orig_pet, raw_pet, proc_pet
%     these specify the full paths of root directories containing data containers
%
% Optional Parameters
%  'VisitIDs': cell array of VisitIDs to be included in StudyInfo
%    if not supplied, include all VisitIDs
%    {default = []}
%  'SubjIDs': cell array of SubjIDs to be included in StudyInfo
%    if not supplied, include all SubjIDs
%    {default = []}
%  'checkflag': [0|1] check each visit for required containers, etc.
%    and get session ID; otherwise, only get VisitID from orig dirs
%    {default = 1}
%  'qcflag': [0|1] only include subjects with StudyInfo.QC = 1
%    {default = 0}
%  'modality': limit returned StudyInfo to sessions of this modality only
%    allowed values: 'MRI', 'MEG', 'PET'
%    {default = []}
%  'required_containers': cell array of container names required to exist
%    e.g. 'raw', 'proc', 'proc_dti', 'proc_bold',
%         'fsurf', 'fsico', 'meg_raw', 'meg_proc', etc.
%    if speficied, all entries will be automatically added to 'required_rootdirs'
%    {default = []}
%  'required_rootdirs': cell array of root directories required to exist
%    e.g. 'raw', 'proc', 'proc_dti', 'proc_bold',
%         'fsurf', 'fsico', 'meg_raw', 'meg_proc', etc.
%    {default = []}
%  'ico': icosahedral order number (e.g 1-7)
%    used to specify names of fsico containers
%    {default = 7}
%
% Output:
%   StudyInfo: struct array of study information
%     created from Container names (if not supplied as input)
%     or filtered by checking for valid scans and rejected subjects (QC=0)
%   RootDirs: struct containing empty values for fields not found in input RootDirs
%
% Created:  03/31/09 by Don Hagler
% Last Mod: 04/03/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms_filter = {...
  'VisitIDs',[],[],...
  'SubjIDs',[],[],...
  'checkflag',true,[false true],...
  'qcflag',false,[false true],...
  'modality',[],{'MRI','MEG','PET'},...
  'required_containers',[],[],...
  'required_rootdirs',[],[],...
  'ico',7,[1:7],...
  'numvec_tags',[],[],...
...
  'default_numvec_tags',{'DCM_RejectSeries'...
       'DTIScanNums' 'DTIScanNums2' 'DTIScanNums3'...
       'BOLDScanNums' 'DTI_resolution' 'DTI_nvoxels' 'BOLD_resolution'...
       'BOLD_nvoxels' 'FP_snums' 'FP_revflags' 'FE_snums' 'FE_revflags'...
       'GLM_snums' 'PROC_valid_event_codes' 'RF_pol_snums' 'RF_ecc_snums' ...
       'RCSE_r_offset_range' 'RCSE_th_offset_range'},[],...    
  'rawContainerStem','MRIRAW',[],...
  'raw_pat','MRIRAW_(?<VisitID>[^]+)_(?<StudyDate>\d{8})\..+',[],...
  'ContainerTypes',{'orig' 'raw' 'proc' 'proc_dti' 'proc_bold'...
    'fsurf' 'long' 'fsico'...
    'raw_asl' 'proc_asl'...
    'orig_meg' 'raw_meg' 'proc_meg'...
    'orig_pet','raw_pet','proc_pet'},[],...
  'standard_tags',{'SubjID' 'VisitID' 'MEG_VisitID' 'STRUCT_VisitID'...
      'GLM_ROI_VisitID' 'PET_VisitID' 'StudyDate' 'SessID' 'VisitNumber'...      
      'orig' 'raw' 'proc' 'proc_dti' 'proc_bold' 'fsurf' 'long'...
      'fsico' 'fsico1' 'fsico2' 'fsico3' 'fsico4' 'fsico5' 'fsico6'...
      'raw_asl' 'proc_asl'...
      'orig_meg' 'raw_meg' 'proc_meg'...
      'valid_event_codes' 'orig_pet' 'raw_pet' 'proc_pet' 'QC' 'Modality'...
      'DCM_RejectSeries'...
      'DTIScanNums' 'DTIScanNums2' 'DTIScanNums3' 'DTI_resolution'...
      'DTI_nvoxels' 'BOLDScanNums' 'BOLD_resolution' 'BOLD_nvoxels'...
      'FP_snums' 'FP_revflags' 'FE_snums' 'FE_revflags' 'GLM_snums'},[],...
  'str_tags',{'SubjID' 'VisitID' 'MEG_VisitID' 'STRUCT_VisitID' ...
              'GLM_ROI_VisitID' 'PET_VisitID'},[],...
};


% extra 'standard_tags':
%      'RawQC' 'AsegQC' 'AsegQCNotes' 'SurfQC' 'SurfQCNotes' 'PETRegQC'...
%      'PETRegQCNotes'},[],...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin,parms_filter);

% make sure RootDirs has required fields (and add empty values for others)
RootDirs = MMIL_Check_RootDirs(RootDirs,parms.required_rootdirs);

% restrict to specified VisitIDs
StudyInfo = restrict_info(StudyInfo,parms.VisitIDs,'VisitID');

% restrict to specified SubjIDs
StudyInfo = restrict_info(StudyInfo,parms.SubjIDs,'SubjID');

% if no StudyInfo supplied, generate from orig rootdirs; initialize fields
StudyInfo = init_StudyInfo(StudyInfo,RootDirs,parms);

if parms.checkflag
  for s=1:length(StudyInfo)
    % check values for some fields, set others
    StudyInfo(s) = check_fields(StudyInfo(s),parms);

    % find each type of container for given VisitIDs
    StudyInfo(s) = get_containers(StudyInfo(s),RootDirs,parms);

    % get VisitID, StudyDate, and SessID from container name
    StudyInfo(s) = get_SessID(StudyInfo(s),RootDirs,parms);

    % check for missing containers
    StudyInfo(s) = check_containers(StudyInfo(s),parms,i);

    % calculate Age from DOB and StudyDate
    if isfield(StudyInfo(s),'DOB')
      if ~isempty(StudyInfo(s).StudyDate) & ~isempty(StudyInfo(s).DOB)
        StudyInfo(s).Age = ...
          (datenum(num2str(StudyInfo(s).StudyDate),'yyyymmdd') -...
           datenum(num2str(StudyInfo(s).DOB),'yyyymmdd'))/365.25;
      else
        StudyInfo(s).DOB = [];
        StudyInfo(s).Age = [];
      end;
    end; 
  end;
end;

if parms.qcflag % only include subjects with QC=1
  i_QC = find([StudyInfo.QC]==1);
elseif parms.checkflag % only include subjects with all required containers
  i_QC = find([StudyInfo.QC]~=-100);
else
  i_QC = [1:length(StudyInfo)];
end;
StudyInfo = StudyInfo(i_QC);

% include only studies with correct modality
if ~isempty(parms.modality)
  StudyInfo = StudyInfo(strcmp(upper({StudyInfo.Modality}),parms.modality));
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(varargin,parms_filter)
  parms = mmil_args2parms(varargin,parms_filter);

  if isempty(parms.numvec_tags)
    parms.numvec_tags = parms.default_numvec_tags;
  else
    if ~iscell(parms.numvec_tags)
      parms.numvec_tags = {parms.numvec_tags};
    end;
    parms.numvec_tags = reshape(parms.numvec_tags,[1,numel(parms.numvec_tags)]);
    parms.numvec_tags = cat(2,parms.numvec_tags,parms.default_numvec_tags);
  end;

  if parms.ico==7
    parms.fsico_type = 'fsico';
  else
    parms.fsico_type = sprintf('fsico%d',parms.ico);
  end;

  % add required_containers to required_rootdirs
  if ~isempty(parms.required_containers)
    parms.required_rootdirs = union(parms.required_rootdirs,parms.required_containers);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [StudyInfo,RootDirs] = init_StudyInfo(StudyInfo,RootDirs,parms)
  if isempty(StudyInfo)
    j = 1;
    dir_tags = {'orig','orig_meg','orig_pet'};
    visit_tags = {'VisitID','MEG_VisitID','PET_VisitID'};
    for d=1:length(dir_tags)
      dir_tag = dir_tags{d};
      visit_tag = visit_tags{d};
      rootdir = mmil_getfield(RootDirs,dir_tag,[]);
      if isempty(rootdir), continue; end;
      dirlist = dir(rootdir);
      % get rid of . and .. in dirlist
      n = cellfun('isempty',regexp({dirlist.name},'^\.'));
      dirlist = dirlist(n);
      dirlist = {dirlist.name};
      for i=1:length(dirlist)
        StudyInfo(j).(dir_tag) = dirlist{i};
        % replace any '.' in VisitID with '_'
        %  (causes problems when inserted in job names and for ADNI)
        StudyInfo(j).(visit_tag) = regexprep(dirlist{i},'\.','_');
        j = j + 1;
      end
    end;
  end;
  % initialize all standard fields
  for i=1:length(StudyInfo)
    for f=1:length(parms.standard_tags)
      StudyInfo(i).(parms.standard_tags{f}) = ...
        mmil_getfield(StudyInfo(i),parms.standard_tags{f},[]);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StudyInfo = restrict_info(StudyInfo,IDs,IDname)
  if ~isempty(IDs) &&...
     (isempty(StudyInfo) || isfield(StudyInfo,IDname))
    if isempty(StudyInfo)
      ind_add = [1:length(IDs)];
    else
      all_IDs = {StudyInfo.(IDname)};
      ind_keep = find(ismember(all_IDs,IDs));
      [tmp,ind_add] = setdiff(IDs,all_IDs);
      if ~isempty(ind_keep)
        % keep StudyInfo for matching IDs
        StudyInfo = StudyInfo(ind_keep);
      else
        StudyInfo = [];
      end;
      if ~isempty(ind_add)
        fprintf('%s: WARNING: %s not found in VisitInfo: %s\n',...
          mfilename,IDname,sprintf('%s ',IDs{ind_add}));
      end;  
    end;
    % add missing IDs to StudyInfo
    for i=1:length(ind_add)
      StudyInfo(end+1).(IDname) = IDs{ind_add(i)};
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = check_fields(info,parms)
  % assume we include study if not specified
  if isempty(info.QC), info.QC = 1; end;
  % set modality unless already specified, default to 'MRI'
  if isempty(info.Modality)
    if ~isempty(info.MEG_VisitID)
      info.Modality = 'MEG';
    elseif ~isempty(info.PET_VisitID)
      info.Modality = 'PET';
    else
      info.Modality = 'MRI';
    end;
  end;
  % make sure ID's are strings, not numbers
  for s=1:length(parms.str_tags)
    tag = parms.str_tags{s};
    if ~isempty(info.(tag)) && isnumeric(info.(tag))
      info.(tag) = num2str(info.(tag));
    end;
  end;
  % set VisitID or SubjID appropriately if not set
  switch info.Modality
    case 'MRI'
      if isempty(info.SubjID), info.SubjID = info.VisitID; end;
      if isempty(info.VisitID), info.VisitID = info.SubjID; end;
    case 'MEG'
      if isempty(info.SubjID), info.SubjID = info.MEG_VisitID; end;
      if isempty(info.MEG_VisitID), info.MEG_VisitID = info.SubjID; end;
    case 'PET'
      if isempty(info.SubjID), info.SubjID = info.PET_VisitID; end;
      if isempty(info.PET_VisitID), info.PET_VisitID = info.SubjID; end;
      % unless otherwise specfied assume that PET and matching MRI
      %   visits have same VisitID
      if isempty(info.VisitID), info.VisitID = info.PET_VisitID; end;
  end;
  if isempty(info.STRUCT_VisitID), info.STRUCT_VisitID = info.VisitID; end;
  if isempty(info.VisitNumber), info.VisitNumber = 1; end;
  % convert certain fields to from strings to numeric vectors
  info = mmil_structarr_str2num(info,parms.numvec_tags);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = get_containers(info,RootDirs,parms)
  % decide whether we need to get containers
  getflag = 0;
  for t=1:length(parms.ContainerTypes)
    ContainerType = parms.ContainerTypes{t};
    if isempty(info.(ContainerType)), getflag = 1; end;
  end;
  if ~getflag, return; end;
  for t=1:length(parms.ContainerTypes)
    ContainerType = parms.ContainerTypes{t};
    DirType = ContainerType;
    if strcmp(ContainerType,'fsico'), ContainerType = parms.fsico_type; end;
    if isempty(info.(ContainerType)) && ~isempty(RootDirs.(DirType))
      if ismember(DirType,{'fsurf','fsico'})
        VisitID = info.STRUCT_VisitID;
        if isempty(VisitID), continue; end;
      elseif ~isempty(regexp(ContainerType,'long'))
        VisitID = info.SubjID;
        if isempty(VisitID), continue; end;
      elseif ~isempty(regexp(ContainerType,'meg'))
        VisitID = info.MEG_VisitID;
        if isempty(VisitID), continue; end;
      elseif ~isempty(regexp(ContainerType,'pet'))
        VisitID = info.PET_VisitID;
        if isempty(VisitID), continue; end;
      else
        VisitID = info.VisitID;
      end;
      if ~isempty(VisitID)
        [tmp,info.(ContainerType)] = ...
          MMIL_Get_Container(RootDirs,VisitID,ContainerType);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = get_SessID(info,RootDirs,parms)
  if ~isempty(info.raw) && exist([RootDirs.raw '/' info.raw]) &&...
     (isempty(info.VisitID) || isempty(info.StudyDate))
    n = regexp(char(info.raw),parms.raw_pat,'names');
    if isempty(n)
      fprintf('%s: WARNING: Container %s does not have expected pattern\n',...
        mfilename,info.raw);
      if isempty(info.VisitID), info.VisitID = info.raw; end;
    else
      if isempty(info.VisitID), info.VisitID = n.VisitID; end;
      if isempty(info.StudyDate), info.StudyDate = str2num(n.StudyDate); end;
    end;
    if isempty(info.SubjID)
      info.SubjID = info.VisitID;
    end;
  end;
  if isempty(info.SessID)
    if ~isempty(info.raw)
      n = regexp(info.raw,[parms.rawContainerStem '_(?<sessid>.+)'],'names');
      if ~isempty(n), info.SessID = n.sessid; end;    
    else
      info.SessID = info.VisitID;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = check_containers(info,parms,studynum)
  % check for missing containers
  for j=1:length(parms.required_containers)
    ContainerType = parms.required_containers{j};
    if strcmp(ContainerType,'fsico'), ContainerType = parms.fsico_type; end;
    if isempty(info.(ContainerType))
      switch info.Modality
        case 'MRI'
          VisitID = info.VisitID;
          label = 'VisitID';
        case 'MEG'
          VisitID = info.MEG_VisitID;
          label = 'MEG_VisitID';
        case 'PET'
          VisitID = info.PET_VisitID;
          label = 'PET_VisitID';
      end;
      if isempty(VisitID)
        VisitID = num2str(studynum);
        label = 'study';
      end;
      fprintf('%s: WARNING: missing %s Container for %s %s\n',...
        mfilename,ContainerType,label,VisitID);
      info.QC = -100;
      break;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

