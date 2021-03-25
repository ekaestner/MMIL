function VisitInfo = MMIL_Update_VisitInfo(ProjID,varargin)
%function VisitInfo = MMIL_Update_VisitInfo(ProjID,[options])
%
% Purpose: update VisitInfo file for a project, get SubjID from VisitID
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and VisitInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%
% Optional Input:
%  'add_missing_flag': [0|1] add to VisitInfo
%     visits in RootDirs.orig not found in existing csv files
%     {default = 1}
%  'SubjID_pattern': regexp pattern used for getting SubjID
%     {default = '(?<SubjID>^[^_]+)_\w+'}
%
% Created:  11/21/13 by Don Hagler
% Last Mod: 02/15/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VisitInfo = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
fprintf('%s: checking input...\n',mfilename);
parms = check_input(ProjID,varargin);

% initialize VisitInfo struct
fprintf('%s: initializing VisitInfo...\n',mfilename);
VisitInfo = init_VisitInfo(parms);

% get SubjID from VisitID
fprintf('%s: getting SubjIDs from VisitIDs...\n',mfilename);
VisitInfo = get_SubjIDs(parms,VisitInfo);

% sort rows by VisitID
fprintf('%s: sorting rows by VisitID...\n',mfilename);
VisitInfo = sort_VisitInfo(parms,VisitInfo);

% write VisitInfo csv file
%save_flag = ...
%  input('save updated VisitInfo? [y/n] ','s');
%if ~strcmp(save_flag,'y'), return; end;
save_VisitInfo(parms,VisitInfo);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
  ... % optional
    'add_missing_flag',true,[false true],...
    'SubjID_pattern','(?<SubjID>^[^_]+)_\w+',[],...
  ... % hidden
    'RootDirs',[],[],...
    'StudyInfo',[],[],...
    'VisitIDs',[],[],...
    'SubjIDs',[],[],...
    'user',[],[],...
  ...
    'info_tags',{'VisitIDs','SubjIDs','StudyInfo','RootDirs','ignore_VisitInfo_flag','user','numvec_tags'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  % set ignore_csv_flag
  parms.ignore_VisitInfo_flag = parms.add_missing_flag;
  % use ProjID to get StudyInfo and RootDirs
  args = mmil_parms2args(parms,parms.info_tags);
  try
    [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Quick_Check_ProjID(ProjID,args{:});
  catch me
    fprintf('%s: WARNING: ProjInfo check failed: %s\n',mfilename,me.message);
    ProjInfo = [];
  end;
  % check StudyInfo, get VisitIDs
  if length(parms.StudyInfo)==0
    error('no valid subjects');
  end;
  parms.VisitIDs = {parms.StudyInfo.VisitID};
  % set VisitInfo file name
  parms.fname_VisitInfo = sprintf('%s/ProjInfo/%s/%s_VisitInfo.csv',...
    parms.RootDirs.home,ProjID,ProjID);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitInfo = init_VisitInfo(parms)
  % create VisitInfo csv file if does not exist
  if ~exist(parms.fname_VisitInfo,'file')
    headers = {'SubjID','DCM_RejectSeries'};
    data = cell(length(parms.VisitIDs),length(headers));
    mmil_write_csv(parms.fname_VisitInfo,data,'firstcol_label','VisitID',...
      'row_labels', parms.VisitIDs, 'col_labels', headers);
  end;
  VisitInfo = mmil_csv2struct(parms.fname_VisitInfo);
  % check for numeric VisitIDs
  for v=1:length(VisitInfo)
    if isnumeric(VisitInfo(v).VisitID)
      VisitInfo(v).VisitID = num2str(VisitInfo(v).VisitID);
    end;
  end;
  % add new subjects to VisitInfo
  if parms.add_missing_flag
    VisitIDs = ... new visit IDs to add
      parms.VisitIDs(~ismember(parms.VisitIDs,{VisitInfo.VisitID}));
    for v = 1:length(VisitIDs)
      VisitInfo(end+1).VisitID = VisitIDs{v};
    end
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitInfo = get_SubjIDs(parms,VisitInfo)
  nvisits = length(VisitInfo);
  for v=1:nvisits
    if ~isfield(VisitInfo,'SubjID'), VisitInfo(v).SubjID = []; end;
    if isempty(VisitInfo(v).SubjID)
      % get SubjID from VisitID using regexp
      n = regexp(VisitInfo(v).VisitID,parms.SubjID_pattern,'names');
      if ~isempty(n)
        VisitInfo(v).SubjID = n.SubjID;
        if isfield(n,'VisitNumber')
          VisitInfo(v).VisitNumber = str2num(n.VisitNumber);
        end;
      else
        fprintf('%s: WARNING: VisitID %s does not match SubjID_pattern\n',...
          mfilename,VisitInfo(v).VisitID);
        VisitInfo(v).SubjID = VisitInfo(v).VisitID;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitInfo = sort_VisitInfo(parms,VisitInfo)
  VisitIDs = {VisitInfo.VisitID};
  [tmp,ind_sort] = sort(VisitIDs);
  VisitInfo = VisitInfo(ind_sort);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_VisitInfo(parms,VisitInfo)
  tmpinfo = squeeze(struct2cell(VisitInfo))';
  tmpnames = fieldnames(VisitInfo)';
  ind_SubjID = find(strcmp(tmpnames,'SubjID'));
  ind_VisitID = find(strcmp(tmpnames,'VisitID'));
  ind_other = setdiff([1:length(tmpnames)],[ind_SubjID,ind_VisitID]);
  new_order = [ind_SubjID,ind_VisitID,ind_other];
  tmpinfo = tmpinfo(:,new_order);
  tmpnames = tmpnames(new_order);
  mmil_write_csv(parms.fname_VisitInfo,tmpinfo,'col_labels',tmpnames);
return;

