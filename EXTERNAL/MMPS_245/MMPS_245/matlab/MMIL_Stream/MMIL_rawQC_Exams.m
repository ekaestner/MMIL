function MMIL_rawQC_Exams(ProjID,varargin)
%function MMIL_rawQC_Exams(ProjID,[options])
%
% Usage:
%  MMIL_rawQC_Exams(ProjID,'key1', value1,...);
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%      orig, raw, proc, proc_dti, proc_bold
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%
% Optional Parameters:
%  'ContainerTypes': cell array of container types to check
%     {default = {'proc' 'proc_dti' 'proc_bold'}}
%  'qccont_flag': [0|1] put output in separate QC containers
%    if RootDirs.qc not specified, will be set to 0
%    {default = 1}
%  'reviewerID': string added to output qcinfo files
%     allows for multiple, independent reviews
%     {default = []}
%  'alt_reviewerIDs': cell array of strings
%     with reviewerIDs of other reviewers
%     {default = []}
%  'max_nrev': maximum number of reviewers per series
%     if alt_reviewerIDs not specified, all reviewers contribute to nrev
%     if alt_reviewerIDs is specified, only those reviewers contribute to nrev
%    {default = 2}
%  'outdir': output directory
%    relative to ContainerPath (must be same as used for autoQC)
%    if empty, will write output to ContainerPath
%    {default = []}
%  'type_sort_flag': [0|1] review exams sorted by container type
%     otherwise, review exams sorted by visit
%    {default = 0}
%  'date_range': vector of minimum and maximum date to review
%     e.g. [20161101,20161201]
%     if not supplied, review all
%    {default = []}
%  'precheck_flag': [0|1] perform check at start
%      to  exclude previously QC'd visits for each type
%    {default = 1}
%  'fname_qcinfo': file name of qcinfo summary csv
%      to exclude previously QC'd, skipping precheck
%    {default = []}
%  'revdisp_flag': [0|1|2] review series with disparity between reviewers
%     only applicable if fname_qcinfo is supplied
%     with forceflag = 1 if reviewing with disparity
%     0: do not review because of disparity
%     1: review series with or without disparity
%     2: review series only with disparity
%    {default = 0}
%  'verbose': [0|1] display status messages and warnings
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  03/20/16 by Don Hagler
% Prev Mod: 05/12/17 by Don Hagler
% Last Mod: 08/25/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);

% check whether visits have already been QC'd and exclude them
if ~isempty(parms.fname_qcinfo)
  parms = load_qcinfo(parms);
  parms.precheck_flag = 1;
elseif parms.precheck_flag
  parms = precheck_visits(parms);
else
  parms.visit_lists = cell(parms.ntypes,1);
  for k=1:parms.ntypes
    parms.visit_lists{k} = [1:parms.nvisits];
  end;
end;

% exclude visits outside date range
if ~isempty(parms.date_range)
  parms = precheck_dates(parms);
end;

if parms.type_sort_flag
  % loop over container types, then visits
  for k=1:parms.ntypes
    vlist = parms.visit_lists{k};
    nvisits = length(vlist);
    if parms.verbose
      fprintf('%s: reviewing %d visits for type %s...\n',...
        mfilename,nvisits,parms.ContainerTypes{k});
    end;
    for i=1:nvisits
      s = vlist(i);
      review_exam(s,k,parms);
    end;
  end
else
  % loop over visits, then container types
  vlist = parms.visit_lists{1};
  for k=2:parms.ntypes
    [tmp_vlist,ind_vlist] = setdiff(parms.visit_lists{k},vlist);
    if ~isempty(tmp_vlist)
      [tmp,ind_sort] = sort(ind_vlist);
      tmp_vlist = tmp_vlist(ind_sort);
      vlist = cat(1,vlist,tmp_vlist);
    end;
  end;
  nvisits = length(vlist);
  if parms.verbose
    fprintf('%s: reviewing all types for %d visits...\n',...
      mfilename,nvisits);
  end;
  for i=1:nvisits
    s = vlist(i);
    for k=1:parms.ntypes
      review_exam(s,k,parms);
    end;
  end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms = mmil_args2parms(options,{...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
  ...
    'ContainerTypes',{'proc' 'proc_dti' 'proc_bold'},{'proc' 'proc_dti' 'proc_bold'},...
    'qccont_flag',true,[false true],...
    'reviewerID',[],[],...
    'alt_reviewerIDs',[],[],...
    'max_nrev',2,[1,100],...
    'outdir',[],[],...
    'type_sort_flag',false,[false true],...
    'date_range',[],[],...
    'precheck_flag',true,[false true],...
    'fname_qcinfo',[],[],...
    'revdisp_flag',0,[0:2],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % hidden
    'raw_pat','MRIRAW_(?<VisitID>[^]+)_(?<StudyDate>\d{8})\..+',[],...
  ...
    'rawQC_tags',{'reviewerID','alt_reviewerIDs','max_nrev',...
                  'outdir','verbose','forceflag'},[],...
    'info_tags',{'RootDirs','ignore_VisitInfo_flag','user','numvec_tags'},[],...
  });
  parms.ntypes = length(parms.ContainerTypes);

  % check date_range
  if ~isempty(parms.date_range)
    if ~isnumeric(parms.date_range)
      error('date_range must be numeric');
    end;
    if length(parms.date_range)~=2
      error('date_range must have two elements');
    end;
    parms.date_range = sort(parms.date_range);
    parms.datenum_min  = datenum(num2str(parms.date_range(1)),'yyyymmdd');
    parms.datenum_max  = datenum(num2str(parms.date_range(2)),'yyyymmdd');
  end;

  % get StudyInfo and RootDirs
  if isempty(parms.StudyInfo)
    if parms.verbose
      fprintf('%s: getting StudyInfo for %s...\n',mfilename,ProjID);
    end;
    args = mmil_parms2args(parms,parms.info_tags);
    [parms.StudyInfo,parms.RootDirs] = MMIL_Quick_StudyInfo(ProjID,args{:});
  elseif parms.verbose
    fprintf('%s: getting list of VisitIDs from user supplied StudyInfo...\n',...
      mfilename);
  end;
  if isempty(parms.RootDirs)
    [ProjInfo,parms.RootDirs] = MMIL_Get_ProjInfo(ProjID);
  end;
  
  % check for qc field in RootDirs
  if parms.qccont_flag && isempty(mmil_getfield(parms.RootDirs,'qc'))
    parms.qccont_flag = 0;
    if parms.verbose
      fprintf('%s: WARNING: no qc RootDir specified, so setting qccont_flag = 0\n',...
        mfilename);
    end;
  end;  

  % get list of VisitIDs
  parms.VisitIDs = {parms.StudyInfo.VisitID};
  parms.nvisits = length(parms.VisitIDs);
  if parms.verbose
    fprintf('%s: %d visits and %d types to review\n',...
      mfilename,parms.nvisits,parms.ntypes);
  end;

  % check qcinfo file
  if ~isempty(parms.fname_qcinfo)
    if ~exist(parms.fname_qcinfo,'file')
      error('file %s not found',parms.fname_qcinfo);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_qcinfo(parms);
  if parms.verbose
    fprintf('%s: loading qcinfo file to exclude previously QC''d visits...\n',...
      mfilename);
    tic;
  end;
  parms.visit_lists = cell(parms.ntypes,1);
  qcinfo = mmil_csv2struct(parms.fname_qcinfo);
  if parms.verbose, toc; end;
  tag = 'QC';
  if ~isempty(parms.reviewerID)
    tag = [tag '_' parms.reviewerID];
  end;
  if ~isfield(qcinfo,tag)
    fprintf('%s: WARNING: QC column %s not found in info file %s\n',...
      mfilename,tag,parms.fname_qcinfo);
    ind_excl = [];
  else
    QC = {qcinfo.(tag)};
    ind_excl = find(~cellfun(@isempty,QC));
  end;
  % exclude visits with sufficient review
  if parms.max_nrev>0
    if ~isempty(parms.alt_reviewerIDs)
      % check each altID
      nalt = length(parms.alt_reviewerIDs);
      nrev = zeros(length(qcinfo),1);
      for i=1:nalt
        altID = parms.alt_reviewerIDs{i};
        tag = ['QC_' altID];
        if ~isfield(qcinfo,tag)
          fprintf('%s: WARNING: QC column %s not found in info file %s\n',...
            mfilename,tag,parms.fname_qcinfo);
        else
          QC_alt = {qcinfo.(tag)};
          ind_alt = find(~cellfun(@isempty,QC_alt));
          nrev(ind_alt) = nrev(ind_alt) + 1;
        end;
      end;
      ind_excl = union(ind_excl,find(nrev >= parms.max_nrev));
    end;
  end;
  parms.qcinfo_VisitIDs = {qcinfo.VisitID};
  % include visits with disparity
  if parms.revdisp_flag
    parms.revdisp = [qcinfo.revdisp];
    ind_revdisp = find(parms.revdisp);
    % review visits with or without disparity
    ind_excl = setdiff(ind_excl,ind_revdisp);
  end;
  if parms.revdisp_flag==2
    if parms.verbose
      fprintf('%s: excluding all but %d previously QC''d visits with disparity...\n',...
        mfilename,length(ind_revdisp));
    end;
    visitlist = parms.qcinfo_VisitIDs(ind_revdisp);
    [tmp,ind] = intersect(parms.VisitIDs,visitlist);
  else
    % exclude previously QC'd visits
    visitlist = parms.qcinfo_VisitIDs(ind_excl);
    if parms.verbose
      fprintf('%s: excluding %d previously QC''d visits...\n',...
        mfilename,length(visitlist));
      if parms.revdisp_flag==1
        fprintf('%s: including %d previously QC''d visits with disparity...\n',...
          mfilename,length(ind_revdisp));
      end;
    end;
    [tmp,ind] = setdiff(parms.VisitIDs,visitlist);
  end;
  for k=1:parms.ntypes
    parms.visit_lists{k} = ind;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = precheck_visits(parms)
  if parms.verbose
    fprintf('%s: performing precheck to exclude previously QC''d visits...\n',...
      mfilename);
    tic;
  end;
  parms.visit_lists = cell(parms.ntypes,1);
  for k=1:parms.ntypes
    if isempty(parms.reviewerID)
      suffix = '_raw_qcinfo';
    else
      suffix = sprintf('_raw_%s_qcinfo',parms.reviewerID);
    end;
    if parms.qccont_flag
      cmd = sprintf('find %s -name "*%s.mat"',...
        parms.RootDirs.qc,suffix);
    else
      cmd = sprintf('find %s -name "*%s.mat"',...
        parms.RootDirs.(parms.ContainerTypes{k}),suffix);
    end;
    [s,r] = unix(cmd);
    if s
      error('cmd %s failed:\n%s',cmd,r);
    end;
    if ~isempty(r)
      % separate unix output into lines
      c=textscan(r,'%s','delimiter','\n');
      flist = c{1};
      % get directory name from full path name of mat file
      dirlist = cellfun(@fileparts, flist, 'UniformOutput',false);
      % remove path
      proclist = cellfun(@(x) regexprep(x,'\w+/',''),dirlist, 'UniformOutput',false);
      % remove /
      proclist = cellfun(@(x) regexprep(x,'/',''),proclist, 'UniformOutput',false);
      % remove PROC dir stem
      if parms.qccont_flag
        visitlist = cellfun(@(x) regexprep(x,'\w+QC_',''),proclist, 'UniformOutput',false);
      else
        visitlist = cellfun(@(x) regexprep(x,'\w+PROC_',''),proclist, 'UniformOutput',false);
      end;
      % remove date and time stamp
      visitlist = cellfun(@(x) regexprep(x,'_\d{8}\.\d{6}_1$',''),visitlist, 'UniformOutput',false);
      if parms.verbose
        fprintf('%s: excluding %d previously QC''d visits...\n',...
          mfilename,length(visitlist));
      end;
      % remove previously QC'd visits from VisitIDs
      [tmp,ind] = setdiff(parms.VisitIDs,visitlist);
      parms.visit_lists{k} = ind;
    else
      parms.visit_lists{k} = [1:parms.nvisits];
    end;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = precheck_dates(parms)
  if parms.verbose
    fprintf('%s: performing precheck to exclude outside of date range...\n',...
      mfilename);
    tic;
  end;
  vlist = unique(cat(1,parms.visit_lists{:}));
  VisitIDs = parms.VisitIDs(vlist);
  k = regexp(VisitIDs,'_(?<StudyDate>\d{8})_','names');
  StudyDates = cellfun(@(x) mmil_getfield(x,'StudyDate',00000000),k,'UniformOutput',false);
  ind_nonempty = find(~cellfun(@isempty,StudyDates));
  if ~isempty(ind_nonempty)
    datenums_study = datenum(StudyDates(ind_nonempty),'yyyymmdd');
    ind_excl = find(datenums_study < parms.datenum_min |...
                    datenums_study > parms.datenum_max);
    visitlist = VisitIDs(ind_nonempty(ind_excl));
    if parms.verbose
      fprintf('%s: excluding %d visits outside date range...\n',...
        mfilename,length(visitlist));
    end;
    [tmp,ind] = setdiff(VisitIDs,visitlist);
    if isempty(ind)
      for k=1:parms.ntypes
        parms.visit_lists{k} = [];
      end;
    else
      tmp_StudyDates = StudyDates(ind);
      [tmp,ind_sort] = sort(tmp_StudyDates);
      ind_keep = vlist(ind(ind_sort));
      for k=1:parms.ntypes
        [ind_visits,ind_v,ind_k] = intersect(parms.visit_lists{k},ind_keep);
        [tmp,ind_sort] = sort(ind_k);
        parms.visit_lists{k} = ind_visits(ind_sort);
      end;
    end;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: no longer used
function valid_date_flag = check_date(s,parms)
  valid_date_flag = 0;
  VisitID = parms.StudyInfo(s).VisitID;
  ContainerType = 'raw';
  [ContainerPath,ContainerDir] = ...
    MMIL_Get_Container(parms.RootDirs,VisitID,ContainerType);
  if ~isempty(ContainerDir)
    n = regexp(ContainerDir,parms.raw_pat,'names');
    if isempty(n)
      fprintf('%s: WARNING: Container %s does not have expected pattern\n',...
        mfilename,ContainerDir);
      return;
    else
      datenum_study = datenum(n.StudyDate,'yyyymmdd');
      if datenum_study >= parms.datenum_min &&...
         datenum_study <= parms.datenum_max
        valid_date_flag = 1;
      else
        fprintf('%s: skipping %s because outside date_range\n',...
          mfilename,VisitID);
      end;      
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function review_exam(s,k,parms)
  VisitID = parms.VisitIDs{s};
  ContainerType = parms.ContainerTypes{k};
  if parms.verbose
    fprintf('%s: preparing to review exam for %s %s...\n',...
    	mfilename,VisitID,ContainerType);
  end;
  [ContainerPath,ContainerDir] =...
    MMIL_Get_Container(parms.RootDirs,VisitID,ContainerType);
  if isempty(ContainerPath)
    if parms.verbose
      fprintf('%s: container not found for %s %s... skipping\n',...
      	mfilename,VisitID,ContainerType);
    end;
    return;
  end;
  if parms.revdisp_flag
    q = find(strcmp(VisitID,parms.qcinfo_VisitIDs));
    if ~isempty(q)
      if parms.revdisp(q)
        parms.forceflag = 1;
        if parms.verbose
          fprintf('%s: setting forceflag=1 because of reviewer disparity for %s...\n',...
            mfilename,VisitID);
        end;
      end;
    end;
  end;
  if parms.qccont_flag
    parms.outdir = sprintf('%s/%s',parms.RootDirs.qc,...
                                   regexprep(ContainerDir,'PROC','QC'));
  end;
  if ~isempty(ContainerPath)
    if parms.verbose
      fprintf('%s: reviewing exam for %s %s...\n',...
      	mfilename,VisitID,ContainerType);
    end;
    args = mmil_parms2args(parms,parms.rawQC_tags);
    MMIL_rawQC_Exam(ContainerPath,args{:});
  end;
return;

