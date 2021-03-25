function qcinfo = MMIL_Summarize_postQC(ProjID,varargin)
% function qcinfo = MMIL_Summarize_postQC(ProjID,[options])
%
% Usage:
%  qcinfo = MMIL_Summarize_postQC(ProjID,'key1', value1,...);
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
%      fsurf, fsico
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%
% Optional Parameters:
%  'ContainerTypes': cell array of container types to convert
%     {default = {'fsurf' 'proc_dti' 'proc_bold'}}
%  'qccont_flag': [0|1] put output in separate QC containers
%    if RootDirs.qc not specified, will be set to 0
%    {default = 1}
%  'reviewerIDs': cell array of reviewer ID strings included in qcinfo files
%     allows for multiple, independent reviews
%     {default = []}
%  'qc_outfix': string  added to output qcinfo files
%     allows for reviews with different criteria
%     {Default = 'post'}
%  'fname_out': output file name
%    if not supplied, will be use outdir, ProjID, and outfix
%    {default = []}
%  'outdir': output directory
%    if not supplied, will be /home/{user}/MetaData/{ProjID}
%    {default = []}
%  'outfix': suffix attached to output filename
%    {default = 'post_qcinfo'}
%  'type_sort_flag': [0|1] review exams sorted by container type
%     otherwise, review exams sorted by visit
%    {default = 0}
%  'verbose': [0|1] display status messages and warnings
%    {default: 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 1}
%
% Output:
%   qcinfo: struct array containing postQC information
%     for each visit
%
% Created:  10/27/17 by Don Hagler
% Prev Mod: 10/27/17 by Don Hagler
% Last Mod: 10/31/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qcinfo = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);

if ~exist(parms.fname_out,'file') || parms.forceflag
  if parms.type_sort_flag
    % loop over container types, then visits
    for k=1:parms.ntypes
      for s=1:parms.nvisits
        tmpinfo = get_qcinfo(s,k,parms);
        qcinfo = concat_qcinfo(qcinfo,tmpinfo);
      end;
    end
  else
    % loop over visits, then container types
    for s=1:parms.nvisits
      for k=1:parms.ntypes
        tmpinfo = get_qcinfo(s,k,parms);
        qcinfo = concat_qcinfo(qcinfo,tmpinfo);
      end;
    end
  end;
  if isempty(qcinfo)
    fprintf('%s: WARNING: qcinfo is empty\n',mfilename);
    return;
  end;
  % save compiled qcinfo to output file
  mmil_struct2csv(qcinfo,parms.fname_out);
elseif nargout
  qcinfo = mmil_csv2struct(parms.fname_out);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms = mmil_args2parms(options,{...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
  ...
    'ContainerTypes',{'fsurf' 'proc_dti' 'proc_bold'},{'fsurf' 'proc_dti' 'proc_bold'},...
    'qccont_flag',true,[false true],...
    'reviewerIDs',[],[],...
    'qc_outfix','post',[],...
    'fname_out',[],[],...
    'outdir',[],[],...
    'outfix','qcinfo',[],...
    'type_sort_flag',false,[false true],...
    'verbose',true,[false true],...
    'forceflag',true,[false true],...
  ... % hidden
    'info_tags',{'RootDirs','ignore_VisitInfo_flag','user','numvec_tags'},[],...
    'stype_list',{'MPR','FLASHhi','DTI','BOLD'},[],...
  ... % qualitative scores
    'dMRI_ftype_list',{'regT1','B0warp','imgQual'},[],...
    'fMRI_ftype_list',{'regT1','B0warp','imgQual'},[],...
    'fsurf_ftype_list',{'motion','PialUnder','PialOver',...
                        'WMUnder','WMOver','Artifact'},[],...
  });
  parms.ntypes = length(parms.ContainerTypes);

  % get StudyInfo and RootDirs
  if isempty(parms.StudyInfo)
    if parms.verbose
      fprintf('%s: getting StudyInfo for %s...\n',mfilename,ProjID);
      tic;
    end;
    args = mmil_parms2args(parms,parms.info_tags);
    [parms.StudyInfo,parms.RootDirs] = MMIL_Quick_StudyInfo(ProjID,args{:});
    if parms.verbose, toc; end;
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
    fprintf('%s: %d visits and %d types to summarize\n',...
      mfilename,parms.nvisits,parms.ntypes);
  end;

  % set output file name
  if isempty(parms.fname_out)
    if isempty(parms.outdir)
      parms.outdir = sprintf('%s/MetaData/%s',parms.RootDirs.home,ProjID);
    end;
    parms.fname_out = sprintf('%s/%s_%s.csv',...
      parms.outdir,ProjID,parms.outfix);
  else
    parms.outdir = fileparts(parms.fname_out);
    if isempty(parms.outdir), parms.outdir = pwd; end;
  end;
  mmil_mkdir(parms.outdir);
  
  % check reviewerIDs
  if ~iscell(parms.reviewerIDs)
    parms.reviewerIDs = {parms.reviewerIDs};
  end;
  parms.nrev = length(parms.reviewerIDs);
  for r=1:parms.nrev
    if isnumeric(parms.reviewerIDs{r})
      parms.reviewerIDs{r} = num2str(parms.reviewerIDs{r});
    end;
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stypes = get_stypes(ContainerPath,parms)
  stypes = [];
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
  if errcode || isempty(ContainerInfo), return; end;
  stypes = fieldnames(ContainerInfo.ScanInfo);
  ind = find(~cellfun(@isempty,struct2cell(ContainerInfo.ScanInfo)));
  stypes = stypes(ind);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = get_qcinfo(s,k,parms)
  qcinfo = [];
  VisitID = parms.StudyInfo(s).VisitID;
  ContainerType = parms.ContainerTypes{k};
  if parms.verbose
    fprintf('%s: checking for %s data for %s...\n',mfilename,ContainerType,VisitID);
    tic;
  end;
  [ContainerPath,ContainerDir,ContainerRootDir] = ...
    MMIL_Get_Container(parms.RootDirs,VisitID,ContainerType);
  if parms.verbose, toc; end;
  if ~isempty(ContainerDir) && exist(ContainerPath,'dir')
    if parms.verbose
      fprintf('%s: checking %s QC info for %s...\n',mfilename,ContainerType,VisitID);
      tic;
    end;
    k = 0; % scan type index
    % get scan types from ContainerInfo
    stypes = get_stypes(ContainerPath,parms);
    for i=1:length(stypes)
      stype = stypes{i};
      if ~ismember(stype,parms.stype_list), continue; end;
      k = k + 1;
      % set list of qualitative scores
      switch ContainerType
        case 'proc_dti'
          parms.ftype_list = parms.dMRI_ftype_list;
        case 'proc_bold'
          parms.ftype_list = parms.fMRI_ftype_list;
        case 'fsurf'
          parms.ftype_list = parms.fsurf_ftype_list;
        otherwise
          parms.ftype_list = [];
      end;
      % use separate indir if qccont_flag = 1
      if parms.qccont_flag
        if ~strcmp(ContainerType,'fsurf')
          indir = sprintf('%s/%s',parms.RootDirs.qc,...
                          regexprep(ContainerDir,'PROC','QC'));
        else
          indir = sprintf('%s/%s',parms.RootDirs.qc,...
                          regexprep(ContainerDir,'FSURF','FSQC'));
        end;
      else
        indir = ContainerPath;
      end;
      % load post qcinfo
      post_qcinfo = load_post_qcinfo(indir,stype,parms);
      qcinfo(k).VisitID = VisitID;
      qcinfo(k).ContainerType = ContainerType;
      qcinfo(k).stype = stype;
      % add post QC info
      qcinfo(k).nrev = 0;
      qcinfo(k).revdisp = 0;
      if ~isempty(post_qcinfo)
        qcinfo(k).QC = 1;
        QC_last = [];
        for r=1:parms.nrev
          reviewerID = parms.reviewerIDs{r};
          if ~isempty(reviewerID)
            field_infix = ['_' reviewerID];
          else
            field_infix = '';
          end;
          fname = ['QC' field_infix];
          qcinfo(k).(fname) = mmil_getfield(post_qcinfo,fname);
          if ~isempty(qcinfo(k).(fname))
            qcinfo(k).nrev = qcinfo(k).nrev + 1;
            QC = qcinfo(k).(fname);
            if ~QC, qcinfo(k).QC = 0; end;
            if ~isempty(QC_last) && QC~=QC_last
              qcinfo(k).revdisp = 1;
            end;
            QC_last = QC;
          end;            
          fname = ['notes' field_infix];
          qcinfo(k).(fname) = mmil_getfield(post_qcinfo,fname);
          if ~isempty(qcinfo(k).(fname))
            qcinfo(k).(fname)  = regexprep(qcinfo(k).(fname),',',';');
          end;
          % add qualitative scores
          for f=1:length(parms.ftype_list)
            ftype = parms.ftype_list{f};
            fname = (['QU' field_infix '_' ftype]);
            qcinfo(k).(fname) = mmil_getfield(post_qcinfo,fname);
          end;
        end;
      else
        qcinfo(k).QC = [];
        for r=1:parms.nrev
          reviewerID = parms.reviewerIDs{r};
          if ~isempty(reviewerID)
            field_infix = ['_' reviewerID];
          else
            field_infix = '';
          end;
          qcinfo(k).(['QC' field_infix]) = [];
          qcinfo(k).(['notes' field_infix]) = [];
          % add qualitative scores
          for f=1:length(parms.ftype_list)
            ftype = parms.ftype_list{f};
            fname = (['QU' field_infix '_' ftype]);
            qcinfo(k).(fname) = [];
          end;
        end;
      end;
    end;
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = load_post_qcinfo(indir,stype,parms)
  qctype = parms.qc_outfix;
  qcinfo = [];
  for r=1:parms.nrev
    reviewerID = parms.reviewerIDs{r};
    if ~isempty(reviewerID)
      field_infix = ['_' reviewerID];
    else
      field_infix = '';
    end;
    outfix = [qctype field_infix];
    outfix = [outfix '_qcinfo'];
    fname_qcinfo = sprintf('%s/%s_%s.mat',indir,stype,outfix);
    if ~exist(fname_qcinfo,'file')
      if parms.verbose
        fprintf('%s: WARNING: qcinfo file %s not found\n',...
          mfilename,fname_qcinfo);
      end;
      continue;
    end;
    if parms.verbose
      fprintf('%s: loading qcinfo file %s...\n',...
        mfilename,fname_qcinfo);
    end;
    post_qcinfo = [];
    load(fname_qcinfo);
    if isempty(post_qcinfo)
      if parms.verbose
        fprintf('%s: WARNING: qcinfo is empty for %s (deleting)\n',...
          mfilename,fname_qcinfo);
      end;
      delete(fname_qcinfo);
      continue;
    end;
    if ~isempty(post_qcinfo)
      % get qualitative scores
      ftype_list = {};
      if isfield(post_qcinfo,'QU')
        ftype_list = parms.ftype_list;
      end; 
      QU = {};
      for f=1:length(ftype_list)
        ftype = ftype_list{f};
        g = find(strcmp(ftype,post_qcinfo.ftype_list));
        if ~isempty(g)
          qu_score = post_qcinfo.QU{g};
          if ischar(qu_score), qu_score = str2num(qu_score); end;
          QU{f} = qu_score;
        end;
      end;
      qcinfo.(['QC' field_infix]) = post_qcinfo.QC;
      qcinfo.(['notes' field_infix]) = post_qcinfo.notes;
      % add QU scores
      if ~isempty(QU)
        for f=1:length(ftype_list)
          ftype = ftype_list{f};
          qcinfo.(['QU' field_infix '_' ftype]) = QU{f};
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = concat_qcinfo(qcinfo,tmpinfo)
  if isempty(tmpinfo), return; end;
  if isempty(qcinfo)
    qcinfo = tmpinfo;
  else
    % reconcile fieldnames (add empty fields to match)
    [qcinfo,tmpinfo] = reconcile_fields(qcinfo,tmpinfo);
    qcinfo = cat(2,qcinfo,tmpinfo);    
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,B] = reconcile_fields(A,B)
  fnamesA = fieldnames(A);
  fnamesB = fieldnames(B);
  fnames = setdiff(fnamesB,fnamesA);
  for f=1:length(fnames)
    A(1).(fnames{f}) = [];
  end;
  fnames = setdiff(fnamesA,fnamesB);
  for f=1:length(fnames)
    B(1).(fnames{f}) = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

