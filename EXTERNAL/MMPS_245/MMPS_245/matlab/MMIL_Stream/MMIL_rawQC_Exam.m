function MMIL_rawQC_Exam(ContainerPath,varargin)
%function MMIL_rawQC_Exam(ContainerPath,[options])
%
% Purpose: Perform manual quality review of raw data
%
% Usage:
%  MMIL_autoQC_Exam(ContainerPath,'key1', value1,...);
%
% Required Input:
%  ContainerPath: full path of MRIPROC, BOLDPROC, or DTIPROC directory
%
% Optional Paramters :
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
%    provide full path or relative to ContainerPath
%    if empty, will write output to ContainerPath
%    {default = []}
%  'verbose': [0|1] display status messages and warnings
%    {default: 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  03/18/16 by Don Hagler
% Prev Mod: 01/17/17 by Don Hagler
% Last Mod: 05/16/17 by Don Hagler
%

%% todo: add "hidden" sMRI_ftype_list, etc. to header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ContainerPath,varargin);

[parms,errcode] = get_info(parms);
if errcode || ~parms.ntypes, return; end;

mmil_mkdir(parms.outdir);

% display images, prompt for quality review
for i=1:parms.ntypes
  stype = parms.stypes{i};
  repeat_flag = true;
  while repeat_flag
    repeat_flag = review_images(stype,parms);
  end;
  % save manual qcinfo to text file
  report_qcinfo(stype,parms);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ContainerPath,options)
  parms = mmil_args2parms(options,{...
    'indir',ContainerPath,[],...
  ...
    'reviewerID',[],[],...
    'alt_reviewerIDs',[],[],...
    'max_nrev',2,[1,100],...
    'outdir',[],[],...
    'export_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden
    'sMRI_ftype_list',{'motion'},[],...
    'dMRI_ftype_list',[],[],...
    'fMRI_ftype_list',[],[],...
    'sMRI_stype_list',{'MPR','FLASHhi','XetaT2'},[],...
    'dMRI_stype_list',{'DTI'},[],...
    'fMRI_stype_list',{'BOLD'},[],...
    'planestrings',{'HOR','SAG','COR'},[],...
    'infix_list',{'tv','SAG'},[],...
  });

  if isempty(parms.outdir)
    parms.outdir = parms.indir;
  elseif mmil_isrelative(parms.outdir)
    parms.outdir = [parms.indir '/' parms.outdir];
  end;

  % construct outfix
  parms.outfix = 'raw';
  if ~isempty(parms.reviewerID)
    if isnumeric(parms.reviewerID)
      parms.reviewerID = num2str(parms.reviewerID);
    end;
    parms.outfix = [parms.outfix '_' parms.reviewerID];
  end;
  parms.outfix = [parms.outfix '_qcinfo'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = get_info(parms)
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(parms.indir);
  if errcode || isempty(ContainerInfo), return; end;
  stypes = fieldnames(ContainerInfo.ScanInfo);
  ind = find(~cellfun(@isempty,struct2cell(ContainerInfo.ScanInfo)));
  parms.stypes = stypes(ind);
  parms.ntypes = length(ind);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function repeat_flag = review_images(stype,parms)
  repeat_flag = false;

  % check if output file exists
  fname_qcinfo = sprintf('%s/%s_%s.mat',...
    parms.outdir,stype,parms.outfix);
  if exist(fname_qcinfo,'file') && ~parms.forceflag
    if parms.verbose
      fprintf('%s: %s raw qcinfo in %s already exists\n',...
        mfilename,stype,parms.outdir);
    end;
    return;
  end;
    
  % compare number of reviewers to max_nrev
  if parms.max_nrev > 0 && ~parms.forceflag
    if ~isempty(parms.alt_reviewerIDs)
      % find number of mat files matching the pattern
      nrev = 0;
      for i=1:length(parms.alt_reviewerIDs)
        altID = parms.alt_reviewerIDs{i};
        fname_qcinfo = sprintf('%s/%s_raw_%s_qcinfo.mat',...
          parms.outdir,stype,altID);
        if  exist(fname_qcinfo,'file')
          nrev = nrev + 1;
        end;
      end;
    else
      % find all mat files matching the pattern
      flist = dir(sprintf('%s/%s_raw_*_qcinfo.mat',parms.outdir,stype));
      nrev = length(flist);
    end;
    if nrev >= parms.max_nrev
      if parms.verbose
        fprintf('%s: %s raw qcinfo in %s already exists (%d reviewers)\n',...
          mfilename,stype,parms.outdir,nrev);
      end;
      return;
    end;
  end;
  
  % review image quality for each file
  fname_qcinfo = sprintf('%s/%s_%s.mat',...
    parms.outdir,stype,parms.outfix);
  if ~exist(fname_qcinfo,'file') || parms.forceflag

    % load auto qcinfo
    auto_qcinfo = load_auto_qcinfo(stype,parms);
    nscans = length(auto_qcinfo);
    if ~nscans, return; end;

    % review each scan
    raw_qcinfo = [];
    for i=1:nscans
      fstem = auto_qcinfo(i).fstem;
      if parms.verbose
        fprintf('%s: reviewing images for %s in %s...\n',...
          mfilename,fstem,parms.outdir);
      end;

      % display images for review
      errcode = display_images(fstem,parms,0);
      if errcode, continue; end;

      % prompt for quality review
      QC = input_QC_init(fstem);
      
      % display more images
      if QC==-1
        errcode = display_images(fstem,parms,1);
        if errcode, continue; end;

        % prompt for quality review
        QC = input_QC_final(fstem);
      end;

      % prompt for quality scores on multiple factors
      QU = {}; ftype_list = {};
      switch stype
        case parms.sMRI_stype_list
          ftype_list = parms.sMRI_ftype_list;
        case parms.dMRI_stype_list
          ftype_list = parms.dMRI_ftype_list;
        case parms.fMRI_stype_list
          ftype_list = parms.fMRI_ftype_list;
      end;
      for f=1:length(ftype_list)
        ftype = ftype_list{f};
        QU{f} = input_QU(ftype);
      end;

      % prompt for notes (required if QC = 0)
      notes = input_notes(QC);

      raw_qcinfo(i).fstem = fstem;
      raw_qcinfo(i).QC = QC;
      raw_qcinfo(i).QU = QU;
      raw_qcinfo(i).ftype_list = ftype_list;
      raw_qcinfo(i).notes = notes;
    end;
    if isempty(raw_qcinfo), return; end;
    
    % provide summary of selections
    review_selections(raw_qcinfo);

    % confirm selections
    [repeat_flag,save_flag] = confirm_selections();
    if repeat_flag || ~save_flag, return; end;
        
    % save qcinfo
    save(fname_qcinfo,'raw_qcinfo');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = load_auto_qcinfo(stype,parms)
  auto_qcinfo = [];
  fname_qcinfo = sprintf('%s/%s_auto_qcinfo.mat',parms.outdir,stype);
  if ~exist(fname_qcinfo,'file')
    if parms.verbose
      fprintf('%s: WARNING: %s auto qcinfo file in %s not found\n',...
        mfilename,stype,parms.outdir);
    end;
    return;
  end;
  load(fname_qcinfo);
  if parms.verbose && isempty(auto_qcinfo)
    fprintf('%s: WARNING: %s auto qcinfo in %s is empty\n',...
      mfilename,stype,parms.outdir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = display_images(fstem,parms,more_flag)
  errcode = 0;
  if ~exist('more_flag','var') || isempty(more_flag), more_flag = 0; end;
  
  % check tif files exist
  fnamelist = {};
  if more_flag
    flist = dir(sprintf('%s/%s*.tif',parms.outdir,fstem));
    fnamelist = {flist.name};
  else
    for i=1:length(parms.infix_list)
      infix = parms.infix_list{i};
      flist = dir(sprintf('%s/%s*%s*.tif',parms.outdir,fstem,infix));
      fnamelist = union(fnamelist,{flist.name});
    end;
  end;
  nimages = length(fnamelist);
  if ~nimages
    if parms.verbose
      fprintf('%s: WARNING: images for %s not found (was autoQC run?)\n',...
        mfilename,fstem);
    end;
    errcode = 1;
    return;
  end;

  % display tif files
  cmd = 'eog';
  for j=1:nimages
    cmd = sprintf('%s %s/%s',cmd,parms.outdir,fnamelist{j});
  end;
  [s,r] = mmil_unix(cmd);
  if s
    fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function QC = input_QC_init(fstem)
  QC = [];
  retry_flag = true;
  while retry_flag
    retry_flag = false;
    QC = input(sprintf('accept %s? [1 | 0 | m]: ',fstem),'s');
    switch QC
      case '0'
        QC = 0;
      case '1'
        QC = 1;
      case {'m','M'}
        QC = -1;
      otherwise
        retry_flag = true;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function QU = input_QU(ftype)
  QU = [];
  retry_flag = true;
  while retry_flag
    retry_flag = false;
    QU = input(sprintf('rate %s? [0 | 1 | 2 | 3]: ',ftype),'s');
    switch QU
      case {'0','1','2','3'}
        QU = str2num(QU);
      otherwise
        retry_flag = true;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function QC = input_QC_final(fstem)
  QC = [];
  retry_flag = true;
  while retry_flag
    retry_flag = false;
    QC = input(sprintf('accept %s? [1 | 0]: ',fstem),'s');
    switch QC
      case '0'
        QC = 0;
      case '1'
        QC = 1;
      otherwise
        retry_flag = true;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function notes = input_notes(QC)
  notes = [];
  retry_flag = true;
  while retry_flag
    if QC
      notes = input('notes (optional): ','s');
    else
      notes = input('notes (required): ','s');
    end;
    if ~isempty(notes) || QC
      retry_flag = false;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function review_selections(raw_qcinfo)
  % display overall QC
  QC_flags = [raw_qcinfo.QC];
  ind_accept = find(QC_flags==1);
  ind_reject = find(QC_flags==0);
  if ~isempty(ind_accept)
    fstemlist = {raw_qcinfo(ind_accept).fstem};
    fprintf('accepted scans: %s\n',...
      sprintf('%s ',fstemlist{:}));
  end;
  if ~isempty(ind_reject)
    fstemlist = {raw_qcinfo(ind_reject).fstem};
    fprintf('rejected scans: %s\n',...
      sprintf('%s ',fstemlist{:}));
  end;

  % display qualitative QC scores
  nscans = length(raw_qcinfo);
  for i=1:nscans
    fstem = raw_qcinfo(i).fstem;
    for f=1:length(raw_qcinfo(i).ftype_list)
      ftype = raw_qcinfo(i).ftype_list{f};
      QU = raw_qcinfo(i).QU{f};
      fprintf('%s: %s = %d\n',fstem,ftype,QU);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [repeat_flag,save_flag] = confirm_selections()
  repeat_flag = false;
  save_flag = true;
  retry_flag = true;
  while retry_flag
    retry_flag = false;
    check_flag = input(sprintf('confirm, repeat, or skip? [c | r | s]: '),'s');    
    switch check_flag
      case 'c'
        fprintf('selection confirmed\n');
      case 'r'
        fprintf('selection discarded, repeating...\n');
        repeat_flag = true;
      case 's'
        fprintf('selection discarded, skipping...\n');
        save_flag = false;
      otherwise
        retry_flag = true;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function report_qcinfo(stype,parms)
  fname_qcinfo_mat = sprintf('%s/%s_%s.mat',parms.outdir,stype,parms.outfix);
  if ~exist(fname_qcinfo_mat,'file'), return; end;
  fname_qcinfo_txt = sprintf('%s/%s_%s.txt',parms.outdir,stype,parms.outfix);
  if ~exist(fname_qcinfo_txt,'file') || parms.forceflag
    load(fname_qcinfo_mat);
    % write summary of qcinfo to text file
    fid = fopen(fname_qcinfo_txt,'wt');
    if fid<0, error('failed to open %s for writing',fname_qcinfo_txt); end;
    for i=1:length(raw_qcinfo)
      fprintf(fid,'%s\n',raw_qcinfo(i).fstem);
      % display overall QC
      fprintf(fid,'  QC = %d\n',raw_qcinfo(i).QC);
      % display qualitative QC scores
      for f=1:length(raw_qcinfo(i).ftype_list)
        ftype = raw_qcinfo(i).ftype_list{f};
        QU = raw_qcinfo(i).QU{f};
        fprintf(fid,'  %s = %d\n',ftype,QU);
      end;
      % display notes
      if ~isempty(raw_qcinfo(i).notes)
        fprintf(fid,'notes = %s\n',raw_qcinfo(i).notes);
      end;
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = get_fstem(fname)
  [tmp,fstem] = fileparts(fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = get_scannum(fname,stype)
  fstem = get_fstem(fname);
  n = regexp(fstem,[stype '(?<snum>\d+)'],'names');
  if isempty(n)
    error('failed to get scan number from file name %s',fname);
  end;
  s = str2num(n.snum);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

