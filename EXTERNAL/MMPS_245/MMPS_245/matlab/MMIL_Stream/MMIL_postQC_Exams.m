function MMIL_postQC_Exams(ProjID,varargin)
%function MMIL_postQC_Exams(ProjID,[options])
%
% Purpose: Perform post-processing quality review
%   for diffusion MRI data
%     check registration to T1, distortion, DTI or RSI measures,
%     and fiber tract segmentation results
%   for functional MRI data
%     check registration to T1, distortion, and image quality
%   for structural MRI data
%     check FreeSurfer cortical surface reconstruction
%
% Usage:
%  MMIL_postQC_Exams(ProjID,'key1', value1,...);
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
%    It may contain qc field
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'ContainerTypes': cell array of container types to check
%     {default = {'fsurf' 'proc_dti' 'proc_bold'}}
%
% Optional Parameters for dMRI data:
%  'dMRI_image_type': type of image used for dMRI QC
%     only applies to proc_dti container type
%     allowed values: 'T2','FA','ND'
%     'T2': average b=0 image
%     'FA': fractional anisotropy from DTI
%     'ND': neurite density from RSI
%    {default = 'FA'}
%  'dMRI_fseg_flag': [0|1] whether to include fiber segmentation overlay
%     only applies to proc_dti container type
%    {default = 0}
%  'DTI_infix': processing output file infix
%    {default = 'corr_regT1'}
%  'DTI_revflag': [0|1|2] forward, reverse, or both phase-encode polarity
%    {default = 2}
%  'DTI_flex_flag': [0|1] DTI_flex scans included in tensor fit
%    {default = 0}
%  'DT_outfix': string attached to DTcalc output files
%    {default = []}
%  'RSI_outfix': string attached to RSIcalc output files
%    {default = []}
%  'RSI_min_bval': minimum b-value used for RSI calculations
%    {default = 2000}
%
% Optional Parameters for dMRI data:
%  'BOLD_infix': processing output file infix
%    {default = 'corr_resBOLD'}
%
% Optional Parameters:
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
%  'qc_outfix': string  added to output qcinfo files
%     allows for reviews with different criteria
%     {Default = 'post'}
%  'qccont_flag': [0|1] put output in separate QC containers
%    if RootDirs.qc not specified, will be set to 0
%    {default = 1}
%  'outdir': output directory (if qcont_flag = 0)
%    relative to ContainerPath; if empty, will write output to ContainerPath
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
%    {default: 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  07/11/16 by Don Hagler
% Prev Mod: 10/30/17 by Don Hagler
% Last Mod: 11/20/17 by Don Hagler
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
    'ContainerTypes',{'fsurf' 'proc_dti' 'proc_bold'},{'fsurf' 'proc_dti' 'proc_bold'},...
  ...
    'dMRI_image_type','FA',{'T2','FA','ND'},...
    'dMRI_fseg_flag',false,[false true],...
    'DTI_infix','corr_regT1',[],...
    'DTI_revflag',2,{0,1,2},...
    'DTI_flex_flag',true,[false,true],...
    'DT_outfix',[],[],...
    'RSI_outfix',[],[],...
    'RSI_min_bval',2000,[],...
    'BOLD_infix','corr_resBOLD',[],...
  ...
    'reviewerID',[],[],...
    'alt_reviewerIDs',[],[],...
    'max_nrev',2,[1,100],...
    'qc_outfix','post',[],...
    'qccont_flag',true,[false true],...
    'outdir',[],[],...
    'type_sort_flag',false,[false true],...
    'date_range',[],[],...
    'precheck_flag',true,[false true],...
    'fname_qcinfo',[],[],...
    'revdisp_flag',0,[0:2],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % qualitative scores
    'dMRI_ftype_list',{'regT1','B0warp','imgQual'},[],...
    'fMRI_ftype_list',{'regT1','B0warp','imgQual'},[],...
    'fsurf_ftype_list',{'motion','PialOver','WMUnder',...
                        'inhomogeneity','Artifact'},[],...
  ...
    'fname_tcl',[],[],...
    'fname_tcl_fsurf','qc_recon.tcl',[],...
    'fname_tcl_FA','qc_FA.tcl',[],...
    'fname_tcl_T2','qc_T2.tcl',[],...
    'fname_tcl_fseg','qc_fseg.tcl',[],...
    'stype_list',{'MPR','FLASHhi','DTI','BOLD'},[],...
    'fname_colorlut',[],[],...
    'fsurf_opts','nu.mgz lh.white -aux brainmask.mgz -aux-surface rh.white -segmentation aparc+aseg.mgz',[],...
    'fsurf_needed_files',{'/mri/aparc+aseg.mgz'...
                          '/mri/nu.mgz' ...
                          '/mri/brainmask.mgz'...
                          '/surf/rh.pial'...
                          '/surf/lh.pial'...
                          '/surf/rh.white'...
                          '/surf/lh.white'},[],...
    'DT_outdir','DTcalc',[],...
    'RSI_outdir','RSIcalc',[],...
  ... % tags
    'postQC_tags',{'reviewerID','alt_reviewerIDs','max_nrev','qc_outfix',...
                   'outdir','resT1_flag','dMRI_image_type','dMRI_fseg_flag',...
                   'DTI_infix','DTI_revflag','DTI_flex_flag',...
                   'DT_outfix','RSI_outfix','RSI_min_bval',...
                   'BOLD_infix','verbose','forceflag',...
                   'dMRI_ftype_list','fMRI_ftype_list','fsurf_ftype_list',...
                   'fname_tcl','fname_tcl_fsurf','fname_tcl_FA',...
                   'fname_tcl_T2','fname_tcl_fseg',...
                   'stype_list','fname_colorlut',...
                   'fsurf_opts','fsurf_needed_files',...
                   'DT_outdir','RSI_outdir'},[],...
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
      suffix = sprintf('_%s_qcinfo',parms.qc_outfix);
    else
      suffix = sprintf('_%s_%s_qcinfo',parms.qc_outfix,parms.reviewerID);
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
        if strcmp(parms.ContainerTypes{k},'fsurf')
          visitlist = cellfun(@(x) regexprep(x,'\w+FSURF_',''),proclist, 'UniformOutput',false);
        else
          visitlist = cellfun(@(x) regexprep(x,'\w+PROC_',''),proclist, 'UniformOutput',false);
        end;
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
  k = regexp(VisitIDs,'_(?<StudyDate>\d{8})','names');
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

function review_exam(s,k,parms)
  VisitID = parms.VisitIDs{s};
  ContainerType = parms.ContainerTypes{k};
  if parms.verbose
    fprintf('%s: preparing to review exam for %s %s...\n',...
    	mfilename,VisitID,ContainerType);
  end;
  [ContainerPath,ContainerDir,ContainerRootDir] = ...
    MMIL_Get_Container(parms.RootDirs,VisitID,ContainerType);
  if isempty(ContainerDir) || ~exist(ContainerPath)
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
    if ismember(ContainerType,{'proc_dti' 'proc_bold'})
      outdir = regexprep(ContainerDir,'PROC','QC');
    else
      outdir = regexprep(ContainerDir,'FSURF','FSQC');
    end;      
    parms.outdir = sprintf('%s/%s',parms.RootDirs.qc,outdir);
  end;
  if parms.verbose
    fprintf('%s: reviewing exam for %s %s...\n',...
      mfilename,VisitID,ContainerType);
  end;
  args = mmil_parms2args(parms,parms.postQC_tags);
  MMIL_postQC_Exam(ContainerPath,args{:});
return;

