function MMIL_postQC_dMRI_Exam(ContainerPath,varargin)
%function MMIL_postQC_dMRI_Exam(ContainerPath,[options])
%
% Purpose: Perform post-processing quality review for dMRI data
%   check registration to T1, distortion, DTI or RSI measures,
%   and fiber tract segmentation results
%
% Usage:
%  MMIL_postQC_dMRI_Exam(ContainerPath,'key1', value1,...);
%
% Required Input:
%  ContainerPath: full path of DTIPROC directory
%
% Optional Paramters :
%  'infix': processing output file infix
%     {default = 'corr_regT1'}
%  'revflag': [0|1|2] forward, reverse, or both phase-encode polarity
%     {default = 2}
%  'resT1_flag': [0|1] display images in T1 resolution
%     {default = 1}
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
% Last Mod: 07/08/16 by Don Hagler
%

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
  if ~strcmp(stype,'DTI'), continue; end;
  repeat_flag = true;
  while repeat_flag
    repeat_flag = review_images(stype,parms);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ContainerPath,options)
  parms = mmil_args2parms(options,{...
    'indir',ContainerPath,[],...
  ...
    'infix','corr_regT1',[],...
    'revflag',2,{0,1,2},...
    'resT1_flag',true,[false true],...
    'outdir',[],[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ...
    'ftype_list',{'regT1','B0warp','imgQual'},{'regT1','B0warp','imgQual','fseg'},...
    'reg_tags',{'infix','revflag'},[],...
    'fname_tcl','qc_dMRI.tcl',[],...
  });

  if isempty(parms.outdir)
    parms.outdir = parms.indir;
  elseif mmil_isrelative(parms.outdir)
    parms.outdir = [parms.indir '/' parms.outdir];
  end;

  %% todo: maybe set fname_tcl depending on options that do not yet exist (DTI FA, RSI F2, fseg)

  % check tcl file
  cmd_options = [];
  if mmil_isrelative(parms.fname_tcl)
    tcl_dir_list = {[getenv('HOME') '/tcl'],[getenv('MMPS_DIR') '/tcl']};
    for i=1:length(tcl_dir_list)
      fname_tmp = sprintf('%s/%s',tcl_dir_list{i},parms.fname_tcl);
      if exist(fname_tmp,'file')
        parms.fname_tcl = fname_tmp;
        break;
      end;
    end;
  end;
  if ~exist(parms.fname_tcl,'file')
    fprintf('%s: WARNING: file %s not found',mfilename,parms.fname_tcl);
    parms.fname_tcl = [];
  end;
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

  % review image quality for each file
  fname_qcinfo = sprintf('%s/%s_postQC_qcinfo.mat',parms.outdir,stype);
  if ~exist(fname_qcinfo,'file') || parms.forceflag
  
    % load reginfo
    [reginfo,errcode] = load_reginfo(parms);
    if errcode, return; end;
    
    % review images
    manual_qcinfo = [];
    if parms.verbose
      fprintf('%s: reviewing %s images in %s...\n',...
        mfilename,stype,parms.indir);
    end;
    
    % display images for review
    errcode = display_images(reginfo,parms);
    if errcode, return; end;

    % prompt for quality scores on multiple factors
    QU = {};
    for f=1:length(parms.ftype_list)
      ftype = parms.ftype_list{f};
      QU{f} = input_QU(ftype);
    end;
    
    % prompt for overall QC review
    QC = input_QC(stype);
    
    % prompt for notes (required if QC = 0)
    notes = input_notes(QC);

    manual_qcinfo.stype = stype;
    manual_qcinfo.ftype_list = parms.ftype_list;
    manual_qcinfo.QU = QU;
    manual_qcinfo.QC = QC;
    manual_qcinfo.notes = notes;

    % provide summary of selections
    review_selections(manual_qcinfo);

    % confirm selections
    [repeat_flag,save_flag] = confirm_selections();
    if repeat_flag || ~save_flag, return; end;
        
    % save qcinfo
    save(fname_qcinfo,'manual_qcinfo');
  else
    if parms.verbose
      fprintf('%s: %s manual qcinfo for %s already exists\n',...
        mfilename,stype,parms.indir);
    end;
  end;

  % save manual qcinfo to text file
  report_qcinfo(stype,parms)
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reginfo,errcode] = load_reginfo(parms)
  args = mmil_parms2args(parms,parms.reg_tags);
  [reginfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(parms.indir,args{:});
  if errcode
    fprintf('%s: WARNING: RegInfo not found in %s\n',...
      mfilename,ContainerPath);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = display_images(reginfo,parms)
  errcode = 0;

  % determine names of dMRI and sMRI files
  if parms.resT1_flag
    [fdir_T2,fstem_T2,fext_T2] = fileparts(reginfo.fname_T2);
    fname_T2 = sprintf('%s/%s_resT1%s',fdir_T2,fstem_T2,fext_T2);
    fname_T1 = reginfo.fname_T1;
  else
    fname_T2 = reginfo.fname_T2;
    fname_T1 = sprintf('%s/T1_resT2_DTI.mgz',parms.indir);
  end;

  %% todo: option to use DTI or RSI meas instead

  % define command string
  cmd = sprintf('tkmedit -f %s -aux %s',fname_T2,fname_T1);
  if ~isempty(parms.fname_tcl)
    cmd = sprintf('%s -tcl %s',cmd,parms.fname_tcl);
  end;

  %% todo: option to show fseg

  % execute command and catch errors    
  [s,r] = mmil_unix(cmd);
  if s,
    fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function QC = input_QC(stype)
  QC = [];
  retry_flag = true;
  while retry_flag
    retry_flag = false;
    QC = input(sprintf('accept %s? [1 | 0]: ',stype),'s');
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

function QU = input_QU(ftype)
  QU = [];
  retry_flag = true;
  while retry_flag
    retry_flag = false;
    QU = input(sprintf('rate %s? [0 | 1 | 2 | 3]: ',ftype),'s');
    if ~ismember(str2num(QU),[0:3])
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

function review_selections(manual_qcinfo)
  QC_flag = manual_qcinfo.QC;
  if QC_flag
    fprintf('accepted %s results\n',manual_qcinfo.stype);
  else
    fprintf('rejected %s results\n',manual_qcinfo.stype);
  end;
  %% todo: also display reg, B0uw, etc.
%  keyboard
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
  fname_qcinfo_mat = sprintf('%s/%s_postQC_qcinfo.mat',parms.outdir,stype);
  fname_qcinfo_txt = sprintf('%s/%s_postQC_qcinfo.txt',parms.outdir,stype);
  if ~exist(fname_qcinfo_txt,'file') || parms.forceflag
    load(fname_qcinfo_mat);
    % write summary of qcinfo to text file
    fid = fopen(fname_qcinfo_txt,'wt');
    if fid<0, error('failed to open %s for writing',fname_qcinfo_txt); end;
    for i=1:length(manual_qcinfo)
      fprintf(fid,'%s\n',manual_qcinfo(i).stype);
      %% todo: output other QC scores (reg, B0uw, etc.)
      fprintf(fid,'  QC = %d\n',manual_qcinfo(i).QC);
      if ~isempty(manual_qcinfo(i).notes)
        fprintf(fid,'  notes = %s\n',manual_qcinfo(i).notes);
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

