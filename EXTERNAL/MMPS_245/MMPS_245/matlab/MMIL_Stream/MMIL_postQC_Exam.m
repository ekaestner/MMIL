function MMIL_postQC_Exam(ContainerPath,varargin)
%function MMIL_postQC_Exam(ContainerPath,[options])
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
%  MMIL_postQC_Exam(ContainerPath,'key1', value1,...);
%
% Required Input:
%  ContainerPath: full path of procesed container
%   valid types: DTIPROC, BOLDPROC, or FSRECON
%
% Optional Parameters for dMRI data:
%  'dMRI_image_type': type of image used for dMRI QC
%     only applies to proc_dti container type
%     allowed values: 'T2','FA','ND'
%     'T2': mean b=0 image
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
%
% Optional Parameters for fMRI data:
%  'BOLD_infix': processing output file infix
%    {default = 'corr_resBOLD'}
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
%  'qc_outfix': string  added to output qcinfo files
%     allows for reviews with different criteria
%     {Default = 'post'}
%  'outdir': output directory
%    provide full path or relative to ContainerPath
%    if empty, will write output to ContainerPath
%    {default = []}
%  'resT1_flag': [0|1] display dMRI or fMRI images in T1 resolution
%     {default = 1}
%  'verbose': [0|1] display status messages and warnings
%    {default: 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  03/18/16 by Don Hagler
% Prev Mod: 10/30/17 by Don Hagler
% Last Mod: 11/20/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: allow T2 from MRIPROC
%%       (requires implementation of registration to T1 for 3D T2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ContainerPath,varargin);

[parms,errcode] = get_info(parms);
if errcode || ~parms.ntypes, return; end;

% set fname_tcl, depending on ContainerType and options
parms = select_tcl_script(parms);

% display images, prompt for quality review
for i=1:parms.ntypes
  stype = parms.stypes{i};
  if ~ismember(stype,parms.stype_list), continue; end;
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
    'reviewerID',[],[],...
    'alt_reviewerIDs',[],[],...
    'max_nrev',2,[1,100],...
    'qc_outfix','post',[],...
    'outdir',[],[],...
    'resT1_flag',true,[false true],...
    'dMRI_image_type','FA',{'T2','FA','N2','ND','F2','FD'},...
    'dMRI_fseg_flag',false,[false true],...
    'DTI_infix','corr_regT1',[],...
    'DTI_revflag',2,{0,1,2},...
    'DTI_flex_flag',true,[false,true],...
    'DT_outfix',[],[],...
    'RSI_outfix',[],[],...
    'BOLD_infix','corr_resBOLD',[],...
    'verbose',false,[false true],...
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
  ...
    'stype_list',{'MPR','FLASHhi','DTI','BOLD'},[],...
  ...
    'flex_flag',true,[false,true],...
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
    'DT_min_bval',2000,[],...
    'RSI_outdir','RSIcalc',[],...
    'RSI_min_bval',2000,[],...
  ...
    'DT_fstem_tags',{'snums','infix','revflag','min_bval','flex_flag',...
                     'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
    'RSI_fstem_tags',{'snums','infix','revflag','min_bval','flex_flag',...
                     'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
  });

  if isempty(parms.outdir)
    parms.outdir = parms.indir;
  elseif mmil_isrelative(parms.outdir)
    parms.outdir = [parms.indir '/' parms.outdir];
  end;
  
  % construct outfix
  parms.outfix = parms.qc_outfix;
  if ~isempty(parms.reviewerID)
    if isnumeric(parms.reviewerID)
      parms.reviewerID = num2str(parms.reviewerID);
    end;
    parms.outfix = [parms.outfix '_' parms.reviewerID];
  end;
  parms.outfix = [parms.outfix '_qcinfo'];

  if parms.dMRI_fseg_flag
    if ~ismember('fseg',parms.dMRI_ftype_list)
      parms.dMRI_ftype_list{end+1} = 'fseg';
    end;
    if isempty(parms.fname_fscolorlut)
      MMPS_parms = getenv('MMPS_PARMS');
      parms.fname_fscolorlut = [MMPS_parms '/MMIL_FSColorLUT.txt'];
    end;
    if ~exist(parms.fname_fscolorlut,'file')
      error('FreeSurfer color lookup table file %s not found',...
        parms.fname_fscolorlut);
    end;
  end;

  parms.FS_version = str2num(getenv('FREESURFER_VER'));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set fname_tcl
function parms = select_tcl_script(parms)
  if isempty(parms.fname_tcl)
    switch parms.ContainerType
      case 'DTIPROC'
        if parms.dMRI_fseg_flag
          fname_tcl = parms.fname_tcl_fseg;
        elseif ismember(parms.dMRI_image_type,{'FA','ND'})
          fname_tcl = parms.fname_tcl_FA;
        else
          fname_tcl = parms.fname_tcl_T2;
        end;
      case 'BOLDPROC'
        fname_tcl = parms.fname_tcl_T2;
      case 'FSURF'
        fname_tcl = parms.fname_tcl_fsurf;
    end;
  else
    fname_tcl = parms.fname_tcl;
  end;
  % find tcl file
  cmd_options = [];
  if mmil_isrelative(fname_tcl)
    tcl_dir_list = {[getenv('HOME') '/tcl'],[getenv('MMPS_DIR') '/tcl']};
    for i=1:length(tcl_dir_list)
      fname_tmp = sprintf('%s/%s',tcl_dir_list{i},fname_tcl);
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
  parms.ContainerType = ContainerInfo.ContainerType;
  switch parms.ContainerType
    case 'DTIPROC'
      parms.ftype_list = parms.dMRI_ftype_list;
    case 'BOLDPROC'
      parms.ftype_list = parms.fMRI_ftype_list;
    case 'FSURF'
      parms.ftype_list = parms.fsurf_ftype_list;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function repeat_flag = review_images(stype,parms)
  repeat_flag = false;

  % check if output file exists
  fname_qcinfo = sprintf('%s/%s_%s.mat',...
    parms.outdir,stype,parms.outfix);
  if exist(fname_qcinfo,'file') && ~parms.forceflag
    if parms.verbose
      fprintf('%s: %s %s qcinfo for %s already exists\n',...
        mfilename,stype,parms.qc_outfix,parms.indir);
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
        fname_qcinfo = sprintf('%s/%s_%s_%s_qcinfo.mat',...
          parms.outdir,stype,parms.qc_outfix,altID);
        if  exist(fname_qcinfo,'file')
          nrev = nrev + 1;
        end;
      end;
    else
      % find all mat files matching the pattern
      flist = dir(sprintf('%s/%s_%s_*_qcinfo.mat',parms.outdir,stype,parms.qc_outfix));
      nrev = length(flist);
    end;
    if nrev >= parms.max_nrev
      if parms.verbose
        fprintf('%s: %s %s qcinfo for %s already exists (%d reviewers)\n',...
          mfilename,stype,parms.qc_outfix,parms.indir,nrev);
      end;
      return;
    end;
  end;

  mmil_mkdir(parms.outdir);

  % review image quality
  fname_qcinfo = sprintf('%s/%s_%s.mat',...
    parms.outdir,stype,parms.outfix);
  if ~exist(fname_qcinfo,'file') || parms.forceflag

    % load reginfo
    switch parms.ContainerType
      case 'DTIPROC'
        [reginfo,errcode] = load_DTI_reginfo(parms);
        if errcode, return; end;
      case 'BOLDPROC'
        [reginfo,errcode] = load_BOLD_reginfo(parms);
        if errcode, return; end;
      otherwise
        reginfo = [];
    end;
  
    % review images
    post_qcinfo = [];
    if parms.verbose
      fprintf('%s: reviewing %s images in %s...\n',...
        mfilename,stype,parms.indir);
    end;
    
    % display images for review
    [pid,errcode] = display_images(reginfo,parms);
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

    post_qcinfo.stype = stype;
    post_qcinfo.ftype_list = parms.ftype_list;
    post_qcinfo.QU = QU;
    post_qcinfo.QC = QC;
    post_qcinfo.notes = notes;

    % provide summary of selections
    review_selections(post_qcinfo);

    % confirm selections
    [repeat_flag,save_flag] = confirm_selections();

    % close image viewing program
    kill_process(pid);
    
    if repeat_flag || ~save_flag, return; end;

    % save qcinfo
    save(fname_qcinfo,'post_qcinfo');
  else
    if parms.verbose
      fprintf('%s: %s %s qcinfo for %s already exists\n',...
        mfilename,stype,parms.qc_outfix,parms.indir);
    end;
  end;

  % save manual qcinfo to text file
  report_qcinfo(stype,parms)
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kill_process(pid)
  cmd = sprintf('kill %d',pid);
  [s,r] = unix(cmd);
  if s
    fprintf('%s: WARNING: failed to close image viewer with pid %d\n',...
      mfilename,pid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reginfo,errcode] = load_DTI_reginfo(parms)
  tparms = [];
  tparms.revflag = parms.DTI_revflag;
  tparms.infix = parms.DTI_infix;
  args = mmil_parms2args(tparms);
  [reginfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo(parms.indir,args{:});
%  if errcode
%    fprintf('%s: WARNING: RegInfo not found in %s\n',...
%      mfilename,parms.indir);
%  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reginfo,errcode] = load_BOLD_reginfo(parms)
  tparms = [];
  tparms.infix = parms.BOLD_infix;
  args = mmil_parms2args(tparms);
  [reginfo,fname_reg,errcode] = BOLD_MMIL_Load_RegInfo(parms.indir,args{:});
%  if errcode
%    fprintf('%s: WARNING: RegInfo not found in %s\n',...
%      mfilename,parms.indir);
%  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pid,errcode] = display_images(reginfo,parms)
  pid = []; errcode = 0;
  switch parms.ContainerType
    case 'DTIPROC'
      [pid,errcode] = display_dMRI_images(reginfo,parms);
    case 'BOLDPROC'
      [pid,errcode] = display_fMRI_images(reginfo,parms);
    case 'FSURF'
      [pid,errcode] = display_fsurf_images(reginfo,parms);
    otherwise
      error('unsupported ContainerType: %s',parms.ContainerType);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function [pid,errcode] = display_dMRI_images(reginfo,parms)
  pid = []; errcode = 0;
  % determine names of dMRI and sMRI files
  if parms.resT1_flag
    [fdir_T2,fstem_T2,fext_T2] = fileparts(reginfo.fname_T2);
    fname_T2 = sprintf('%s/%s_resT1%s',fdir_T2,fstem_T2,fext_T2);
    fname_T1 = reginfo.fname_T1;
  else
    fname_T2 = reginfo.fname_T2;
    fname_T1 = sprintf('%s/T1_resT2_DTI.mgz',parms.indir);
  end;
  % check whether to use DTI or RSI measure
  switch parms.dMRI_image_type
    case 'FA'
      parms.infix = parms.DTI_infix;
      parms.revflag = parms.DTI_revflag;
      parms.min_bval = parms.DT_min_bval;
      parms.outfix = parms.DT_outfix;
      parms.outdir = parms.DT_outdir;
      args = mmil_parms2args(parms,parms.DT_fstem_tags);
      [fstem,parms.snums] = DTI_MMIL_Set_DT_fstem(parms.indir,args{:});
      if isempty(fstem)
        errcode = 1;
        return;
      end;
      if parms.resT1_flag
        fname_T2 = sprintf('%s_%s_resT1.mgz',fstem,parms.dMRI_image_type);
      else
        fname_T2 = sprintf('%s_%s.mgz',fstem,parms.dMRI_image_type);
      end;
    case {'N2','ND','F2','FD'}
      parms.infix = parms.DTI_infix;
      parms.revflag = parms.DTI_revflag;
      parms.min_bval = parms.RSI_min_bval;
      parms.outfix = parms.RSI_outfix;
      parms.outdir = parms.RSI_outdir;
      args = mmil_parms2args(parms,parms.RSI_fstem_tags);
      [fstem,parms.snums] = DTI_MMIL_Set_RSI_fstem(parms.indir,args{:});
      if isempty(fstem)
        errcode = 1;
        return;
      end;
      if parms.resT1_flag
        fname_T2 = sprintf('%s_%s_f0_resT1.mgz',fstem,parms.dMRI_image_type);
      else
        fname_T2 = sprintf('%s_%s.mgz',fstem,parms.dMRI_image_type);
      end;
  end;
  % check whether files exist
  if ~exist(fname_T1,'file')
    fprintf('%s: WARNING: file %s not found\n',...
      mfilename,fname_T1);
    errcode = 1;
    return;
  end;
  if ~exist(fname_T2,'file')
    fprintf('%s: WARNING: file %s not found\n',...
      mfilename,fname_T2);
    errcode = 1;
    return;
  end;
  % set fseg file name
  if parms.dMRI_fseg_flag
    if parms.resT1_flag
      fname_fseg = sprintf('%s/AtlasTrack/fseg_resT1_xcg.mgz',parms.indir);
    else
      fname_fseg = sprintf('%s/AtlasTrack/fseg_resDTI_xcg.mgz',parms.indir);
    end;
    if ~exist(fname_fseg,'file')
      fprintf('%s: WARNING: file %s not found\n',...
        mfilename,fname_T2);
      errcode = 1;
      return;
    end;
  end;  
  % define command string
  cmd = sprintf('tkmedit -f %s -aux %s',fname_T2,fname_T1);
  % add fseg to command
  if parms.dMRI_fseg_flag
    cmd = sprintf('%s -segmentation %s %s',cmd,fname_fseg,parms.fname_fscolorlut);
  end;
  % add tcl script
  if ~isempty(parms.fname_tcl)
    cmd = sprintf('%s -tcl %s',cmd,parms.fname_tcl);
  end;
  % execute command and catch errors    
  [pid,errcode] = exec_cmd(cmd,'tkmedit.bin',fname_T2);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pid,errcode] = display_fMRI_images(reginfo,parms)
  pid = []; errcode = 0;
  % determine names of fMRI and sMRI files
  if parms.resT1_flag
    [fdir_T2,fstem_T2,fext_T2] = fileparts(reginfo.fname_T2);
    fname_T2 = sprintf('%s/%s_resT1%s',fdir_T2,fstem_T2,fext_T2);
    fname_T1 = reginfo.fname_T1;
  else
    fname_T2 = reginfo.fname_T2;
    fname_T1 = sprintf('%s/T1_resT2_DTI.mgz',parms.indir);
  end;
  % define command string
  cmd = sprintf('tkmedit -f %s -aux %s',fname_T2,fname_T1);
  % add tcl script
  if ~isempty(parms.fname_tcl)
    cmd = sprintf('%s -tcl %s',cmd,parms.fname_tcl);
  end;
  % execute command and catch errors    
  [pid,errcode] = exec_cmd(cmd,'tkmedit.bin',fname_T2);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pid,errcode] = display_fsurf_images(reginfo,parms)
  pid = []; errcode = 0;
  % check for color lut file
  fname_lut = [getenv('FREESURFER_HOME') '/FreeSurferColorLUT.txt'];
  if ~exist(fname_lut,'file')
    error('file %s not found',fname_lut);
  end;
  % get root and container dir
  [FSRootDir,FSContainerDir,fext] = fileparts(parms.indir);
  FSContainerDir = [FSContainerDir,fext];
  % check that recon is complete
  [status,message] = MMIL_Get_FSReconStatus(parms.indir,parms.FS_version);
  if status ~= 2 & status ~= 5
    fprintf('%s: WARNING: skipping %s (recon incomplete)\n',...
      mfilename,FSContainerDir);
    errcode = 1;
    return;
  end;
  % check that required files exist
  all_exist = 1;
  for nf = 1:length(parms.fsurf_needed_files);
    fname_test = sprintf('%s/%s',parms.indir,parms.fsurf_needed_files{nf});
    if ~exist(fname_test,'file')
      all_exist = 0;
      break;
    end
  end
  if ~all_exist
    fprintf('%s: WARNING: skipping %s (missing required files)\n',...
      mfilename,FSContainerDir);
    errcode = 1;
    return;
  end
  % define command string
  cmd = sprintf('setenv SUBJECTS_DIR %s; tkmedit %s %s %s',...
    FSRootDir,FSContainerDir,parms.fsurf_opts,fname_lut);  
  % add tcl script
  if ~isempty(parms.fname_tcl)
    cmd = sprintf('%s -tcl %s',cmd,parms.fname_tcl);
  end;
  % execute command and catch errors    
  [pid,errcode] = exec_cmd(cmd,'tkmedit.bin',FSContainerDir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pid,errcode] = exec_cmd(cmd,cmd_name,cmd_tag)
  pid = []; errcode = 0;
  cmd = [cmd ' &'];
  [s,r] = mmil_unix(cmd);
  if s
    fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
    errcode = 1;
    return;
  end;
  % get information about running process
  cmd = sprintf('ps -ef | grep ''%s'' | grep ''%s''',cmd_name,cmd_tag);
  [s,r] = mmil_unix(cmd);
  if s
    fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
    errcode = 1;
    return;
  end;
  % extract pid
  R = mmil_splitstr(r,'\n');
  for i=1:length(R)
    x = R{i};
    if ~isempty(regexp(x,'grep')), continue; end;
    k = regexp(x,'[\w\+]\s+(?<pid>\d+)\s+\d+\s+','names');
    if ~isempty(k)
      pid = cat(1,pid,str2num(k.pid));
    end;
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
    switch QU
      case {'0','1','2','3'}
        QU = str2num(QU);
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

function review_selections(post_qcinfo)
  % display qualitative QC scores: reg, B0uw, etc.
  for f=1:length(post_qcinfo.ftype_list)
    ftype = post_qcinfo.ftype_list{f};
    QU = post_qcinfo.QU{f};
    fprintf('%s = %d\n',ftype,QU);
  end;
  % display overall QC
  QC_flag = post_qcinfo.QC;
  if QC_flag
    fprintf('accepted %s results\n',post_qcinfo.stype);
  else
    fprintf('rejected %s results\n',post_qcinfo.stype);
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
    for i=1:length(post_qcinfo)
      fprintf(fid,'%s\n',post_qcinfo(i).stype);
      % display overall QC
      fprintf(fid,'  QC = %d\n',post_qcinfo(i).QC);
      % display qualitative QC scores
      for f=1:length(post_qcinfo(i).ftype_list)
        ftype = post_qcinfo(i).ftype_list{f};
        QU = post_qcinfo(i).QU{f};
        fprintf(fid,'  %s = %d\n',ftype,QU);
      end;
      % display notes
      if ~isempty(post_qcinfo(i).notes)
        fprintf(fid,'notes = %s\n',post_qcinfo(i).notes);
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

