function [cmd,err] = MMIL_Freesurfer_Recon_Exam(ContainerDir,RootDirs,varargin)
%function [cmd,err] = MMIL_Freesurfer_Recon_Exam(ContainerDir,RootDirs,[options])
%
% Required Parameters:
%  ContainerDir: Input MRIPROC Container Dir
%  RootDirs: a struct which must contain the fields proc and fsurf
%    specifying the locations of data
%
% Optional Parameters:
%  'ProjID': project name added to ContainerInfo
%    {default = []}
%  'FS_version': which version of Freesurfer to use (e.g. 305, 450, 510)
%    {default = 530}
%  'FS_full_recon_flag': [0|1] whether to do full recon or only volume aseg
%    {default = 1}
%  'FS_rescale_orig_flag': [0|1] whether to allow orig.mgz to be rescaled
%    when converting to uchar
%    if 0, orig.mgz is a uchar copy of rawavg.mgz
%       as long as highest value in input image is 255,
%       orig.mgz will not be rescaled
%    {default = 1}
%  'FS_talairach_flag': [0|1] whether to do talairach registration
%    {default = 1}
%  'FS_nu_flag': [0|1[ whether to run nonuniformity intensity correction
%    if 0, nu.mgz is simply a copy of orig.mgz
%    {default = 1}
%  'FS_nu3T_flag': [0|1] whether to use nu settings
%    optimized for 3T
%  'FS_labelV1_flag': [0|1] whether to label V1 using sulcal patterns
%    {default = 1}
%  'FS_surfsegedit_flag': [0|1] whether to edit aseg with surfaces
%    {default = 0}
%  'FS_gca_dir': full path containing substitute GCA atlas files (for aseg)
%    {default = []}
%  'FS_gca': name of substitute GCA atlas file (for aseg)
%    {default = []}
%  'FS_gca_skull': name of substitute GCA atlas file with skull (for aseg)
%    {default = []}
%  'FS_T2_flag': [0|1|2] whether to use T2-weighted volume for pial surface
%    requires FS_version >= 5.3
%    0: do not use T2-weighted volume
%    1: use T2 if available
%    2: use T2 or quit if not available
%    {default = 0}
%  'FS_fstem_T2': T2-weighted image file stem
%    {default = 'T2w_res'}
%  'STRUCT_T1type': which type of T1 series ('MPR' or 'hiFA') to use
%    0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%    {default = 2}
%  'STRUCT_BEMflag' - [0|1] copy BEM surfaces to FreeSurfer bem directory
%    optimized for 3T (version >= 500 only)
%    {default = 0}
%  'old_rootdir': root directory with existing FreeSurfer recons
%      Use this to import edited wm.mgz and brainmask.mgz files
%    {default = []}
%  'BEMtype': type of BEM surfaces ('T1','aseg','PD','PD_NFT')
%    {default = 'T1'}
%  'touchonly_flag': [0|1] whether to create touchfiles and return empty cmd
%    {default = 0}
%  'forceflag' - [0|1] whether to create script even if recon is complete
%    {default = 0}
%
% Output:
%   cmd: csh command to run for FreeSurfer recon
%   err: string containing any error messages
%
% Created:  01/01/07 by Don Hagler
% Prev Mod: 10/26/17 by Feng Xue (To allow tar as input)
% Last Mod: 11/19/17 by Feng Xue
%

%% todo: remake flags as parameters to remake all subjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmd = []; err = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'ProjID',[],[],...
  'FS_version',530,[],...
  'FS_recon_flag',true,[false true],...
  'FS_full_recon_flag',true,[false true],...
  'FS_rescale_orig_flag',true,[false true],...
  'FS_talairach_flag',true,[false true],...
  'FS_nu_flag',true,[false true],...
  'FS_nu3T_flag',false,[false true],...
  'FS_labelV1_flag',false,[false true],...
  'FS_surfsegedit_flag',false,[false true],...
  'FS_gca',[],[],...
  'FS_gca_skull',[],[],...
  'FS_gca_dir',[],[],...
  'FS_T2_flag',0,[0:2],...
  'FS_fstem_T2','T2w_res',[],...
  'STRUCT_T1type',2,[0:3],...
  'STRUCT_BEMflag',false,[false true],...
  'old_rootdir',[],[],...
  'BEMtype','T1',{'T1','aseg','PD','PD_NFT'},...
  'touchonly_flag',false,[false true],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(ContainerDir)
  err = sprintf('ContainerDir can not be empty\n');
  return;
else
  ContainerPath = sprintf('%s/%s',RootDirs.proc,ContainerDir);
  ContainerPath_tar = [];
end
switch exist(ContainerPath)
  case 0
    err = 'ContainerPath does not exist';
    return
  case 2
    fname_tar = ContainerDir;
    ContainerDir = regexprep(fname_tar,'\.tar$','');
    ContainerPath_tar = ContainerPath;
    ContainerPath = sprintf('%s/%s',RootDirs.proc,ContainerDir);
    if ~exist(sprintf('%s/ContainerInfo.mat',ContainerPath),'file');
      cmd = sprintf('tar xf %s/%s -C %s %s/ContainerInfo.mat', RootDirs.proc, fname_tar, RootDirs.proc, ContainerDir);
      [status,result] = unix(cmd);
      if status
        disp(result);
        err = sprintf('cmd %s failed',cmd);
        return
      end;
    end
end

[ContainerInInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode, return; end;

ContainerInfo.SourceDir = ContainerDir;
ContainerInfo.ContainerType = 'FSURF';
ContainerInfo.ContainerUID = '1'; % Hardcode for now, need system call to obtain unique container ID
ContainerInfo.VisitID = ContainerInInfo.VisitID;
ContainerInfo.StudyDate = ContainerInInfo.StudyDate;
ContainerInfo.StudyTime = ContainerInInfo.StudyTime;
ContainerInfo.ScanInfo = ContainerInInfo.ScanInfo;
ContainerInfo.ContainerCreationDate = datestr(now);
ContainerInfo.FS_version = parms.FS_version;

ContainerOutDir = sprintf('FSURF_%s_%s.%s_%s',...
  ContainerInfo.VisitID,ContainerInfo.StudyDate,...
  ContainerInfo.StudyTime,ContainerInfo.ContainerUID);
ContainerOutPath = sprintf('%s/%s',RootDirs.fsurf,ContainerOutDir);
if exist(sprintf('%s/ContainerInfo.mat',ContainerOutPath),'file') || exist(sprintf('%s/touch/fs.finish.all.touch',ContainerOutPath),'file'), cmd = []; return; end;
if exist(sprintf('%s.tar',ContainerOutPath),'file')
  cmd = sprintf('tar tf %s --occurrence --wildcards "%s/touch/fs.finish.all.touch"', ContainerOutPath,ContainerOutDir );
  [status, ~] = unix(cmd);
  if ~status, cmd = []; return; end;
end
% choose T1 file
if isempty(ContainerPath_tar)
  [fname_T1,errcode] = MMIL_Choose_T1(ContainerPath,'T1type',parms.STRUCT_T1type);
else
  [fname_T1,errcode] = MMIL_Choose_T1(ContainerPath_tar,'T1type',parms.STRUCT_T1type);
end
if errcode
  err = sprintf('hiFA or MPR not found in %s\n',ContainerDir);
  return;
end;
ContainerInfo.inputfile = fname_T1;

% choose T2 file
if parms.FS_T2_flag
  if isempty(ContainerPath_tar)
    [fname_T2,errcode] = MMIL_Choose_T2(ContainerPath,'fstem_T2',parms.FS_fstem_T2);
  else
    [fname_T2,errcode] = MMIL_Choose_T2(ContainerPath_tar,'fstem_T2',parms.FS_fstem_T2);
  end
  if errcode
    err = sprintf('%s not found in %s\n',parms.FS_fstem_T2,ContainerDir);
    return;
  end;
  ContainerInfo.inputfile_T2 = fname_T2;
else
  fname_T2 = [];
end;

if ~isempty(parms.old_rootdir)
  fname_bm = sprintf('%s/%s/mri/brainmask.mgz',...
    parms.old_rootdir,ContainerOutDir);
  if ~exist(fname_bm,'file')
    fprintf('%s: WARNING: input brain mask file %s not found\n',...
      mfilename,fname_bm);
    fname_bm = [];
  end

  fname_wm = sprintf('%s/%s/mri/wm.mgz',...
    parms.old_rootdir,ContainerOutDir);
  if ~exist(fname_wm,'file')
    fprintf('%s: WARNING: input white matter mask file %s not found\n',...
      mfilename,fname_wm);
    fname_wm = [];
  end
else
  fname_bm = [];
  fname_wm = [];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmd = fs_recon(ContainerOutDir,...
  'subjdir',RootDirs.fsurf,...
  'version',parms.FS_version,...
  'nocort_flag',~parms.FS_full_recon_flag,...
  'talairach_flag',parms.FS_talairach_flag,...
  'rescale_orig_flag',parms.FS_rescale_orig_flag,...
  'nu_flag',parms.FS_nu_flag,...
  'nu3T_flag',parms.FS_nu3T_flag,...
  'labelV1_flag',parms.FS_labelV1_flag,...
  'surfsegedit_flag',parms.FS_surfsegedit_flag,...
  'gca_dir',parms.FS_gca_dir,...
  'gca',parms.FS_gca,...
  'gca_skull',parms.FS_gca_skull,...
  'fname_T1',fname_T1,...
  'fname_T2',fname_T2,...
  'fname_wm',fname_wm,...
  'fname_bm',fname_bm,...
  'touchonly_flag',parms.touchonly_flag,...
  'forceflag',parms.forceflag);

if isempty(cmd)
  fprintf('%s: freesurfer recon for %s aleady finished...skipping\n',...
    mfilename,ContainerInfo.VisitID);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: creating freesurfer recon script for %s...\n',...
  mfilename,ContainerInfo.VisitID);
tmp_cmd = sprintf('#!/bin/csh\n\n');
tmp_cmd = sprintf('%ssetenv SUBJECTS_DIR %s\n',tmp_cmd,RootDirs.fsurf);
tmp_cmd = sprintf('%ssource $PUBSH/bin/SetUpFreeSurfer.csh %d\n',...
  tmp_cmd,parms.FS_version);

tmp_cmd = sprintf('%s\n\n%s\n\n',tmp_cmd,cmd);
cmd = tmp_cmd;

% check status with matlab function
cmd = sprintf('%smatlab -nojvm -nosplash -r "try, MMIL_Get_FSReconStatus(''%s'',%d,1);',...
  cmd,ContainerOutPath,parms.FS_version);
if parms.STRUCT_BEMflag
  % create BEM surfaces from aseg
  cmd = sprintf('%s MMIL_Make_BEM_Surfs(''%s'',''%s'');',...
    cmd,ContainerPath,ContainerOutPath);
  % copy BEM surfaces to FreeSurfer recon dir
  cmd = sprintf('%s MMIL_Copy_BEM_Surfs(''%s'',''%s'',''%s'');',...
    cmd,ContainerPath,ContainerOutPath,parms.BEMtype);
end;
cmd = sprintf('%s catch, disp(lasterr), end, exit"\n',cmd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save ContainerInfo
if ~exist(ContainerOutPath,'dir')
  mkdir(ContainerOutPath);
end
errcode = MMIL_Save_ContainerInfo(ContainerOutPath,ContainerInfo);
if errcode ~= 0, return; end;

% add ProjID and MMPSVER to ContainerInfo
[errcode,warncode] = MMIL_Stamp_Container(ContainerOutPath,parms.ProjID);

