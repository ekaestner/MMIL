function MMIL_IcoResamp_FSRecon_Exam(ContainerDir,RootDirs,varargin)
%function MMIL_IcoResamp_FSRecon_Exam(ContainerDir,RootDirs,[options])
%
% Required Parameters:
%  ContainerDir: Input FSURF Container Dir
%  RootDirs: a struct which must contain the fields fsurf and fsico
%    specifying the container root dirs
%
% Optional Parameters:
%  'ProjID': project name added to ContainerInfo
%    {default = []}
%  'FS_version': FreeSurfer version number (e.g. 302, 305, 450)
%    {default = 450}
%  'ico': icosahedral order number (1-7)
%    {default = 7}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  08/01/08 by Don Hagler
% Last Mod: 09/29/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmd = []; err = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'ProjID',[],[],...
  'FS_version',450,[],...
  'ico',7,[1:7],...
  'forceflag',false,[false true],...
...
  'source_typestr','FSURF',[],...
  'dest_typestr','FSICO',[],...  
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(RootDirs,'fsico') | isempty(RootDirs.fsico)
  RootDirs.fsico = RootDirs.fsurf;
end;

parmdir = getenv('MMPS_PARMS');
fname_identity = [parmdir '/identity.xfm'];

fprintf('%s: FreeSurfer version: %s\n',mfilename,num2str(parms.FS_version));
fprintf('%s: matlab version: %s\n',mfilename,version);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ContainerPath = [RootDirs.fsurf '/' ContainerDir];
[ContainerInInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode, return; end;

ContainerInfo.SourceDir = ContainerDir;
ContainerInfo.ContainerType = parms.dest_typestr;
ContainerInfo.ContainerUID = '1'; % Hardcode for now, need system call to obtain unique container ID
ContainerInfo.VisitID = ContainerInInfo.VisitID;
ContainerInfo.StudyDate = ContainerInInfo.StudyDate;
ContainerInfo.StudyTime = ContainerInInfo.StudyTime;
ContainerInfo.ScanInfo = ContainerInInfo.ScanInfo;
ContainerInfo.ContainerCreationDate = datestr(now);
ContainerInfo.FS_version = parms.FS_version;

ContainerOutDir = regexprep(ContainerDir,parms.source_typestr,parms.dest_typestr);
ContainerOutPath = [RootDirs.fsico '/' ContainerOutDir];

% if exists, delete or quit
runflag = 1;
if exist(ContainerOutPath,'dir')
  if parms.forceflag
    fprintf('%s: deleting existing ico directory %s...\n',...
      mfilename,ContainerOutPath);
    cmd = ['rm -rf ' ContainerOutPath];
    fprintf('%s: cmd =\n%s\n',mfilename,cmd);
    [status,result] = mmil_unix(cmd);
    if status
      error('cmd %s failed:\n%s',cmd,result);
    end;
  else
    fprintf('%s: ico directory %s already exists\n',...
      mfilename,ContainerOutPath);
    runflag = 0;
  end;
end;

if runflag
  % replace existing talairach file with identity matrix
  fname_talairach = sprintf('%s/%s/mri/transforms/talairach.xfm',...
    RootDirs.fsurf,ContainerDir);
  fname_talairach_orig = sprintf('%s.orig',fname_talairach);
  if ~exist(fname_talairach_orig,'file')
    if exist(fname_talairach,'file')
      cmd = sprintf('cp %s %s',fname_talairach,fname_talairach_orig);
      fprintf('%s: cmd =\n%s\n',mfilename,cmd);
      [status,result] = mmil_unix(cmd);
      if status
        error('command failed:\n%s\n',result);
      end
    end
    cmd = sprintf('cp %s %s',fname_identity,fname_talairach);
    fprintf('%s: cmd =\n%s\n',mfilename,cmd);
    [status,result] = mmil_unix(cmd);
    if status
      error('command failed:\n%s\n',result);
    end
  end

  % create _ico version of subject resampled to spherical average space
  fprintf('%s: resampling FreeSurfer subject to spherical average space...\n',mfilename);
  cmd = sprintf('setenv SUBJECTS_DIR %s\n',RootDirs.fsurf);
  cmd = sprintf('%ssource $PUBSH/bin/SetUpFreeSurfer.csh %d\n',...
    cmd,parms.FS_version);
  if parms.FS_version >= 500
    cmd = sprintf('%smake_average_subject --sd-out %s --subjects %s --out %s --ico %d',...
      cmd,RootDirs.fsico,ContainerDir,ContainerOutDir,parms.ico);
  elseif parms.FS_version >= 400
    cmd = sprintf('%smake_average_subject_v4 --sd-out %s --subjects %s --out %s --ico %d',...
      cmd,RootDirs.fsico,ContainerDir,ContainerOutDir,parms.ico);
  else
    cmd = sprintf('%smake_average_subject_v3 --sd-out %s --subjects %s --out %s --ico %d',...
      cmd,RootDirs.fsico,ContainerDir,ContainerOutDir,parms.ico);
  end;
  fprintf('%s: cmd =\n%s\n',mfilename,cmd);
  [status,result] = mmil_unix(cmd);
  if status
    error('command failed:\n%s\n',result);
  end

  % restore original talairach file
  if exist(fname_talairach_orig,'file')
    fprintf('%s: restoring original talairach transform...\n',...
      mfilename);
    cmd = sprintf('mv %s %s',fname_talairach_orig,fname_talairach);
    fprintf('%s: cmd =\n%s\n',mfilename,cmd);
    [status,result] = mmil_unix(cmd);
    if status
      error('command failed:\n%s\n',result);
    end
  end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = MMIL_Save_ContainerInfo(ContainerOutPath,ContainerInfo);
if errcode ~= 0, return; end;

% add ProjID and MMPSVER to ContainerInfo
[errcode,warncode] = MMIL_Stamp_Container(ContainerOutPath,parms.ProjID);

