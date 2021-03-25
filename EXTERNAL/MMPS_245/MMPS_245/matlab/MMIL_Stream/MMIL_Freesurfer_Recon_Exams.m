function MMIL_Freesurfer_Recon_Exams(ProjID,varargin)
%function MMIL_Freesurfer_Recon_Exams(ProjID,[options])
%
% Required Parameteres:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%      must contain this field:
%        SubjID
%      may contain these fields:
%        VisitID or StudyDate
%        VisitNumber
%        proc
%      if proc is unspecified, will look for Containers
%        with VisitID, or with SubjID and StudyDate
%       (will choose first one if more than one)
%      if StudyInfo is empty, use all exams found in RootDirs.proc
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: raw, proc
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'FS_version': which version of Freesurfer to use (e.g. 305, 450, 510)
%    if empty, will use FREESURER_VER environment variable
%    {default = []}
%  'FS_recon_flag': [0|1] whether to run FreeSurfer for this project
%    {default = 1}
%  'FS_full_recon_flag': [0|1] whether to do full recon or only volume aseg
%    {default = 1}
%  'FS_followup_flag': [0|1] whether to do recons on sessions with VisitNumber>1
%     0 = baseline only, 1 = baseline and follow-ups (done independently)
%     StudyInfo must have VisitNumber field to limit to baseline only
%    {default = 0}
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
%    {default = 0}
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
%  'STRUCT_BEMflag': [0|1] copy BEM surfaces to FreeSurfer bem directory
%    optimized for 3T (version >= 500 only)
%    {default = 0}
%  'BEMtype': type of BEM surfaces ('T1','PD','aseg')
%    {default = 'T1'}
%  'old_rootdir': root directory with existing FreeSurfer recons
%      Use this to import edited wm.mgz and brainmask.mgz files
%    {default = []}
%  'batchname': name of output batchdir
%    {default = 'MMIL_Freesurfer_Recon_Exams'}
%  'logflag': [0|1] whether to create log file
%    {default = 1}
%  'logfile': [0|1] name of output log file
%     if not specified, will create log file in RootDirs.home/logs
%      with name = ProjID_MMIL_Freesurfer_Recon_Exams.log
%    {default = []}
%  'touchonly_flag': [0|1] whether to create touchfiles only
%    {default = 0}
%  'forceflag': [0|1] whether to create script even if recon is complete
%    {default = 0}
%
% Created:  04/03/09 by Don Hagler
% Prev Mod: 05/11/17 by Don Hagler
% Last Mod: 10/26/17 by Feng Xue (To allow tar as input)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = { ...
  'ProjID',ProjID,[],...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'FS_version',[],[],...
  'FS_recon_flag',true,[false true],...
  'FS_full_recon_flag',true,[false true],...
  'FS_followup_flag',false,[false true],...
  'FS_rescale_orig_flag',true,[false true],...
  'FS_talairach_flag',true,[false true],...
  'FS_nu_flag',true,[false true],...
  'FS_nu3T_flag',false,[false true],...
  'FS_labelV1_flag',false,[false true],...
  'FS_surfsegedit_flag',false,[false true],...
  'FS_gca_dir',[],[],...
  'FS_gca',[],[],...
  'FS_gca_skull',[],[],...
  'FS_T2_flag',0,[0:2],...
  'FS_fstem_T2','T2w_res',[],...
  'STRUCT_T1type',2,[0:3],...
  'STRUCT_BEMflag',false,[false true],...
  'old_rootdir',[],[],...
  'BEMtype','T1',{'PD','T1','aseg'},...
  'batchname','MMIL_Freesurfer_Recon_Exams',[],...
  'logflag',true,[false true],...
  'logfile',[],[],...
  'touchonly_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'required_rootdirs',{'proc','fsurf'},[],...
... % for backward compatibility with ProjInfo.csv files
  'FS_ReconFlag',[],[],...
  'FS_FullReconFlag',[],[],...
};
parms = mmil_args2parms(varargin,parms_filter);

% excl_tags are fields that should not be passed to MMIL_Process_Exam
excl_tags = {'RootDirs' 'StudyInfo' 'batchname' ...
  'logflag' 'logfile' ...
  'FS_recon_flag' 'FS_followup_flag'...
  'FS_ReconFlag' 'FS_FullReconFlag'...
  'required_rootdirs'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

% set FreeSurfer version if not already set
if isempty(parms.FS_version)
  parms.FS_version = str2num(getenv('FREESURFER_VER'));
end;

% for backward compatibility
if ~isempty(parms.FS_ReconFlag)
  parms.FS_recon_flag = parms.FS_ReconFlag;
end;
if ~isempty(parms.FS_FullReconFlag)
  parms.FS_full_recon_flag = parms.FS_FullReconFlag;
end;

if ~parms.FS_recon_flag
  error('Freesurfer recon not allowed for project %s',ProjID);
end;

tags = setdiff(fieldnames(parms),excl_tags);
args = mmil_parms2args(parms,tags);

% check for VisitNumber
if ~parms.FS_followup_flag
  if ~isfield(StudyInfo,'VisitNumber')
    StudyInfo(1).VisitNumber = [];
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output batch directory
if isempty(parms.batchname)
  parms.batchname = mfilename;
end;
if ~isempty(parms.ProjID)
  parms.batchname = [parms.ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s\n',batchdir);
  fprintf('cmd = %s',cmd);
  [status,result] = unix(cmd);
  if status
    fprintf('%s: WARNING: cmd %s failed:\n%s',mfilename,cmd,result);
  end;
end;
mmil_mkdir(batchdir);

% create log file
if parms.logflag
  if isempty(parms.logfile)
    logdir = sprintf('%s/logs',RootDirs.home);
    mmil_mkdir(logdir);
    parms.logfile = [logdir '/' parms.batchname '.log'];
  end;
  flog = fopen(parms.logfile,'wt');
  if flog==-1
    error('unable to create log file %s',parms.logfile);
  end;
end;

% create jobs
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
fid = fopen(scriptlistfname,'w'); fclose(fid);
j = 1;
for i=1:length(StudyInfo)
  SubjID = StudyInfo(i).SubjID;
  VisitID = StudyInfo(i).VisitID;
  ContainerDir = StudyInfo(i).proc;
  tarball = dir(sprintf('%s/*%s*.tar',RootDirs.proc,VisitID));
  if isempty(tarball)
    if isempty(ContainerDir)
      if parms.logflag
        fprintf(flog,'WARNING: skipping recon for %s with missing proc container\n',...
          VisitID);
      end;
      continue;
    end
  end;

  if ~parms.FS_followup_flag &...
     ~isempty(StudyInfo(i).VisitNumber) & StudyInfo(i).VisitNumber~=1
    if parms.logflag
      fprintf(flog,'WARNING: skipping recon for %s with VisitNumber = %d\n',...
        ContainerDir,StudyInfo(i).VisitNumber);
    end;
    continue;
  end;
  if ~strcmp(VisitID,StudyInfo(i).STRUCT_VisitID)
    if parms.logflag
      fprintf(flog,'WARNING: skipping recon for %s with STRUCT_VisitID = %s\n',...
        ContainerDir,StudyInfo(i).STRUCT_VisitID);
    end;
    continue;
  end;
  if ~isempty(tarball)
    [cmd,err] = ...
      MMIL_Freesurfer_Recon_Exam(tarball.name,RootDirs,args{:});
  else
    [cmd,err] = ...
      MMIL_Freesurfer_Recon_Exam(ContainerDir,RootDirs,args{:});
  end
  if ~isempty(err)
    if parms.logflag
      fprintf(flog,'ERROR: unable to run freesurfer recon for %s:\n%s\n',...
        ContainerDir,err);
    end;
    continue;
  end;
  if isempty(cmd)
    if parms.logflag
      fprintf(flog,'WARNING: recon already complete for %s\n',...
        ContainerDir);
    end;
    continue;
  else
    if parms.logflag
      fprintf(flog,'SUCCESS: Creating recon script for %s\n',...
        ContainerDir);
    end;
  end;
  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = sprintf('%s/%s.csh',batchdir,jobID);
  fid = fopen(jobfname,'w');
  fprintf(fid,'%s',cmd);
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end;
if parms.logflag
  fclose(flog);
end;

MMIL_Check_Usage(RootDirs.fsurf);

fprintf('%%%% Now login to a cluster and run this:\n',parms.batchname);
if parms.FS_version>=500
  fprintf('    qcshjobs2 %s\n',parms.batchname);
else
  fprintf('    qcshjobs %s\n',parms.batchname);
end;

