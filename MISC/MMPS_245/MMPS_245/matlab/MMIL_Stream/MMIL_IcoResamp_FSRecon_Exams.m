function MMIL_IcoResamp_FSRecon_Exams(ProjID,varargin)
%function MMIL_IcoResamp_FSRecon_Exams(ProjID,[options])
%
% Required Parameteres:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%    Use this option to limit processing to specific subjects
%    Should have a field 'VisitID' that indicates the name of the
%       input data directory (i.e. in orig data directory)
%    If VisitID field is missing, will use SubjID field instead
%    May also specify subject-specific parameter values
%     but note that these will override ProjInfo and command line input
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: fsurf, fsico
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'FS_version': FreeSurfer version number (e.g. 302, 305, 450)
%    If empty, will look for $FREESURFER_VER
%    If not found, will use MMIL_Freesurfer_Recon_Exam default (450)
%    {default = []}
%  'ico': icosahedral order number (1-7)
%    {default = 7}
%  'batchname': name of output batchdir
%    {default = 'MMIL_IcoResamp_FSRecon_Exams'}
%  'logflag': [0|1] whether to create log file
%    {default = 1}
%  'logfile': [0|1] name of output log file
%    if not specified, will create log file in RootDirs.home/logs
%     with name = {ProjID}_MMIL_Freesurfer_Recon_Exams.log
%    {default = []}
%  'check_complete_flag': [0|1] whether to require that recon is complete
%    {default = 1}
%  'forceflag': [0|1] whether to overwrite even if IcoResamp is complete
%    {default = 0}
%
% Created:  04/05/09 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'ProjID',ProjID,[],...
...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'FS_version',[],[],...
  'ico',7,[1:7],...
  'batchname',[],[],...
  'logflag',true,[false true],...
  'logfile',[],[],...
  'check_complete_flag',true,[false true],...
  'forceflag',false,[false true],...
... % undocumented
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'required_rootdirs',{'fsurf','fsico'},[],...
  'source_typestr','FSURF',[],...
  'dest_typestr','FSICO',[],...  
...
  'resamp_tags',{'ProjID','FS_version','ico','forceflag',...
                 'source_typestr','dest_typestr'},[],...
};
parms = mmil_args2parms(varargin,parms_filter);

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

if parms.ico~=7
  parms.dest_typestr = sprintf('%s%d',parms.dest_typestr,parms.ico);
end;

% set FreeSurfer version
if isempty(parms.FS_version)
  parms.FS_version = str2num(getenv('FREESURFER_VER'));
end;

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
    error('cmd %s failed:\n%s',cmd,result);
  end;
end;
mmil_mkdir(batchdir);

% create log file
if parms.logflag
  if isempty(parms.logfile)
    logdir = sprintf('%s/logs',RootDirs.home);
    if ~exist(logdir,'dir'), mmil_mkdir(logdir); end;
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
fpaths = {};
for i=1:length(StudyInfo)
  SubjID = StudyInfo(i).SubjID;
  VisitID = StudyInfo(i).VisitID;
  ContainerDir = StudyInfo(i).fsurf;
  ContainerPath = [RootDirs.fsurf '/' ContainerDir];
  % do not create multiple jobs for same fpath
  if ismember(ContainerPath,fpaths), continue; end;
  fpaths{end+1} = ContainerPath;
  % skip FS recon is missing
  if isempty(ContainerDir)
    tmp_SubjID = StudyInfo(i).SubjID;
    if isempty(tmp_SubjID), tmp_SubjID = num2str(i); end; 
    fprintf('%s: WARNING: skipping subj %s (recon missing)\n',...
      mfilename,tmp_SubjID);
    if parms.logflag
      fprintf(flog,'WARNING: skipping subj %s (recon missing)\n',...
        tmp_SubjID);
    end;
    continue;
  end;
  % check status of recon
  if parms.check_complete_flag
    [status,message] = MMIL_Get_FSReconStatus(ContainerPath,parms.FS_version);
    if status ~= 2 & status ~= 5
      fprintf('%s: WARNING: skipping %s (recon incomplete)\n',...
        mfilename,ContainerDir);
      if parms.logflag
        fprintf(flog,'WARNING: skipping %s (recon incomplete)\n',...
          ContainerDir);
      end;
      continue;
    end;
  end;
  % check if ico recon already exists and needs to be forced
  ContainerDir_ico = regexprep(ContainerDir,...
    parms.source_typestr,parms.dest_typestr);
  ContainerPath_ico = [RootDirs.fsico '/' ContainerDir_ico];
  tmp_parms = parms;

  if exist(ContainerPath_ico,'dir') & ~parms.forceflag
    % check if and when key files from IcoResamp were created
    % check if native space recon has been remade recently
    for h=1:length(parms.hemilist)
      fname = sprintf('%s/surf/%s.white.avg.area.mgh',...
        ContainerPath_ico,parms.hemilist{h});
      tmp = dir(fname);
      if isempty(tmp)
        fprintf('%s: NOTE: forcing rerun for %s (IcoResamp incomplete)\n',...
          mfilename,ContainerDir);
        if parms.logflag
          fprintf(flog,'NOTE: forcing rerun for %s (IcoResamp incomplete)\n',...
            ContainerDir);
        end;
        tmp_parms.forceflag = 1;
        break;
      else
        % look for .remade* files in native space container
        flist = dir(sprintf('%s/touch/remade*',ContainerPath));
        for i=1:length(flist)
          if flist(i).datenum>tmp.datenum % compare dates
            fprintf('%s: NOTE: forcing rerun for %s (recon was remade on %s)\n',...
              mfilename,ContainerDir,datestr(flist(i).datenum,'dd-mmm-yyyy'));
            if parms.logflag
              fprintf(flog,'NOTE: forcing rerun for %s (recon was remade on %s)\n',...
                ContainerDir,datestr(flist(i).datenum,'dd-mmm-yyyy'));
            end;
            tmp_parms.forceflag = 1;
            break;
          end;
        end;
      end;
    end;
    if ~tmp_parms.forceflag, continue; end; % already complete
  end;
  fprintf('%s: SUCCESS: creating job for %s\n',mfilename,ContainerDir);
%  if isempty(ContainerDir), keyboard; end;
  if parms.logflag
    fprintf(flog,'SUCCESS: creating job for %s\n',ContainerDir);
  end;

  % create script
  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = [batchdir '/' jobID '.m'];
  matfname = [batchdir '/' jobID '.mat'];
  save(matfname,'RootDirs');
  fid = fopen(jobfname,'w');
  fprintf(fid,'load(''%s'')\n',matfname);
  fprintf(fid,'MMIL_IcoResamp_FSRecon_Exam(''%s'',RootDirs,...\n',ContainerDir);
  mmil_write_tags(fid,parms.resamp_tags,tmp_parms);
  fprintf(fid,');\n');
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end;
fclose(flog);

if j>1
  fprintf('%%%% Now login to a cluster and run this:\n',parms.batchname);
  fprintf('    qmatjobs %s\n',parms.batchname);
else
  fprintf('%s: no jobs created\n',mfilename);
end;

