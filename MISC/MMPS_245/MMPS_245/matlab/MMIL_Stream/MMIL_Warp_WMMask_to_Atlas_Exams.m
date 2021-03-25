function MMIL_Warp_WMMask_To_Atlas_Exams(ProjID,varargin)
%function MMIL_Warp_WMMask_To_Atlas_Exams(ProjID,[options])
%
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
%    must include these fields: fsurf
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'FS_version': FreeSurfer version number (e.g. 450, 530)
%    if empty, will look for $FREESURFER_VER
%    used to check whether recon is complete
%    {default = []}
%  'batchname': name of output batchdir
%    {default = 'MMIL_Warp_WMMask_To_Atlas_Exams'}
%  'logflag': [0|1] whether to create log file
%    {default = 1}
%  'logfile': [0|1] name of output log file
%    if not specified, will create log file in RootDirs.home/logs
%     with name = ProjID_MMIL_Freesurfer_Recon_Exams.log
%    {default = []}
%  'check_complete_flag': [0|1] whether to require that recon is complete
%    {default = 1}
%  'outdir': output directory relative to FreeSurfer container dir
%    {default = 'atlas'}
%  'forceflag': [0|1] overwrite existing files
%    {default = 0}
%
% Created:  10/21/14 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'ProjID',ProjID,[],...
...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'FS_version',[],[],...
  'batchname',[],[],...
  'logflag',true,[false true],...
  'logfile',[],[],...
  'check_complete_flag',true,[false true],...
  'outdir','atlas',[],...
  'forceflag',false,[false true],...
... % undocumented
  'required_rootdirs',{'fsurf'},[],...
  'test_files',{'wmmask.mgz','nu_dctReg2Atlas.mat','nu_atlas.mgz',...
                'wmmask_atlas.mgz','orig_atlas.mgz',...
                'wm_atlas.mgz','brainmask_atlas.mgz'},[],...
  'ribbon_file','ribbon.mgz',[],...
  'warp_tags',{'outdir','forceflag'},[],...
};
parms = mmil_args2parms(varargin,parms_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if ~isempty(ProjInfo)
  if isfield(ProjInfo,'FS_version') && isempty(parms.FS_version)
    parms.FS_version = ProjInfo.FS_version;
  end;
end;

% set FreeSurfer version
if isempty(parms.FS_version)
  parms.FS_version = str2num(getenv('FREESURFER_VER'));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create jobs
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
fid = fopen(scriptlistfname,'w'); fclose(fid);
j = 1;
fpaths = {};
for s=1:length(StudyInfo)
  VisitID = StudyInfo(s).VisitID;
  FSContainerDir = StudyInfo(s).fsurf;
  FSContainerPath = [RootDirs.fsurf '/' FSContainerDir];
  % do not create multiple jobs for same fpath
  if ismember(FSContainerPath,fpaths), continue; end;
  fpaths{end+1} = FSContainerPath;

  % skip FS recon is missing
  if isempty(FSContainerDir)
    fprintf('%s: WARNING: skipping subj %s (recon missing)\n',...
      mfilename,VisitID);
    if parms.logflag
      fprintf(flog,'WARNING: skipping subj %s (recon missing)\n',...
        VisitID);
    end;
    continue;
  end;

  % check status of recon
  if parms.check_complete_flag
    [status,message] = MMIL_Get_FSReconStatus(FSContainerPath,parms.FS_version);
    if status ~= 2 & status ~= 5
      fprintf('%s: WARNING: skipping %s (recon incomplete)\n',...
        mfilename,VisitID);
      if parms.logflag
        fprintf(flog,'WARNING: skipping %s (recon incomplete)\n',...
          VisitID);
      end;
      continue;
    end;
  end;

  % check input file
  fname_ribbon = sprintf('%s/mri/%s',FSContainerPath,parms.ribbon_file);
  n_ribbon = dir(fname_ribbon);
  if isempty(n_ribbon)
    fprintf('%s: WARNING: skipping subj %s (ribbon file missing)\n',...
      mfilename,VisitID);
    if parms.logflag
      fprintf(flog,'WARNING: skipping subj %s (recon missing)\n',...
        VisitID);
    end;
    continue;
  end;
    
  % check if output already exists and needs to be forced
  full_outdir = sprintf('%s/%s',FSContainerPath,parms.outdir);
  tmp_parms = parms;
  if exist(full_outdir,'dir') & ~parms.forceflag
    % check if output files are older than ribbon file
    all_exist = 1;
    for i=1:length(parms.test_files)
      fname_test = sprintf('%s/%s',full_outdir,parms.test_files{i});
      n_test = dir(fname_test);
      if isempty(n_test)
        all_exist = 0;
        break;
      elseif n_ribbon(1).datenum > n_test(1).datenum % compare dates
        fprintf('%s: WARNING: forcing rerun for %s (ribbon file created on %s)\n',...
          mfilename,VisitID,datestr(n_ribbon(1).datenum,'dd-mmm-yyyy'));
        if parms.logflag
          fprintf(flog,'WARNING: forcing rerun for %s (ribbon file created on %s)\n',...
            VisitID,datestr(n_ribbon(1).datenum,'dd-mmm-yyyy'));
        end;
        tmp_parms.forceflag = 1;
        break;
      end;
    end;
    if all_exist && ~tmp_parms.forceflag % already complete
      fprintf('%s: NOTE: skipping %s (warp to atlas complete)\n',...
        mfilename,VisitID);
      if parms.logflag
        fprintf(flog,'NOTE: skipping %s (warp to atlas complete)\n',...
          VisitID);
      end;
      continue;
    end;
  end;

  fprintf('%s: SUCCESS: creating job for %s\n',mfilename,VisitID);
  if parms.logflag
    fprintf(flog,'SUCCESS: creating job for %s\n',VisitID);
  end;

  % create script
  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = [batchdir '/' jobID '.m'];
  fid = fopen(jobfname,'w');
  fprintf(fid,'mmil_warp_wmmask_to_atlas(''%s'',...\n',FSContainerPath);
  mmil_write_tags(fid,parms.warp_tags,tmp_parms);
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
  fprintf('    qmatjobs2 %s\n',parms.batchname);
else
  fprintf('%s: no jobs created\n',mfilename);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

