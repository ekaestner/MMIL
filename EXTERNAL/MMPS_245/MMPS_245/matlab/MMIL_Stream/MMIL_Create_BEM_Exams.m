function MMIL_Create_BEM_Exams(ProjID,varargin)
%function MMIL_Create_BEM_Exams(ProjID,[options])
%
% Purpose: create BEM surfaces and copy to fsurf and fsico containers
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%       (e.g. read from VisitInfo.csv file with MMIL_Read_StudyInfo)
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: raw, proc
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'batchname': name of output batchdir
%      {default = 'MMIL_Create_BEM_Exams'}
%  'BEMtype': source of BEM surfaces (e.g. 'PD','PD_NFT','T1','aseg')
%    for 'aseg' type, 'T1' is used for outer skull and outer scalp
%      and dilated FreeSurfer aseg is used for inner skull
%    for 'PD_NFT' type, surfaces from PD are used in combination with
%      FreeSurfer aseg and NFT mesh generation to create 4-shell BEM
%    {default = 'T1'}
%  'ico': icosahedral order number (1-7)
%    {default = 4}
%  'make_forceflag': [0|1] whether to overwrite existing output
%    in proc container
%    {default = 0}
%  'copy_forceflag': [0|1] whether to overwrite existing output
%    in fsurf and fsico containers
%    {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output
%    overrides both 'make_forceflag' and 'copy_forceflag'
%    {default = 0}
%
% Created:  02/16/11 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchname',[],[],...
  'BEMtype','T1',{'PD','PD_NFT','T1','aseg'},...
  'ico',4,[1:7],...
  'make_forceflag',false,[false true],...
  'copy_forceflag',false,[false true],...
  'forceflag',false,[false true],...
...
  'required_rootdirs',{'proc','fsurf','fsico'},[],...
};

tags = {'BEMtype','ico','make_forceflag','copy_forceflag','forceflag'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% create output batch directory
if isempty(parms.batchname)
  parms.batchname = mfilename;
end;
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
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

% create jobs
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
fid = fopen(scriptlistfname,'w'); fclose(fid);
j = 1;
for i=1:length(StudyInfo)
  SubjID = StudyInfo(i).SubjID;
  VisitID = StudyInfo(i).VisitID;
  ContainerDir = StudyInfo(i).proc;

  % skip if ContainerDir is missing or unspecified
  if isempty(ContainerDir)
    continue;
  end;
  % skip if this session uses a different session's structural scan
  if ~strcmp(VisitID,StudyInfo(i).STRUCT_VisitID)
    continue;
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
  fprintf(fid,'MMIL_Create_BEM_Exam(''%s'',RootDirs,...\n',ContainerDir);
  mmil_write_tags(fid,tags,parms);
  fprintf(fid,');\n');
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end;

% check available disk space
MMIL_Check_Usage(RootDirs.proc);
MMIL_Check_Usage(RootDirs.fsurf);
MMIL_Check_Usage(RootDirs.fsico);

fprintf('%%%% Now login to mmilcluster and run this:\n',parms.batchname);
fprintf('    qmatjobs %s\n',parms.batchname);

