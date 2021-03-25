function MMIL_Convert_Exams(ProjID,varargin)
% function MMIL_Convert_Exams(ProjID,[options])
%
% Usage:
%  MMIL_Convert_Exams(ProjID,'key1', value1,...);
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
%      fsurf, fsico
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'ContainerTypes': cell array of container types to convert
%     {default = {'proc' 'proc_dti' 'proc_bold'}}
%  'batchrootdir': top level directory containing output batch job directories
%     {default = /home/$USER/batchdirs}
%  'batchname': name of output batchdir
%     {default = 'MMIL_Convert_Exams'}
%
% Optional Parameters:
%  'outdir': output directory
%    relative to ContainerPath
%    if empty, will write output to ContainerPath
%    {default = []}
%  'infix': file name infix (e.g. [], 'corr')
%    {default = []}
%  'in_type': input file type
%    supported types: 'nii','mgh','mgz'
%    {default = 'mgz'}
%  'out_type': output file type
%    supported types: 'nii','mgh','mgz'
%    'nii' format can be used by FSL and AFNI
%    {default = 'nii'}
%  'out_orient': output slice orientation
%    if empty or omitted, keep original orientation
%      e.g. 'LPS', 'RAS', etc.
%    for use with FSL, 'LAS' may be preferred
%    {default = []}
%  'verbose': [0|1] display status messages and warnings
%    {default: 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  02/19/16 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'ContainerTypes',{'proc' 'proc_dti' 'proc_bold'},{'proc' 'proc_dti' 'proc_bold'},...
  'batchrootdir',[],[],...
  'batchname','MMIL_Convert_Exams',[],...
...
  'outdir',[],[],...
  'in_type','mgz',{'nii','mgh','mgz'},...
  'out_type','nii',{'nii','mgh','mgz'},...
  'snums_export',[],[],...
  'infix',[],[],...
  'out_orient',[],[],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
...
  'required_containers',{'raw'},[],...
});

excl_tags = {'StudyInfo','RootDirs','batchrootdir','batchname',...
  'ContainerTypes','required_containers'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create job scripts

j = 1;
for s=1:length(StudyInfo)
  % select study info for this subject visit
  VisitID = StudyInfo(s).VisitID;
  % loop over Container types
  for k=1:length(parms.ContainerTypes)
	  ContainerType = parms.ContainerTypes{k};
	  [ContainerPath,ContainerDir,ContainerRootDir] = ...
  	  MMIL_Get_Container(RootDirs,VisitID,ContainerType);
    if ~isempty(ContainerDir) && exist(ContainerPath,'dir')
      % create script
      jstem = regexprep(VisitID,'\^','_');
      jstem = jstem(1:min(20,length(jstem)));
      jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
      jobfname = [batchdir '/' jobID '.m'];
      tags = setdiff(fieldnames(parms),excl_tags);
      mmil_write_script(jobfname,'MMIL_Convert_Exam',{ContainerPath},tags,parms);
      % add to list
      fid = fopen(scriptlistfname,'a');
      fprintf(fid,'%s\n',jobID);
      fclose(fid);
    end;
  end;
end

% check available disk space
for k=1:length(parms.ContainerTypes)
	ContainerType = parms.ContainerTypes{k};
  MMIL_Check_Usage(RootDirs.(ContainerType));
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);

