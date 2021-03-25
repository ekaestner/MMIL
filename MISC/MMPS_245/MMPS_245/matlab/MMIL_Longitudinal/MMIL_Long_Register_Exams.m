function MMIL_Long_Register_Exams(ProjID,varargin)
%function MMIL_Long_Register_Exams(ProjID,[options])
%
% Registers the T1 images for Longitudinal analysis and can also do DTI images
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
% 
% Optional Parameters that specify study specific information:
%  'RootDirs':
%     struct containing the following fields: proc, fsurf, long
%     long_dti is also required if LongDTI_flag is specified.
%     these specify the locations of data
%  'StudyInfo': struct array of study information
%    (e.g. read from csv file with MMIL_Read_StudyInfo)
%    must contain these fields: SubjID, StudyDate, VisitNumber
%    may contain these fields: proc, fsurf
%    if proc and fsurf are unspecified, will look for Containers
%      with SubjID and StudyDate (will choose first one if more than one)
%    if empty, use all subjects found in RootDirs.proc and RootDirs.fsurf
%    {default = []}
%  'batchname': name of output batchdir
%    {default = 'MMIL_Long_Register_Exams'}
%  'qcflag': only include subjects with StudyInfo.QC=1 
%    {default = 1}
%
% Optional Parameters:
%  'baseflag': [0|1] register all visits to visit 1
%    otherwise, register all visits to all earlier visits
%    {default: 1}
%  'resDTI_flag': [0|1] whether to register DTI resolution images
%   otherwise use T1 resolution images
%    {default: 0}
%  'nobiasflag':[0|1] register T1 using Quarc with no bias (forward and reverse)
%    {default: 1}
%  'LongDTI_flag': [0|1] Longitudinal DTI analysis
%    {default:0}
%  'forceflag': [0|1] if output files exist, delete them and then run setup
%    {default:0}
%
% Created:  04/19/10 by Vijay Venkatraman
% Last Mod: 03/11/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'RootDirs',[],[],...
  'StudyInfo',[],[],...
  'batchname',[],[],...
  'qcflag',false,[false true],...
...
  'baseflag',true,[false true],...
  'resDTI_flag',false,[false true],...
  'nobiasflag', true, [false true],...
  'LongDTI_flag',false,[false true],...  
  'forceflag',false,[false true],...
... % undocumented:
  'dirprefix',[],[],...
  'bindir',[],[],...
  'parmdir',[],[],...
  'imgtypes',{'FA'},{'FA','ADC'},...
  'ext','.mgz',{'.mgh','.mgz'},...
  'required_containers',{'proc'},[],...
  'required_rootdirs',{'proc','fsurf'},[],...
  'modality','MRI',[],...
};
parms = mmil_args2parms(varargin,parms_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter StudyInfo or generate StudyInfo struct if none supplied
args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;
OrigStudyInfo = StudyInfo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(parms.imgtypes)
  parms.imgtypes = {parms.imgtypes};
end;

if parms.LongDTI_flag
  % add T1 to use in MMIL_Long_Register_DTI_Exam
  parms.imgtypes= [{'T1'} parms.imgtypes];
end;

if isempty(parms.parmdir)
  if parms.LongDTI_flag
    parms.parmdir = [getenv('MMPS_PARMS') '/LONGDTI'];
  else
    parms.parmdir = [getenv('MMPS_PARMS') '/QUARC'];
  end;
end;

if isempty(parms.bindir)
  parms.bindir = [getenv('MMPS_DIR') '/bin'];
end;

if parms.LongDTI_flag
  parms.required_rootdirs= [parms.required_rootdirs cellstr('long_dti')];
  else
  parms.required_rootdirs= [parms.required_rootdirs cellstr('long')];
end;

parms.dirprefix= 'LONG';
parms.required_containers = union(parms.required_containers,{'fsurf'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify subjects with multiple visits

VisitNumbers = cell2mat({StudyInfo.VisitNumber});
ind = find(VisitNumbers>1);
SubjIDs = unique(({StudyInfo(ind).SubjID}));
if length(SubjIDs)==0
  error('no subjects with followup scans found -- check StudyInfo VisitNumber');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create batchdir and scriptlist
if isempty(parms.batchname)
  parms.batchname = mfilename;
  if parms.LongDTI_flag
    parms.batchname = [parms.batchname '_DTI'];
  end;
end;
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];

if parms.LongDTI_flag
  batchdir = [RootDirs.batch '/' parms.batchname];
end;

scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s/*\n',batchdir);
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

j=1;
for i=1:length(SubjIDs)
  SubjID = SubjIDs{i};

  % check that subject has a visit 1
  ind = find(strcmp(SubjID,{StudyInfo.SubjID}));
  i_baseline = find(cell2mat({StudyInfo(ind).VisitNumber})==1);
  if isempty(i_baseline)
    fprintf('%s: WARNING: no VisitNumber=1 for subj %s\n',...
      mfilename,SubjID);
    continue;
  end;
  % check that longitudinal directory exists
  
  
  if parms.LongDTI_flag
    ContainerPath = sprintf('%s/%s_%s',...
    RootDirs.long_dti,parms.dirprefix,SubjID);
    if ~exist(ContainerPath,'dir')
      fprintf('%s: WARNING: dir %s not found -- run MMIL_Long_Setup_Exams first\n',...
      mfilename,ContainerPath);
    continue;
    end;
  else
    ContainerPath = sprintf('%s/%s_%s',...
    RootDirs.long,parms.dirprefix,SubjID);
    if ~exist(ContainerPath,'dir')
      fprintf('%s: WARNING: dir %s not found -- run MMIL_Long_Setup_Exams first\n',...
      mfilename,ContainerPath);
    continue;
    end;
  end;
  
  % sort studies by VisitNumber
  ind = find(strcmp(SubjID,{StudyInfo.SubjID}));
  tmpStudyInfo = StudyInfo(ind);
  ind = find(cell2mat({tmpStudyInfo.VisitNumber})>=1);
  if length(ind)==0
    error('no visits with valid VisitNumber');
  end;
  tmpStudyInfo = tmpStudyInfo(ind);
  VisitNumbers = cell2mat({tmpStudyInfo.VisitNumber});
  [VisitNumbers,ind] = sort(VisitNumbers);
  tmpStudyInfo = tmpStudyInfo(ind);
  if min(VisitNumbers)>1
    fprintf('%s: WARNING: no VisitNumber=1 for subj %s\n',...
      mfilename,SubjIDs{i});
    continue;
  end;

  if parms.baseflag
    nbase = 1;
  else
    nbase = length(tmpStudyInfo);
  end;

  % create scripts
  for k=1:nbase
    dirA = sprintf('%s/visit%d',ContainerPath,tmpStudyInfo(k).VisitNumber);
    if ~exist(dirA,'dir')      
      fprintf('%s: WARNING: dir %s not found -- run MMIL_Long_Setup_Exams first\n',...
        mfilename,dirA);
      continue;
    end;
    for m=k+1:length(tmpStudyInfo)
      dirB = sprintf('%s/visit%d',ContainerPath,tmpStudyInfo(m).VisitNumber);
      if ~exist(dirB,'dir')      
        fprintf('%s: WARNING: dir %s not found -- run MMIL_Long_Setup_Exams first\n',...
          mfilename,dirB);
        continue;
      end;  
      jstem = regexprep(SubjID,'\^','_');
      jstem = jstem(1:min(20,length(jstem)));
      jobID = sprintf('job_%03d_%s',j,jstem); 
      jobfname= sprintf('%s/%s.m',batchdir,jobID);          
      if parms.LongDTI_flag
        % create scripts for Longitudinal DTI analysis 
        args= {dirA,dirB};
        tags= {'imgtypes','parmdir','resDTI_flag','ext','forceflag'};
        mmil_write_script(jobfname,'MMIL_Long_Register_DTI_Exam',...
          args,tags,parms); 
      else
        % create scripts for Longitudinal T1 analysis
        args={dirA,dirB};
        tags={'nobiasflag','bindir','parmdir','ext','forceflag'};
        mmil_write_script(jobfname,'MMIL_Long_Register_Exam',...
          args,tags,parms); 
      end;
      fid = fopen(scriptlistfname,'a');
      fprintf(fid,'%s\n',jobID);
      fclose(fid);
      j=j+1;
    end;
  end;
end

% check available disk space
MMIL_Check_Usage(RootDirs.long);

fprintf('%%%% Now login to a cluster and run this:\n',mfilename);
if parms.LongDTI_flag
  fprintf('    qmatjobs %s\n',parms.batchname);
else
  fprintf('    qmatjobs3 %s\n',parms.batchname);
end;
