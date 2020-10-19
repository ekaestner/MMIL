function MMIL_Analyze_Long_DTI_Exams(ProjID,varargin)
%function MMIL_Analyze_Long_DTI_Exams(ProjID,varargin)
%
% Usage: MMIL_Analyze_Long_DTI_Exams(ProjID,'key',value,...)
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
%              (e.g. read from csv file with MMIL_Read_StudyInfo)
%               If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%     'proc' and 'fibers'
%     to specify full paths of root directories containing data containers
%     {default = []}
%  'batchrootdir': top level directory containing output batch job directories
%     {default = /home/$USER/batchdirs}
%  'batchname': name of output batchdir
%     if empty, will use mfilename
%     {default = []}
%  'required_rootdirs': List of rootdirs fields which must not be empty
%     {default= ['proc' 'fibers']}
%  'qcflag': only include subjects with StudyInfo.QC=1
%     {default=1}
%
% Optional Parameters to select the DTI data
%  'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%  'min_bval' - minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%  'nob0_flag' - [0|1] whether to exclude b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'min_ndirs': minimum number of directions to be included in tensor fit 
%    (default= 6)
%  'min_nb0': minimum number of b=0 images required for tensor calculations (default=1)
%  'infix': infix of fibers used to warp to atlas (default= 'corr_regT1')
%  'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%
% Optional Parameters:
%   'lesion_flag': [0/1] whether to create the lesion mask or not 
%     {default= 0}
%   'lesion_T1_mask_flag': [0|1] use T1 masks to create lesion mask (otherwise FA)
%     {default = 1}
%   'diff_smooth': smoothing sigma (voxels) applied to difference volumes
%     {default = 10}
%   'xcg_flag': [0|1] exclude CSF and gray-matter from fiber ROIs
%     {default = 1}
%   'masksf_flag': [0|1] exclude voxels with two or more fibers from fiber ROIs
%     {default = 0}
%   'dirprefix': prefix used for the output directory
%     {default= 'LONG'}
%   'resDTI_flag': [0/1] specify whether in resDTI 
%     {default= 0}
%   'fibers': vector of fiber numbers to include in analysis
%     {default= [101:110,115:123,133:138,141:150,1011,1021]}
%   'forceflag': run analysis even if output files exist
%     {default = 0}
%
% Created:  05/10/11 by Vijay Venkatraman  (Original Code from Don Hagler)
% Last Mod: 07/31/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin, { ...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchrootdir',[],[],...
  'batchname',[],[],...
  'required_rootdirs',{'proc','long_dti'},[],...
  'modality',{'MRI'},{'MRI','MEG','PET'},...
  'qcflag',true, [false true],...
... % parameters for selecting the data
  'snums',[],[],...
  'min_bval',1000,[],...
  'flex_flag',false,[false true],...
  'nob0_flag',false,[false true],...
  'min_ndirs',6,[],...
  'min_nb0',1,[],...
  'infix','corr_regT1',[],...
  'revflag',0,[0,1,2],...
...  
  'lesion_flag',false,[false true],...
  'lesion_T1_mask_flag',true,[false true],...
  'diff_smooth',10,[0,Inf],...
  'xcg_flag',true,[false true],...
  'masksf_flag',false,[false true],...
  'fibers',[101:110,115:123,133:138,141:150],[],...
  'resDTI_flag',false,[false true],...
  'dirprefix','LONG',[],...
...
  'forceflag',false,[false true],...  
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get project specific info

tags_excl= {'batchname','batchrootdir','StudyInfo','RootDirs','required_rootdirs','modality','qcflag', ...
  'snums','min_bval','flex_flag','nob0_flag','min_ndirs','min_nb0',...
  'infix','revflag','dirprefix'};
tags= setdiff(fieldnames(parms),tags_excl);

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;
OrigStudyInfo = StudyInfo;

if isempty(parms.batchname)
  parms.batchname = [mfilename];
end;

if ~isempty(ProjID)
 parms.batchname = [ProjID '_' parms.batchname];
end;

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify subjects with multiple visits

VisitNumbers = cell2mat({OrigStudyInfo.VisitNumber});
ind = find(VisitNumbers>1);
SubjIDs = unique(({OrigStudyInfo(ind).SubjID}));
if length(SubjIDs)==0
  error('no subjects with followup scans found -- check StudyInfo VisitNumber');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create batchdir and scriptlist
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'file')
  cmd = sprintf('rm -rf %s/*\n',batchdir);
  fprintf(1,'cmd = %s',cmd);
  system(cmd);
else
  mkdir(batchdir);
end
fid = fopen(scriptlistfname,'w'); fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create job scripts

j=1;
for i=1:length(SubjIDs)
  % check that subject has a visit 1
  ind = find(strcmp(SubjIDs{i},{OrigStudyInfo.SubjID}));
  i_baseline = find(cell2mat({OrigStudyInfo(ind).VisitNumber})==1);
  if isempty(i_baseline)
    fprintf('%s: WARNING: no VisitNumber=1 for subj %s\n',...
      mfilename,SubjIDs{i});
    continue;
  end;
  % check that longitudinal directory exists
  ContainerPath = sprintf('%s/%s_%s',...
    RootDirs.long_dti,parms.dirprefix,SubjIDs{i});
  if ~exist(ContainerPath,'dir')
    fprintf('%s: WARNING: dir %s not found -- run DTI_MMIL_Longitudinal_Setup_Exams first\n',...
      mfilename,ContainerPath);
    continue;
  end;
  % set FSPath
  FSContainerDir = OrigStudyInfo(ind(i_baseline)).fsurf;
  FSPath = [RootDirs.fsurf '/' FSContainerDir];
  if ~exist(FSPath,'dir') || isempty(FSContainerDir)
    FSPath = [];
  end;    
  
  % select study info for this subject
  StudyInfo = OrigStudyInfo(ind);
  
  % create script
  jstem = regexprep(SubjIDs{i},'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem);
  jobfname = sprintf('%s/%s.m',batchdir,jobID);
  matfname = sprintf('%s/%s.mat',batchdir,jobID);
  save(matfname,'StudyInfo','RootDirs');
  fid = fopen(jobfname,'w');
  fprintf(fid,'load(''%s'')\n',matfname);
  fprintf(fid,'MMIL_Analyze_Long_DTI_Exam(StudyInfo,...\n');
  fprintf(fid,'  ''%s'',...\n',ContainerPath);
  fprintf(fid,'  ''FSPath'',''%s'',...\n',FSPath);
  mmil_write_tags(fid,tags,parms);
  fprintf(fid,');\n');
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
  j=j+1;
end
fprintf('%%%% Now login to cluster and run this:\n',mfilename);
fprintf('    qmatjobs2 %s\n',parms.batchname);


