function MMIL_Long_Setup_Exams(ProjID,varargin)
%function MMIL_Long_Setup_Exams(ProjID,[options])
%
% Runs the MMIL_Long_Setup_Exam for Longitudinal analysis
%
% Required Input:
%   ProjID: project ID to run (e.g. 'REC_TEST') Required Input:
% 
% Optional Parameters that specify study specific information:
%  'RootDirs':
%    a struct which must contain the following fields: proc, fsurf, long
%    and may contain the following fields: orig, raw, fsico, fsclean
%    these specify the locations of data
%     {default = []}
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%      must contain these fields:
%        SubjID
%        StudyDate
%        VisitNumber
%      may contain these fields
%        proc, proc_dti, fsurf
%      if proc, proc_dti, and fsurf are unspecified, will look for Containers
%        with SubjID and StudyDate
%       (will choose first one if more than one)
%      if empty, use all subjects found in RootDirs.proc or RootDirs.proc_dti and RootDirs.fsurf
%     {default = []}
%  'batchname': name of output batchdir
%     {default = 'MMIL_Long_Setup_Exams'}
%  'required_containers': List of containers which must not be empty
%     if LongDTI_flag, 'proc_dti' also required
%     {default = {'proc'}}
%  'required_rootdirs': List of rootdirs fields which must not be empty
%     if LongDTI_flag, 'proc_dti' also required
%     {default = {'proc','fsurf','long'}
%  'qcflag': only include subjects with StudyInfo.QC=1
%     {default=1}
%  'modality': limit returned StudyInfo to sessions of this modality only
%     {default = 'MRI'}
%
% Optional Parameters to select the DTI data:
%  'DTI_min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'DTI_flex_flag': [0|1] DTI_flex scans included in tensor fit
%    {default = 0}
%  'DTI_nob0_flag': [0|1] whether to exclude b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'DTI_min_ndirs': minimum number of directions to be included in tensor fit 
%    (default = 6)
%  'DTI_min_nb0': minimum number of b=0 images required for tensor calculations 
%     (default=1)
%  'DTI_infix': infix of fibers used to warp to atlas (default= 'corr_regT1')
%  'DTI_revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%  'DTI_snums_flag': [0|1|2|3] which set of "snums" to use
%     0: use all available scans
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%    {default = 0}
%  'DTI_DT_regT1flag': [0 1 2] whether to resample DT output
%    0=not resampled, 1=registered T1, 2=resampled T1
%    {default = 0}
%
% Optional Parameters:
%  'fibers': List of fibers to be included
%     {default= [101:110,115:123,133:138,141:150]}
%  'LongDTI_flag': [0|1] setup for longitudinal DTI analysis
%     {default = 0}
%  'T1type': which type of T1 series ('MPR' or 'hiFA') to use
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%  'xcg_flag': [0|1] exclude CSF and gray-matter from fiber ROIs
%    {default = 1}
%  'masksf_flag': [0|1] exclude voxels with two or more fibers from fiber ROIs
%    {default = 0}
%  'dirprefix': prefix used for the output directory
%    {default= 'LONG'}
%  'lesion_flag': [0|1] if the study subjects have lesion, the MMIL_Long_Setup_Exam would
%     need fsurf for all visits
%    {default = 0}
%  'outext': extension of output files
%    {default = '.mgz'}
%  'fiberdir_resT1': atlas fibers in T1 resolution
%    {default = 'AtlasTrack/fiber_maps_resT1'}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%  
% Created:  04/20/11 by Vijay Venkatraman (Original by Don Hagler 03/15/09)
% Last Mod: 07/31/16 by Don Hagler
%

% todo: dilate the FA brain mask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_all= {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchname',[],[],...
  'required_containers',{'proc'},[],...
  'required_rootdirs',{'proc','fsurf','long'},[],...
  'modality','MRI',[],...
  'qcflag',false,[false true],...
...
  'DTI_snums_flag',0,[0:3],...
  'DTI_nob0_flag',false,[false true],...
  'DTI_min_ndirs',6,[],...
  'DTI_min_bval',1,[1 Inf],...
  'DTI_flex_flag',false,[false true],...
  'DTI_min_nb0',1,[],...
  'DTI_revflag',0,[0 1 2],...
  'DTI_infix','corr_regT1',[],...
  'DTI_DT_regT1flag',2,[0,1,2],...
...
  'LongDTI_flag',false,[false true],...
  'xcg_flag',true,[false true],...
  'masksf_flag',false,[false true],...
  'fibers',[101:110,115:123,133:138,141:150],[],...
  'dirprefix','LONG',[],...
  'T1type',2,[0 1 2 3]...
  'outext','.mgz',{'.mgh','.mgz'},...
  'fiberdir_resT1','AtlasTrack/fiber_maps_resT1',[],...
  'lesion_flag',false,[false true],...
  'atlasname',[],[],...
  'forceflag',false,[false true]... 
};

parms = mmil_args2parms(varargin,parms_all);

if ~isempty(parms.atlasname)
  parms.fiberdir_resT1 = [parms.fiberdir_resT1 '_' parms.atlasname];
end;

% excl_tags are fields that should not be passed to MMIL_Long_Setup_Exam
excl_tags = {'StudyInfo','RootDirs','batchname','required_containers',...
  'required_rootdirs','modality','qcflag','dirprefix','DTI_snums_flag',...
  'DTI_nob0_flag','DTI_min_ndirs','DTI_min_bval','DTI_flex_flag',...
  'DTI_min_nb0','DTI_revflag',...
  'DTI_infix','DTI_DT_regT1flag'};
tags = setdiff(fieldnames(parms),excl_tags);

% get input from ProjInfo
args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;
OrigStudyInfo = StudyInfo;

if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_all);
end;

if parms.LongDTI_flag
  parms.required_containers = union(parms.required_containers,'proc_dti');
  parms.required_rootdirs = union(parms.required_rootdirs,'proc_dti');
  if isfield(RootDirs,'long_dti') & ~isempty(RootDirs.long_dti)
    RootDirs.long = RootDirs.long_dti;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(ProjInfo,'ReconFlag')
  parms.T1type = ProjInfo.ReconFlag;
elseif isfield(ProjInfo,'STRUCT_T1type')
  parms.T1type = ProjInfo.STRUCT_T1type;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

j = 1;
for i=1:length(SubjIDs)
  SubjID = SubjIDs{i};
  % check that subject has a visit 1
  ind = find(strcmp(SubjID,{OrigStudyInfo.SubjID}));
  i_baseline = find(cell2mat({OrigStudyInfo(ind).VisitNumber})==1);
  if isempty(i_baseline)
    fprintf('%s: WARNING: no VisitNumber=1 for subj %s\n',...
      mfilename,SubjID);
    continue;
  end;

  % select study info for this subject
  StudyInfo = OrigStudyInfo(ind);

  % set output dir
  OutContainerPath = sprintf('%s/%s_%s',...
    RootDirs.long,parms.dirprefix,SubjID);

  % create script
  jstem = regexprep(SubjID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = sprintf('%s/%s.m',batchdir,jobID);
  matfname = sprintf('%s/%s.mat',batchdir,jobID);
  save(matfname,'RootDirs' ,'StudyInfo');
  fid = fopen(jobfname,'w');
  fprintf(fid,'load(''%s'')\n',matfname);
  fprintf(fid,'MMIL_Long_Setup_Exam(RootDirs,StudyInfo,...\n');
  fprintf(fid,'  ''%s'',...\n',OutContainerPath);
  fprintf(fid,'  ''snums_flag'',%d,...\n',parms.DTI_snums_flag);
  fprintf(fid,'  ''nob0_flag'',%d,...\n',parms.DTI_nob0_flag); 
  fprintf(fid,'  ''min_ndirs'',%d,...\n',parms.DTI_min_ndirs); 
  fprintf(fid,'  ''min_bval'',%d,...\n',parms.DTI_min_bval); 
  fprintf(fid,'  ''flex_flag'',%d,...\n',parms.DTI_flex_flag); 
  fprintf(fid,'  ''min_nb0'',%d,...\n',parms.DTI_min_nb0); 
  fprintf(fid,'  ''revflag'',%d,...\n',parms.DTI_revflag); 
  fprintf(fid,'  ''infix'',''%s'',...\n',parms.DTI_infix); 
  fprintf(fid,'  ''DT_regT1flag'',%d,...\n',parms.DTI_DT_regT1flag); 
  mmil_write_tags(fid,tags,parms);
  fprintf(fid,');\n');
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

% check available disk space
MMIL_Check_Usage(RootDirs.long);

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);
