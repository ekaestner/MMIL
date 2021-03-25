function MMIL_Analyze_Long_Exams(ProjID,varargin)
% function MMIL_Analyze_Long_Exams(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
% 
% Optional Parameters that specify project specific information
%  'RootDirs':
%    struct containing the following fields: proc, fsurf, long
%    these specify the locations of data
%  'StudyInfo': struct array of study information
%    (e.g. read from csv file with MMIL_Read_StudyInfo)
%    must contain these fields: SubjID, StudyDate, VisitNumber
%    may contain these fields: proc, fsurf
%    if proc and fsurf are unspecified, will look for Containers
%      with SubjID and StudyDate (will choose first one if more than one)
%    if empty, use all subjects found in RootDirs.proc and RootDirs.fsurf
%    {default = []}
%  'batchname': name of output batchdir
%    {default = 'MMIL_Analyze_Long_Exams'}
%  'qcflag': [0|1] only create jobs for StudyInfo entries with QC=1
%    {default: 0}
%
% Optional Parameters:
%  'baseflag': [0|1] register all visits to visit 1
%    otherwise, register all visits to all earlier visits
%    {default: 1}
%  'sphsmoothsteps': vector of surface smoothing steps (on sphere)
%    slope of FWHM vs. sqrt(N) is ~1.25 for fsaverage
%    (FWHM = full-width-half-max smoothing kernel
%        N = number of smoothing steps)
%     with 705, approx FWHM (mm) = 33.2
%    {default = 705}
%  'mask_midbrain_flag', [0|1] whether to mask out mid brain and other
%     cortical regions marked "unknown" (masking done before smoothing)
%    {default = 0}
%  'outdir': where to place output files in relative to Long ContainerDir
%    {default = 'analysis'}
%  'aseg_flag': use aseg ROIs (non-cortical, volume segmentation)
%    {default = 1}
%  'aseg_roigroups_flag': [0|1] use masks for groups of aseg roi codes
%    includes 'WholeBrain', 'LatVentricles', and 'AllVentricles'
%    {default = 1}
%  'subhippo_flag': use subdivided hippocampal ROIs
%    {default = 1}
%  'aparc_flag': use aparc ROIs (cortical surface parcellation)
%    {default = 1}
%  'nobiasflag':[0|1] register T1 using Quarc with no bias (forward and reverse)
%    {default: 1}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Created:  03/01/10 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'RootDirs',[],[],...
  'StudyInfo',[],[],...
  'batchname',[],[],...
  'qcflag',false,[false true],...
... % other options:
  'baseflag',true,[false true],...
  'sphsmoothsteps',705,[0,Inf],...
  'mask_midbrain_flag',false,[false true],...
  'outdir','analysis',[],...
  'aseg_flag',true,[false true],...
  'aseg_roigroups_flag',true,[false true],...
  'subhippo_flag',true,[false true],...
  'aparc_flag',true,[false true],...
  'nobiasflag',true,[false true],...
  'forceflag',false,[false true],...
... % undocumented:
  'baseline_outdir','analysis',[],... % in FSURF container
  'projdist',1,[-5,5],...
  'projfrac',[],[],...
  'hemilist',{'lh','rh'},{'lh' 'rh'},...
  'sphere_flag',true,[false true],...
  'required_containers',{'proc'},[],...
  'required_rootdirs',{'proc','fsurf','long'},[],...
  'modality','MRI',[],...  
... % not sure whether to keep
  'erode_flag',true,[false true],...
  'erode_nvoxels',1,[1:100],...
  'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79,10001:10003,20001:20003,20009,20010],[1,Inf],...
  'dirprefix','LONG',[],...
...
  'analyze_tags',{'baseflag' 'sphsmoothsteps' 'mask_midbrain_flag'...
                  'outdir' 'forceflag' 'projdist' 'projfrac'...
                  'hemilist' 'aseg_flag' 'aparc_flag' 'sphere_flag'...
                  'subhippo_flag' 'erode_flag' 'erode_nvoxels'...
                  'aseg_roilist' 'dirprefix' 'aseg_roigroups_flag'...
                  'nobiasflag'},[],...
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
end;
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
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

  % create script
  jstem = regexprep(SubjID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = [batchdir '/' jobID '.m'];
  matfname = [batchdir '/' jobID '.mat'];
  save(matfname,'RootDirs' ,'StudyInfo');
  fid = fopen(jobfname,'w');
  fprintf(fid,'load(''%s'')\n',matfname);
  fprintf(fid,'MMIL_Analyze_Long_Exam(RootDirs,StudyInfo,...\n');
  mmil_write_tags(fid,parms.analyze_tags,parms);
  fprintf(fid,');\n');
  fprintf(fid,'exit\n');
  fclose(fid);

  % add to list
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);

