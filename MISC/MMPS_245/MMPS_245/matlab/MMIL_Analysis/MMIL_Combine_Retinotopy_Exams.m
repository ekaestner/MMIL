function MMIL_Combine_Retinotopy_Exams(ProjID,varargin)
%function MMIL_Combine_Retinotopy_Exams(ProjID,[options])
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
%    including these fields: SubjID, VisitID, STRUCT_VisitID,
%      FP_snums, FP_revflags, FP_phase_offset,
%      FE_snums, FE_revflags, FE_phase_offset
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: proc_bold, fsurf (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied, MMIL_ProjInfo.csv is not required
%    {default = []}
%  'qcflag': [0|1] whether to exclude subjects with StudyInfo.QC=0
%    {default = 1}
%  'batchname': name of ougtput batchdir
%    {default = 'MMIL_Combine_Retinotopy_Exams'}
%  'multi_flag': [0|1|2] indicates source of input analyses
%    0: use results from individual scans
%    1: use results from multiple scans within session
%    2: use results from multiple sessions
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0 }
%
% Optional Parameters that specify how analysis was done
%  'infix': string attached to processed BOLD file names
%    {default = 'corr_resBOLD'}
%  'fstats_type': [0|1|2] how output Fourier components should be scaled
%    0: raw, no scaling
%    1: scaled by sqrt(F-ratio)
%    2: scaled by significance values (-log10(p))
%    {default = 2}
%  'resamp_flag': [0|1] whether Fourier results were resampled to 1x1x1
%     before painting
%    { default = 0 }
%  'smoothsteps': smoothing steps on surface after painting
%    {default = 0}
%
% Created:  04/10/10 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'infix',[],[],...
  'batchname','MMIL_Combine_Retinotopy_Exams',[],...
  'qcflag',true,[false true],...
  'multi_flag',0,[0:2],...
  'forceflag',false,[false true],...
...
  'fstats_type',2,[0:2],...
  'resamp_flag',false,[false true],...
  'cxfstatsflag',true,[false true],...
  'smoothsteps',0,[0,1000],...
... % parameters for viewing results
  'tksmooth',2,[0,1000],...
  'fthresh',0,[0,100],...
  'fmid',1.5,[0.1,100],...
  'fslope',1.0,[0.1,100],...
  'view','med',{'lat','med','ven','pos','dor'},...
  'surf','inflated',[],...
...
  'required_containers',{'proc_bold','fsurf'},[],...
  'modality','MRI',[],...
});

tags_excl = {'batchname','qcflag','required_containers','modality',...
  'StudyInfo','RootDirs'};
% fields to be passed to MMIL_Combine_Retinotopy_Exam
tags = setdiff(fieldnames(parms),tags_excl);

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});

if isempty(StudyInfo)
  fprintf('%s: ERROR: empty StudyInfo\n',mfilename);
  return;
end;

OrigStudyInfo = StudyInfo;
All_SubjIDs = {StudyInfo.SubjID};
SubjIDs = unique(All_SubjIDs);

% find pol and ecc sessions
ind_pol = find(~cellfun(@isempty,{StudyInfo.FP_snums}));
ind_ecc = find(~cellfun(@isempty,{StudyInfo.FE_snums}));
if isempty(ind_pol) && isempty(ind_ecc)
  fprintf('%s: WARNING: no sessions with retintotopy data\n',...
    mfilename);
  return;
end;

% create output batch directory
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
[s,msg] = mkdir(batchdir);
if ~s
  error('failed to create batch job directory %s',batchdir);
end;

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

j = 1;
for s=1:length(SubjIDs)
	SubjID = SubjIDs{s};
	ind = find(strcmp(SubjID,All_SubjIDs));
	StudyInfo = OrigStudyInfo(ind);

  % find pol and ecc sessions
  ind_pol = find(~cellfun(@isempty,{StudyInfo.FP_snums}));
  ind_ecc = find(~cellfun(@isempty,{StudyInfo.FE_snums}));
  if isempty(ind_pol) && isempty(ind_ecc)
    fprintf('%s: WARNING: no retintotopy data for %s\n',...
      mfilename,SubjID);
    continue;
  end;

  % create script
  jstem = regexprep(SubjID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = sprintf('%s/%s.m',batchdir,jobID);
  matfname = sprintf('%s/%s.mat',batchdir,jobID);
  save(matfname,'RootDirs' ,'StudyInfo');
  fid = fopen(jobfname,'w');
  fprintf(fid,'load(''%s'')\n',matfname);
  fprintf(fid,'MMIL_Combine_Retinotopy_Exam(RootDirs,StudyInfo,...\n');
  mmil_write_tags(fid,tags,parms);
  fprintf(fid,');\n');
  fprintf(fid,'exit\n');
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);
