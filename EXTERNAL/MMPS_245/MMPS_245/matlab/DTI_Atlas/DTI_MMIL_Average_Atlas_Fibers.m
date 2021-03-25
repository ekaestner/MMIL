function DTI_MMIL_Average_Atlas_Fibers(ProjID,varargin)
%function DTI_MMIL_Average_Atlas_Fibers(ProjID,varargin)
%
% Required Input:
%   ProjID: project ID to run (e.g. 'REC_TEST')
%
% Optional Parameters:
%  'StudyInfo': struct array of study information
%              (e.g. read from csv file with MMIL_Read_StudyInfo)
%               If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%    'proc_dti' and 'fibers'
%     specifying full paths of root directories containing data containers
%    {default = []}
%  'batchname': name of batch directory to contain job scripts
%    {default = [ProjID '_DTI_MMIL_Average_Atlas_Fibers']}
%  'batchrootdir': top level directory containing output batch job directories
%    {default = /home/$USER/batchdirs}
%  'required_containers': List of Containers which must exist for each study
%    {default = {'proc_dti'}}
%  'modality': only include subjects with particular modality
%    {default = 'MRI'}
%  'qcflag': only include subjects with StudyInfo.QC=1
%    {default = 1}
%  'fibers': vector of fiber numbers
%    {default = [101:110,115:123,133:138,141:150]}
%  'inclusion_modes': which subjects to include in the atlas, can be a cell array
%    valid values are those found in 'Group' column in <ProjID>_VisitInfo.csv
%    {default: 'all}
%
% Parameters for the function:
%  'outdir': path for storing the fiber atlas output
%    {default = [getenv('MMPS_DIR') '/atlases/DTI_Atlas/']}
%   tensor_smooth_sigma: smoothing applied to tensors after averaging
%     {default: 5}
%   countflag: [0|1] use fiber counts as input rather than fiber masks
%     {default: 1}
%   first_only_flag: [0|1] use first eigen vector only
%     {default = 1}
%   min_tensor_count: for calculating mean after smoothing
%     {default = 0.5}
%   smf: threshold used for the standardized tensors
%     {default = 10^-6}
%  'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Created:  03/8/11 by Vijay Venkatraman
% Rcnt Mod: 07/20/11 by Vijay Venkatraman
% Last Mod: 09/10/12 by Don Hagler
%

%% todo: age range instead of 'adult' and 'child'
%% todo: check inclusion modes are valid based on values in StudyInfo.Group

%  % Hidden parameters
%   fibers: List of fiber numbers 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return;end;
if isempty(ProjID), error('empty ProjID');end;

parms_filter = {....
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchname',[],[],...
  'batchrootdir',[],[],...
  'required_containers',{'proc_dti'},[],...
  'modality',{'MRI'},{'MRI','MEG','PET'},...
  'qcflag',true, [false true],...
  'fibers',[101:110,115:123,133:138,141:150],[],...
  'inclusion_modes','all',[],...
... % parameters for DTI_MMIL_Average_Atlas_Fiber   
  'outdir',[],[],...
  'tensor_smooth_sigma',5,[],... 
  'countflag',1,[],...
  'first_only_flag',1,[],...
  'min_tensor_count',0.5,[],...
  'smf',10^-6,[],...
  'forceflag',false,[true false],...
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,parms_filter);

% excl_tags are fields that should not be passed to DTI_MMIL_Average_Atlas_Fiber
excl_tags= {'StudyInfo','RootDirs','batchname','batchrootdir',...
    'required_containers','qcflag','modality','fibers','inclusion_modes'};
tags = setdiff(fieldnames(parms),excl_tags);

% get info from ProjInfo
args= MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

StudyInfo = StudyInfo(strcmp({StudyInfo.Modality},'MRI'));
if isempty(StudyInfo), error('empty StudyInfo'); end;
OrigStudyInfo = StudyInfo;
if ~isfield(OrigStudyInfo,'Group')
  OrigStudyInfo(1).Group = '';
end;
for s=1:length(OrigStudyInfo)
  if isempty(OrigStudyInfo(s).Group)
    OrigStudyInfo(s).Group = '';
  end;
end;
Groups = lower({OrigStudyInfo.Group});

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;
if isempty(parms.outdir)
  parms.outdir = [getenv('MMPS_DIR') '/atlases/DTI_Atlas'];   
end;

% check inclusion modes for validity
if ~iscell(parms.inclusion_modes)
  parms.inclusion_modes = {parms.inclusion_modes};
end;
for i=1:length(parms.inclusion_modes)
  if isempty(parms.inclusion_modes{i})
    parms.inclusion_modes{i} = 'all';
  end;
end;
valid_inclusion_modes = union(Groups,{'adult','child','all'});
ind_badmode = find(~ismember(parms.inclusion_modes,valid_inclusion_modes))
if ~isempty(ind_badmode)
  badmode_str = sprintf('''%s'' ',parms.inclusion_modes{ind_badmode});
  error('invalid exclusion modes: %s',badmode_str);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output batch directory
if isempty(parms.batchname), parms.batchname = mfilename; end;
if ~isempty(ProjID), parms.batchname = [ProjID '_' parms.batchname]; end;
parms.batchdir = [RootDirs.batch '/' parms.batchname];
if exist(parms.batchdir,'dir')
  cmd = sprintf('rm -rf %s\n',parms.batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(parms.batchdir);

% create scriptlist
scriptlistfname = sprintf('%s/scriptlist.txt',parms.batchdir);
fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create jobs for each fiber in StudyInfo
j = 1;
for i = 1:length(parms.inclusion_modes)
  inclusion_mode = parms.inclusion_modes{i};
  if isempty(inclusion_mode), inclusion_mode = 'all'; end;
  if ~strcmp(inclusion_mode,'all')
    switch inclusion_mode
     case 'adult'
      ind_keep = find(ismember(Groups,{'patient','control'}));
     case 'child'
      ind_keep = find(ismember(Groups,{'ace','ping'}));
     otherwise
      ind_keep = find(strcmp(Groups,parms.inclusion_mode));
    end;
    StudyInfo = OrigStudyInfo(ind_keep);
  else
    StudyInfo = OrigStudyInfo;
  end;
  nstudies = length(StudyInfo);
  if ~nstudies
    fprintf('%s: WARNING: no subjects for inclusion mode ''%s''\n',..
      mfilename,inclusion_mode);
    continue;
  end;
  tmp_parms = parms;
  tmp_parms.outdir = [parms.outdir '/' inclusion_mode];
  for f = parms.fibers
    % create script
    jobID = sprintf('job_%03d_%s_%03d',j,'fiber',i); j = j+1;
    jobfname = sprintf('%s/%s.m',tmp_parms.batchdir,jobID);
    matfname = sprintf('%s/%s.mat',tmp_parms.batchdir,jobID);
    save(matfname,'RootDirs' ,'StudyInfo');
    fid = fopen(jobfname,'w');
    fprintf(fid,'load(''%s'')\n',matfname);
    fprintf(fid,'DTI_MMIL_Average_Atlas_Fiber(RootDirs,StudyInfo,%d,...\n',f);
    mmil_write_tags(fid,tags,tmp_parms);
    fprintf(fid,');\n');
    fprintf(fid,'exit\n');
    fclose(fid);
    fid = fopen(scriptlistfname,'a');
    fprintf(fid,'%s\n',jobID);
    fclose(fid);
  end
end

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs3 %s\n',parms.batchname);

% check available disk space
MMIL_Check_Usage(RootDirs.outdir);


