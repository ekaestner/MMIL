function DTI_MMIL_Import_DTIStudio_Fibers_Exams(ProjID,varargin)
%function DTI_MMIL_Import_DTIStudio_Fibers_Exams(ProjID,[options])
%
% Required input: 
%   ProjID: Project ID string
%           used to to load ProjInfo and StudyInfo from user's home
%           (e.g. '/home/user/ProjInfo/REC_MMIL_ProjInfo.csv'
%             '/home/user/ProjInfo/ProjID/ProjID_VisitInfo.csv' )
%
% Optional Parameters:
%  'RootDirs': struct that must contain the following fields:
%     'proc_dti' and 'fibers'
%     specifying full paths of root directories containing data containers
%     {default = []}
%  'StudyInfo': struct array of study information 
%              (e.g. read from csv file with MMIL_Read_StudyInfo)
%              If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'batchrootdir': top level directory containing output batch job directories
%    {default = /home/$USER/batchdirs}
%  'batchname': name of output batchdir
%    {default = 'DTI_MMIL_Import_DTIStudio_Fibers_Exams'}
%  'required_rootdirs': List of rootdirs fields which must not be empty
%    {default= {'proc_dti','fibers'}}
%  'qcflag': only include subjects with StudyInfo.QC=1
%    {default=1}
%	 'forceflag': [0|1] overwrite existing output
%    {default=0}
%
% Optional Parameters that determine input diffusion data:
%  'DTI_min_ndirs': minimum number of directions for ContainerPath
%    {default = 6}
%	 'DTI_min_bval': minimum b-value for ContainerPath
%    {default = 1000}
%  'DTI_flex_flag': [0|1] DTI_flex scans included in tensor fit
%    {default = 0}
%	 'DTI_min_nb0': minimum number of b=0 images for ContainerPath
%    {default = 1}
%	 'DTI_snums': scan numbers to be included for ContainerPath
%    {default = []}
%	 'DTI_revflag': whether to include rev phase-encode scans for ContainerPath
%    {default = 0}
%	 'DTI_fiber_infix': DTI processed data infix for ContainerPath
%    {default = 'corr_regT1'}
% 
% Optional Parameters that determine manual fiber data
%  'Fiber_studytype': method to create fiber file path
%    0: use snums and infix
%    1: use min_bval and infix
%    {default = 0}
%  'Fiber_dirname': Name of manual fiber directory
%    {default = 'fibers'}
%  'Fiber_snums': scan numbers to be included for FiberContainerPath
%    {default = []}
%  'Fiber_infix': Fiber ContainerPath infix
%    {default = 'reg_mc_ecc_B0uw_gruw_iso'}
%  'Fiber_min_ndirs': minimum number of directions for FiberContainerPath
%    {default = 6}
%  'Fiber_min_bval': minimum b-value for FiberContainerPath
%    {default = 1000}
%  'Fiber_flex_flag': [0|1] DTI_flex included for FiberContainerPath
%    {default = 0}
%  'Fiber_min_nb0': minimum number of b=0 images for FiberContainerPath
%    {default = 1}
%	 'Fiber_DTIStudio_fpat': manual fibers filename pattern
%    {default = 'Fiber_(?<fnum>\d+)_path.dat'}
%	 'Fiber_revsliceflag': [0|1] reverse slice order to LPS  
%    {default = 1}
%	 'Fiber_permvec': permutation order in case the slices are not axial
%    {default = [1,2,3]}
%	 'Fiber_mask_datatype': datatype for mask output
%    {default = 'uint8'}
%	 'Fiber_count_datatype' datatype for count output
%    {default = 'float'}
%	 'Fiber_countflag': [0|1] output option as mask or count maps
%    {default = 1}
%
% Created:  02/22/11 by Vijay Venkatraman (Original Code from Don Hagler)
% Last Mod: 02/04/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if isempty(ProjID), error('empty ProjID');end;

parms_all = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchrootdir',[],[],...
  'batchname', [],[],...
  ...
  'required_rootdirs',{'proc_dti','fiberdir'},[],...
  'qcflag',true, [false true],...
  ...
  'forceflag',false,[false true],...
	... % DTI ContainerPath parameters
	'DTI_min_ndirs',6,[],...
	'DTI_min_bval',1000,[],...
  'DTI_flex_flag',false,[false true],...
	'DTI_min_nb0',1,[],...
	'DTI_snums',[],[],...
	'DTI_revflag',2,[0 1 2],...
	'DTI_fiber_infix','corr_regT1',[],...
	... % Fiber Container
  'Fiber_studytype',0,[0 1]...
  'Fiber_dirname','fibers',[],...  
  'Fiber_snums',[],[],...
  'Fiber_infix','reg_mc_ecc_B0uw_gruw_iso',[],...
  'Fiber_min_ndirs',6,[],...
  'Fiber_min_bval',1000,[],...
  'Fiber_flex_flag',false,[false true],...
  'Fiber_min_nb0',1,[],...
	'Fiber_DTIStudio_fpat','Fiber_(?<fnum>\d+)_path.dat',[],...
	'Fiber_revsliceflag',1,[],...
	'Fiber_permvec',[1,2,3],[],...
	'Fiber_mask_datatype','uint8',[],...
	'Fiber_count_datatype','float',[],...
  'Fiber_countflag',1,[0 1],...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,parms_all);

% excl_tags are fields that should not be passed to DTI_MMIL_Import_DTIStudio_Fibers_Exam
excl_tags= {'StudyInfo','RootDirs','batchname',...
  'batchrootdir','required_rootdirs','qcflag'};
tags = setdiff(fieldnames(parms),excl_tags);

% Get input from ProjInfo
args= MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_all);
end;

StudyInfo = StudyInfo(strcmp({StudyInfo.Modality},'MRI'));
if isempty(StudyInfo)
  fprintf('%s: ERROR: empty StudyInfo\n',mfilename);
  return;
end;

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output batch directory
if isempty(parms.batchname), parms.batchname = mfilename; end;
if ~isempty(ProjID), parms.batchname = [ProjID '_' parms.batchname]; end;
batchdir = [RootDirs.batch '/' parms.batchname];
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

% create scriptlist
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

% create jobs for each subject in StudyInfo
j = 1;
for i=1:length(StudyInfo)
  tmp_parms = parms;
  ContainerPath = [RootDirs.proc_dti '/' StudyInfo(i).proc_dti];
  if ~exist(ContainerPath,'dir') 
    fprintf('%s: WARNING: ContainerPath %s not found ... skipping\n',...
      mfilename,ContainerPath);
    continue;
  end;
  FiberContainerPath = [RootDirs.fiberdir '/' StudyInfo(i).proc_dti];
  if ~exist(FiberContainerPath,'dir') 
    % possibly has old style of MRIPROCESSED container type
    oldFiberContainerPath = ...
      regexprep(FiberContainerPath,'DTIPROC','MRIPROCESSED');
    if exist(oldFiberContainerPath,'dir')
      FiberContainerPath = oldFiberContainerPath;
    else
      fprintf('%s: WARNING: FiberContainerPath %s not found ... skipping\n',...
        mfilename,FiberContainerPath);
      continue;
    end;
  end;
  tmp_parms.DTI_snums = StudyInfo(ind).DTIScanNums; 

  % create script
  jobID = sprintf('job_%03d_%s',j,VisitID); j = j+1;
  jobfname = [batchdir '/' jobID '.m'];
  mmil_write_script(jobfname,'MMIL_RetFit_BOLD_Exam',...
    {ContainerPath,FiberContainerPath},tags,tmp_parms);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);


