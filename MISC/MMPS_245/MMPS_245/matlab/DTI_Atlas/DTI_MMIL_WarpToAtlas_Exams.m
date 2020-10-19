function DTI_MMIL_WarpToAtlas_Exams(ProjID,varargin)
%function DTI_MMIL_WarpToAtlas_Exams(ProjID,varargin)
%
% Usage:
%  DTI_MMIL_WarpToAtlas_Exams(ProjID,'key1', value1,...);
%
% Required Parameters:
%   ProjID: Project ID string
%           used to to load ProjInfo and StudyInfo from user's home
%           (e.g. '/home/user/ProjInfo/REC_MMIL_ProjInfo.csv'
%             '/home/user/ProjInfo/ProjID/ProjID_VisitInfo.csv' )
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information (see MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%     'proc','proc_dti'
%     specifying full paths of root directories containing data containers
%     {default = []}
%  'batchrootdir': top level directory containing output batch job directories
%     {default = /home/$USER/batchdirs}
%  'batchname': name of output batchdir
%     {default = 'DTI_MMIL_WarpToAtlas'}
%  'required_rootdirs': List of rootdirs fields which must not be empty
%     {default= {'proc','proc_dti'}}
%  'qcflag': only include subjects with StudyInfo.QC=1
%     {default = 1}
%
%  Optional Parameters:
%  'outdir': output directory relative to ContainerPath
%     {default = 'WarpToAtlas'}
%  'STRUCT_T1type': which T1 series to reg the DTI to
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default = 2}
%  'DTI_snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%  'DTI_min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'DTI_flex_flag': [0|1] DTI_flex scans included in tensor fit
%    {default = 0}
%  'DTI_nob0_flag': [0|1] whether to exclude b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'DTI_min_ndirs': minimum number of directions to be included in tensor fit 
%     {default = 6}
%  'DTI_min_nb0': minimum number of B0 images
%     {default = 1}
%  'DTI_fiber_infix': infix of fibers used to warp to atlas
%     {default = 'corr_regT1'}
%  'DTI_revflag': reverse flag option 
%     {default = 0}
%  'DTI_fibers': fiber numbers to be warped to atlas
%     {default=[101:110,115:123,133:138,141:150]}
%  'DTI_measlist': DTI measures to be warped to atlas
%     {default = {FA, b0}}
%  'fiber_countflag': [0|1] whether fiber was imported as count map or mask map
%     {default = 1}
%  'sparse_flag': [0|1] save volumes as sparse mat files
%     {default = 1}
%
% Optional Parameters controlling processing steps:
%  'warp_masks_flag': [0|1] whether to warp brain masks to atlas
%     {default = 1}
%  'warp_tensors_flag': [0|1] whether to warp tensors, eigen vectors,
%     and eigen values to atlas
%     {default = 1}
%  'warp_DTmeas_flag': [0|1] whether to warp DTI measures (FA, etc.) to atlas
%     {default = 1}
%  'warp_fibers_flag': [0|1] whether to warp fibers to atlas
%     {default = 1}
%  'warp_fiber_tensors_flag': [0|1] whether to smooth tensors within fibers
%     and warp to atlas
%     {default = 1}
%  'first_only_flag': [0|1] whether to remove 2nd and 3rd Eigenvectors
%     {default = 1}
%  'atlasfibers_flag': [0|1] whether to use original, manual fibers (0)
%     or new fibers from atlas (1)
%     {default = 0}
%  'forceflag': overwrite existing output
%     {default = 0}
%
% Created:  03/03/11 by Vijay Venkatraman
% Prev Mod: 08/09/13 by Don Hagler
% Last Mod: 06/07/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if isempty(ProjID), error('empty ProjID');end;

parms_all = { ...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchrootdir',[],[],...
  'batchname',[],[],...
  'required_rootdirs',{'proc','proc_dti'},[],...
  'qcflag',true, [false true],...
  'modality',{'MRI'},{'MRI','MEG','PET'},...
... % parameters for selecting the data
  'DTI_min_ndirs',6,[],...
  'DTI_min_bval',1,[],...
  'DTI_flex_flag',false,[false true],...
  'DTI_min_nb0',1,[],...
  'DTI_snums',[],[],...
  'DTI_revflag',0,[0 1 2],...
  'DTI_fiber_infix','corr_regT1',[],...
  'DTI_nob0_flag',false,[false true],...
  'STRUCT_T1type',2,[],... 
  'T1ContainerPath',[],[],...
  'FSContainerPath',[],[],...
  'DTI_measlist',{'FA','b0'},[],...
  'DTI_fibers',[101:110,115:123,133:138,141:150],[],...
  'fiber_countflag',true,[false true],...
  'sparse_flag',true,[false true],...
  'outdir','WarpToAtlas',[],...
... % flags
  'warp_masks_flag',true,[false true],...
  'warp_tensors_flag',true,[false true],...
  'warp_DTmeas_flag',true,[false true],...
  'warp_fibers_flag',true,[false true],...
  'warp_fiber_tensors_flag',true,[false true],...
  'first_only_flag',true,[false true],...
  'atlasfibers_flag',false,[false true],...
...
  'fiber_maps_infix','sm5.00_probV0_countatlas',[],... %% todo: should be constructed from other input parms
  'forceflag',false,[false true],...
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parms = mmil_args2parms(varargin,parms_all);

% excl_tags are fields that should not be passed to DTI_MMIL_WarpToAtlas_Exam
excl_tags= {'StudyInfo','RootDirs','batchname','batchrootdir',...
    'required_rootdirs','qcflag','modality'};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  if isempty(StudyInfo(i).proc_dti)
    fprintf('%s: WARNING: no processed DTI container for %s ... skipping\n',...
      mfilename,StudyInfo(i).VisitID);
    continue;
  end;
  tmp_parms = parms;
  tmp_parms.DTI_snums = StudyInfo(i).DTIScanNums; 
  ContainerPath = [RootDirs.proc_dti '/' StudyInfo(i).proc_dti];
  tmp_parms.T1ContainerPath = [RootDirs.proc '/' StudyInfo(i).proc];
  if isfield(RootDirs,'fsurf') && ~isempty(StudyInfo(i).fsurf)
    tmp_parms.FSContainerPath = [RootDirs.fsurf '/' StudyInfo(i).fsurf];
  end;
  if ~exist(ContainerPath,'dir') 
    fprintf('%s: WARNING: ContainerPath %s not found ... skipping\n',...
      mfilename,ContainerPath);
    continue;
  end;
  if ~exist(tmp_parms.T1ContainerPath,'dir') 
    fprintf('%s: WARNING: T1ContainerPath %s not found ... skipping\n',...
      mfilename,tmp_parms.T1ContainerPath);
    continue;
  end;
  if ~isempty(tmp_parms.FSContainerPath) &&...
     ~exist(tmp_parms.FSContainerPath,'dir') 
    fprintf('%s: WARNING: FSContainerPath %s not found\n',...
      mfilename,tmp_parms.FSContainerPath);
    tmp_parms.FSContainerPath = [];
  end;
  % create script
  jobID = sprintf('job_%03d_%s',j,StudyInfo(i).VisitID); j = j+1;
  jobfname = [batchdir '/' jobID '.m'];
  mmil_write_script(jobfname,'DTI_MMIL_WarpToAtlas_Exam',...
    ContainerPath,tags,tmp_parms);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs3 %s\n',parms.batchname);

