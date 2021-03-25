function MMIL_Analyze_MRI_Exams(ProjID,varargin)
% function MMIL_Analyze_MRI_Exams(ProjID,[options])
%
% Usage:
%  MMIL_Analyze_MRI_Exams(ProjID,'key1', value1,...);
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
%  'batchrootdir': top level directory containing output batch job directories
%     {default = /home/$USER/batchdirs}
%  'batchname': name of output batchdir
%     {default = 'MMIL_Analyze_MRI_Exams'}
%
% Optional Parameters Controlling which Steps to Run:
%  'run_all_flag': run all steps, overriding other control flags
%     {default = 0}
%  'aseg_flag': read freesurfer's aseg.stats, save in mat file
%     {default = 1}
%  'aparc_flag': read freesurfer's aparc.stats, save in mat files
%     {default = 1}
%  'thick_flag': resample thickness to sphere, smooth
%     {default = 1}
%  'sulc_flag': resample sulc to sphere, smooth
%     {default = 1}
%  'area_flag': resample area to sphere (with Jacobian correction), smooth
%     {default = 1}
%  'cortvol_flag': calculate cortical volume from thicknes and area
%     {default = 1}
%  'T1w_aseg_flag': extract average T1-weighted values (nu.mgz) in aseg ROIs
%     {default = 1}
%  'T1w_surf_flag': sample T1-weighted volume to surface,
%       resample to sphere, smooth, compute gray-white contrast,
%       calculate aparc ROI averages
%     {default = 1}
%  'proc_aseg_flag': extract average values in aseg ROIs
%       with input files determined by 'proc_inputlist'
%     {default = 1}
%  'proc_surf_flag': sample volume data to surface,
%       resample to sphere, smooth, compute gray-white contrast,
%       calculate aparc ROI averages
%       with input files determined by 'proc_inputlist'
%     {default = 1}
%  'fuzzy_flag': [0|1] use fuzzy cluster weighted surface ROIs
%     found in atlases/fuzzy_clusters for thickness, area, and T1w
%     applies if one or more of thick_flag, area_flag,
%       cortvol_flag, or T1w_surf_flag = 1
%     {default = 1}
%  'aseg_roigroups_flag': [0|1] create masks for groups of aseg roi codes
%     includes 'WholeBrain', 'LatVentricles', and 'AllVentricles'
%     {default = 0}
%  'subhippo_flag': [0|1] create subdivided hippocampal ROIs
%     {default = 0}
%  'erode_flag': [0|1] create and use eroded ROIs (aseg, groups, hippo)
%     {default = 0}
%
% Optional Parameters:
%  'outdir': output directory relative to FreeSurfer Container
%     {default = 'analysis'}
%  'proc_inputlist': list of input file stems
%       relative to proc_indir e.g. 'MPR_res', 'T2w_res'
%       input files should be in register with orig.mgz
%         and have identical dimensions and vox2ras matrices
%       if empty, proc_aseg_flag and proc_surf_flag will be set to 0
%     {default = []}
%  'proc_outputlist': list of output file stems
%       corresponding to elements of proc_inputlist
%       if empty, will use file stems from proc_inputlist
%         e.g. {'T1w','T2w'}
%       if proc_outputlist contains 'T1w'
%         T1w_aseg_flag and T1w_surf_flag will be set to 0
%     {default = []}
%  'proc_scalefacts': vector of scaling factors for
%       each element of proc_inputlist
%       if empty, will be set to 1 for each
%     {default = []}
%  'aparc_infix': string in aparc annot and stats file
%     e.g. 'aparc', 'aparc.a2009s'
%    {default = 'aparc'}
%  'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     {default = [-0.2,0.2]}
%  'sphsmoothsteps': surface smoothing steps on ico sphere
%     may be vector of values to loop over
%     slope of FWHM vs. sqrt(N) is ~1.13 for fsaverage (v3)
%     (FWHM = full-width-half-max smoothing kernel
%         N = number of smoothing steps)
%      with [176,705,2819], approx FWHM (mm) = 15,30,60
%      with [100,300,3000], approx FWHM (mm) = 11.3, 19.6, 61.9
%     {default = 2819}
%  'mask_midbrain_flag': [0|1] whether to mask out mid brain and other
%     cortical regions marked "unknown" (masking done before smoothing)
%     {default = 0}
%  'fuzzy_dir': input directory for fuzzy cluster ROIs
%     if empty, will use [getenv('MMPS_dir') '/atlases/fuzzy_clusters']
%     {default = []}
%  'fuzzy_fstem': input file stem for fuzzy cluster ROIs
%     with expected names like {fstem}{order}-{hemi}.mgz
%     {default = 'fuzzy'}
%  'fuzzy_order': [0|2|4|12|18] number of fuzzy cluster ROIs
%     note: set of 18 includes combined sets of 2, 4, and 12
%     if order=0, names are like {fstem}-{hemi}.mgz
%     {default = 18}
%  'fuzzy_thresh': threshold applied to fuzzy cluster ROIs
%     {default = 0}
%  'fuzzy_smooth': smoothing applied to surface maps before extracting
%     values from fuzzy clusters
%     {default = 2819}
%  'qcflag': [0|1] whether to only include studies with StudyInfo.QC=1
%     {default = 1}
%  'check_stale_flag': check creation date of output directory
%     and delete if stale (created prior to fs.finish.all.touch)
%     {default = 0}
%  'check_complete_flag': [0|1] whether to require that recon is complete
%    {default = 1}
%  'FS_version': which version of Freesurfer used (e.g. 305, 450, 510, 530)
%    for checking whether recon is complete
%    if empty, will use FREESURER_VER environment variable
%      or get from ContainerInfo
%    {default = []}
%  'verbose': [0|1] display status meassages
%     {default = 0}
%  'forceflag': overwrite existing output
%     {default = 0}
%
% Created:  03/31/09 by Don Hagler
% Prev Mod: 07/31/16 by Don Hagler
% Last Mod: 07/11/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchrootdir',[],[],...
  'batchname','MMIL_Analyze_MRI_Exams',[],...
... % control flags:
  'run_all_flag',false,[false true],...
  'aseg_flag',true,[false true],...
  'aparc_flag',true,[false true],...
  'thick_flag',true,[false true],...
  'sulc_flag',true,[false true],...
  'area_flag',true,[false true],...
  'cortvol_flag',true,[false true],...
  'T1w_aseg_flag',true,[false true],...
  'T1w_surf_flag',true,[false true],...
  'proc_aseg_flag',true,[false true],...
  'proc_surf_flag',true,[false true],...
  'fuzzy_flag',true,[false true],...
  'aseg_roigroups_flag',false,[false true],...
  'subhippo_flag',false,[false true],...
  'erode_flag',false,[false true],...
...
  'outdir','analysis',[],...
  'proc_inputlist',[],[],...
  'proc_outputlist',[],[],...
  'proc_scalefacts',[],[],...
  'aparc_infix','aparc',[],...
  'projdist_list',[-0.2,0.2],[-5,5],...
  'sphsmoothsteps',2819,[0,Inf],...
  'mask_midbrain_flag',false,[false true],...
  'fuzzy_dir',[],[],...
  'fuzzy_fstem','fuzzy',[],...
  'fuzzy_order',18,[0,2,4,12,18],...
  'fuzzy_thresh',0,[0,Inf],...
  'fuzzy_smooth',2819,[0,Inf],...
  'qcflag',true,[false true],...
  'check_stale_flag',false,[false true],...
  'check_complete_flag',true,[false true],...
  'FS_version',[],[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
... % undocumented:
  'sphere_flag',true,[false true],...
  'erode_nvoxels',1,[1:100],...
  'hemilist',{'lh','rh'},{'lh' 'rh'},...
  'aseg_roigroups',[],[],...
  'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79,251:255],[1,Inf],...
  'aseg_aparc_flag',0,[0,1,2],...
  'required_rootdirs',{'fsurf','fsico'},[],...
  'QC_raw',true,[false true],...
  'QC_recon',true,[false true],...
});

excl_tags = {'StudyInfo','RootDirs','batchrootdir','batchname',...
  'qcflag','required_rootdirs','QC_raw','QC_recon'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check proc_inputlist
if isempty(parms.proc_inputlist)
  parms.proc_aseg_flag = 0;
  parms.proc_surf_flag = 0;
end;
if parms.proc_aseg_flag || parms.proc_surf_flag
  ninputs = length(parms.proc_inputlist);
  if ~iscell(parms.proc_inputlist)
    parms.proc_inputlist = {parms.proc_inputlist};
  end;
  % check proc_outputlist
  if isempty(parms.proc_outputlist)
    parms.proc_outputlist = parms.proc_inputlist;
  elseif ~iscell(parms.proc_outputlist)
    parms.proc_outputlist = {parms.proc_outputlist};
  end;
  if length(parms.proc_outputlist) ~= ninputs
    error('number of elements in proc_outputlist does not match proc_inputlist');
  end;
  if isempty(parms.proc_scalefacts)
    parms.proc_scalefacts = ones(ninputs,1);
  elseif length(parms.proc_scalefacts) ~= ninputs
    error('number of elements in proc_scalefacts does not match proc_inputlist');
  end;
  % disable T1w analysis if proc_outputlist contains 'T1w'
  if ismember('T1w',parms.proc_outputlist)
    parms.T1w_aseg_flag = 0;
    parms.T1w_surf_flag = 0;
  end;
end;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;

% set FreeSurfer version if not already set
if isempty(parms.FS_version)
  parms.FS_version = str2num(getenv('FREESURFER_VER'));
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
fpaths = {};
for s=1:length(StudyInfo)
  % select study info for this subject visit
  SubjID = StudyInfo(s).SubjID;
  VisitID = StudyInfo(s).VisitID;
  tparms = parms;

  % skip if no processed container
  if parms.proc_aseg_flag || parms.proc_surf_flag
    tparms.proc_indir = [RootDirs.proc '/' StudyInfo(s).proc];
    if isempty(StudyInfo(s).proc) | ~exist(tparms.proc_indir,'dir')
      fprintf('%s: skipping %s because missing proc container\n',...
        mfilename,VisitID);
      continue;
    end;
  end;

  % skip if no FreeSurfer recon container
  FSContainerPath = [RootDirs.fsurf '/' StudyInfo(s).fsurf];
  if isempty(StudyInfo(s).fsurf) | ~exist(FSContainerPath,'dir')
    fprintf('%s: skipping %s because missing fsurf container\n',...
      mfilename,VisitID);
    continue;
  end;

  % skip job for duplicate container path (STRUCT_VisitID)
  if ismember(FSContainerPath,fpaths), continue; end;
  fpaths{end+1} = FSContainerPath;

  % create script
  jstem = regexprep(SubjID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = [batchdir '/' jobID '.m'];
  tags = setdiff(fieldnames(tparms),excl_tags);
  mmil_write_script(jobfname,'MMIL_Analyze_MRI_Exam',...
    {FSContainerPath},tags,tparms);

  % add to list
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

% check available disk space
MMIL_Check_Usage(RootDirs.fsurf);

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);

