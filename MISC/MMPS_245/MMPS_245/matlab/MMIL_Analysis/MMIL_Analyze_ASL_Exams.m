function MMIL_Analyze_ASL_Exams(ProjID,varargin)
%function MMIL_Analyze_ASL_Exams(ProjID,[options])
%
% Usage:
%  MMIL_Analyze_ASL_Exams(ProjID,'key1', value1,...);
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
%      fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'batchrootdir': top level directory containing output batch job directories
%     {default = /home/$USER/batchdirs}
%  'batchname': name of output batchdir
%     {default = 'MMIL_Analyze_ASL_Exams'}
%
% Optional Parameters:
%   'outdir': analysis output directory
%     full path or relative to OutPath
%     {default = 'analysis'}
%   'aseg_flag': [0|1] extract aseg ROI results
%     {default = 1}
%   'cortsurf_flag': [0|1] paint to cortical surface
%     and extract aparc ROI results
%     {default = 1}
%   'measlist': list of ASL "measures" to extract for each ROI
%      e.g. 'CBF', cerebral blood flow
%     {default = {'CBF'}}
%   'scalefacts': scaling factors applied to each measure in measlist
%     {default = [1]}
%   'minval': minimum value; if NaN, include all voxels
%     {default = 1e-6}
%   'csv_flag': [0|1] output results in csv files
%     {default = 1}
%   'verbose': [0|1] display status meassages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Optional Parameters specific to aseg ROIs:
%   'erode_flag': [0|1] whether to "erode" ROIs
%     by smoothing and thresholding (to reduce edge effects)
%     {default = 1}
%   'erode_nvoxels': number of voxels to erode (integer)
%     {default = 1}
%   'aseg_aparc_flag': [0|1|2] whether to use cortical parcellation ROIs
%     0: aseg only
%     1: aparc only   NOTE: for best results with aparc, use erode_flag=0
%     2: aparc+aseg
%     {default = 0}
%
% Optional Parameters specific to cortical surface:
%   'fnames_aparc': cell array of annotation files (one for each hemisphere)
%     if empty, will use ?h.aparc.annot files in fspath/label
%     if not full path, assumed to be relative to fspath/label
%     {default = []}
%   'fnames_fparc': cell array of annotation files in fsaverage space
%     will be resampled to individual subject space before use
%     {default = []}
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast maps
%     {default = 1}
%   'gwnorm_flag': [0|1] for gray/white contrast, normalize difference by mean
%     {default = 0}
%   'smoothsteps': number of smoothing iterations on individual subject surface
%     {default = 0}
%   'sphere_flag': [0|1[ whether to resample to spherical atlas
%     {default = 0}
%   'sphsmoothsteps': number of smoothing iterations on sphere
%     {default = 0}
%
% Optional Parameters for resampling volumes to atlas space
%   'atlas_warp_flag': [0|1] nonlinearly resample input volumes to atlas space
%     {default = 0}
%   'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%   'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%
% Created:  02/07/13 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchrootdir',[],[],...
  'batchname','MMIL_Analyze_ASL_Exams',[],...
...
  'outdir','analysis',[],...
  'aseg_flag',true,[false true],...
  'cortsurf_flag',false,[false true],...
  'measlist',{'CBF'},[],...
  'scalefacts',[1],[],...
  'csv_flag',true,[false true],...
  'minval',1e-6,[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
... % aseg parameters
  'erode_flag',true,[false true],...
  'erode_nvoxels',1,[1:100],...
  'aseg_aparc_flag',0,[0,1,2],...
... % cortsurf parameters
  'fnames_aparc',[],[],...
  'fnames_fparc',[],[],...
  'projdist_list',[1],[-5,5],...
  'gwnorm_flag',false,[false true],...
  'smoothsteps',0,[0,Inf],...
  'sphere_flag',false,[false true],...
  'sphsmoothsteps',0,[0,Inf],...
... % resample atlas parameters
  'atlas_warp_flag',false,[false true],...
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],...    
... % hidden parameters
  'input_files',{'intermediate/anat+reg+rs+orig.BRIK','CBF+orig.BRIK'},[],...
  'output_files',{'T1_lowres.mgz','CBF.mgz'},[],...
  'resT1_outfix','resT1',[],...
  'fnum',1,[1,Inf],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'fname_colorlut',[],[],...
  'required_rootdirs',{'fsurf','raw_asl','proc_asl'},[],...
  'QC_recon',true,[false true],...
};
parms = mmil_args2parms(varargin,parms_filter);

excl_tags = {'StudyInfo','RootDirs','batchrootdir','batchname',...
  'required_rootdirs','QC_recon'};

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;

if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% create jobs for each subject in StudyInfo
j = 1;
for i=1:length(StudyInfo)
  VisitID = StudyInfo(i).VisitID;
  indir = StudyInfo(i).raw_asl;
  if isempty(indir)
    if parms.verbose
      fprintf('%s: WARNING: missing raw_asl dir for %s\n',...
        mfilename,VisitID);
    end;
    continue;
  end;
  InPath = sprintf('%s/%s',RootDirs.raw_asl,indir);
  if ~exist(InPath,'dir')
    if parms.verbose
      fprintf('%s: WARNING: missing raw_asl dir %s not found\n',...
        mfilename,InPath);
    end;
    continue;
  end;
  fsdir = StudyInfo(i).fsurf;
  if isempty(fsdir)
    if parms.verbose
      fprintf('%s: WARNING: missing fsurf dir for %s\n',...
        mfilename,VisitID);
    end;
    continue;
  end;
  FSContainerPath = sprintf('%s/%s',RootDirs.fsurf,fsdir);
  if ~exist(FSContainerPath,'dir')
    if parms.verbose
      fprintf('%s: WARNING: missing fsurf dir %s not found\n',...
        mfilename,FSContainerPath);
    end;
    continue;
  end;
  OutPath = sprintf('%s/ASLPROC_%s',RootDirs.proc_asl,VisitID);
  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
  jobfname = sprintf('%s/%s.m',batchdir,jobID);
  tags = setdiff(fieldnames(parms),excl_tags);
  mmil_write_script(jobfname,'MMIL_Analyze_ASL_Exam',...
    {InPath,OutPath,FSContainerPath},tags,parms);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end

% check available disk space
MMIL_Check_Usage(RootDirs.proc_asl);

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);


