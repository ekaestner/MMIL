function MMIL_Analyze_BOLD_Exams(ProjID,varargin)
% function MMIL_Analyze_BOLD_Exams(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject including fields:
%    SubjID, VisitID, STRUCT_VisitID, FP_snums, and FE_snums
%    May include several other FP, FE, and GLM options
%    (see MMIL_Analyze_BOLD_Exam)
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: proc_bold, fsurf (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied, MMIL_ProjInfo.csv is not required
%    {default = []}
%  'qcflag': [0|1] whether to exclude subjects with StudyInfo.QC=0
%    {default = 1}
%  'batchname': name of output batchdir
%    {default = 'MMIL_Analyze_BOLD_Exams'}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional Parameters controlling which analysis steps are run:
%  'FP_flag': [0|1] whether to run polar Fourier analyses
%    for those subjects with FP_snums in StudyInfo
%    {default = 1}
%  'FE_flag': [0|1] whether to run eccen Fourier analyses
%    for those subjects with FE_snums in StudyInfo
%    {default = 1}
%  'GLM_flag': [0|1] whether to run GLM analyses
%    for those subjects with GLM_snums in StudyInfo
%    {default = 1}
%  'GLM_ROI_flag': [0|1] whether to run GLM ROI analysis
%    {default = 1}
%
% Optional Parameters specific to Fourier analysis:
%  'fstats_type': [0|1|2] how output Fourier components should be scaled
%    0: raw, no scaling
%    1: scaled by sqrt(F-ratio)
%    2: scaled by significance values (-log10(p))
%    {default = 2}
%  'cxfstatsflag' [0|1] whether to calculate cross-scan complex f-stats
%    (in addition to amplitude f-stats)
%    {default = 1}
%
% Optional Parameters specific to GLM analysis:
%  'GLM_concat_flag': [0|1] concatenate across multiple runs
%     Each scan will also be analyzed individually
%    {default = 1}
%  'GLM_pthresh': probability threshold applied to GLM f-stats for each
%     condition
%    {default = 0.05}
%  'GLM_contrasts_flag': [0|1] calculate glt contrasts between each condition
%    {default = 0}
%  'GLM_iresp_flag': [0|1] output impulse response functions for each condition
%    {default = 0}
%
% Optional Parameters specific to GLM ROI analysis:
%  'GLM_ROI_outstem': string added to output file names (appended to GLM stem)
%    {default = []}
%  'GLM_ROI_dir': name of directory containing label files relative to
%    FreeSurfer ContainerPath, or, if StudyInfo has GLM_ROI_VisitID for a
%    subject, relative to processed BOLD ContainerPath
%    {default = 'label'}
%
% Optional Parameters applicable to all analyses:
%  'infix': input file will be sprintf('BOLD%d_%s.mgz',snum,infix)
%    to specify no infix (i.e. raw data), use 'none'
%    {default = 'corr_resBOLD'}
%  'mc_flag': [0|1] whether within-scan motion correction was done
%    Will use motion.1D files as nuisance regressors
%    {default = 1}
%  'mc_inter_flag': [0|1] whether between-scan motion correction was done
%    Allows for a single reference scan to be used for registration to FS
%    {default = 1}
%  'regFS_flag': [0|1] whether to register BOLD scans to FreeSurfer nu.mgz
%    if registration already done, will not redo
%    necessary for painting to surface (if 0, paint_flag is ignored)
%    {default = 1}
%  'paint_flag': [0|1|2] whether to paint results to surface
%     if 0, only do Fourier or GLM analysis, no painting to surface
%     if 1, paint results to surface after averaging across scans
%       (requires motion correction and register.dat file for reference scan)
%     if 2, paint results to surface before averaging across scans
%       (requires register.dat file for each scan if mc_inter_flag=0)
%    {default = 1}
%  'projdist': dist (mm) to project along normal when painting
%    {default = 1}
%  'projfrac': fractional dist to project along normal when painting
%    {default = 0.5}
%  'projfract_flag': [0|1] whether to use projfrac (1) or projdist (0)
%    {default = 0}
%  'projfrac_avg': vector of [min max del] for averaging multiple samples
%    If empty, use projfrac instead if projfrac_flag=1
%    {default = []}
%  'projdist_avg': vector of [min max del] for averaging multiple samples
%    If empty, use projdist instead
%    {default = []}
%  'resamp_flag': [0|1] whether to resample results to 1x1x1 before painting
%    {default = 0}
%  'force_repaint_flag': [0|1] whether to repaint even if output files exist
%    (do not redo volume calculations)
%    {default = 0}
%  'sphere_flag': [0|1] whether to sample to icosohedral sphere
%    {default = 0}
%
% NOTE: project-specific defaults for any of the
%  'FP', 'FE', or 'GLM' options may be set in MMIL_ProjInfo.csv file
%
% Created:  06/04/09 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchname','MMIL_Analyze_BOLD_Exams',[],...
  'qcflag',true,[false true],...
  'forceflag',false,[false true],...
...
  'FP_flag',true,[false true],...
  'FE_flag',true,[false true],...
  'GLM_flag',true,[false true],...
  'GLM_ROI_flag',true,[false true],...
...
  'FP_datatype','polar',[],...
  'FP_snums',[],[],...
  'FP_revflags',[],[],...
  'FP_skipTRs',0,[0,1000],...
  'FP_phase_offset',0,[-1,1],...
  'FP_phase_offset_postrev',0,[-1,1],...
  'FP_stimfreq',8,[1,1000],...
  'FP_tksmooth',2,[0,1000],...
  'FP_smoothsteps',0,[0,1000],...
  'FP_sphsmoothsteps',0,[0,1000],...
...
  'FE_datatype','eccen',[],...
  'FE_snums',[],[],...
  'FE_revflags',[],[],...
  'FE_skipTRs',0,[0,1000],...
  'FE_phase_offset',0,[-1,1],...
  'FE_phase_offset_postrev',0,[-1,1],...
  'FE_stimfreq',8,[1,1000],...
  'FE_tksmooth',2,[0,1000],...
  'FE_smoothsteps',0,[0,1000],...
  'FE_sphsmoothsteps',0,[0,1000],...
...
  'GLM_snums',[],[],...
  'GLM_concat_flag',true,[false true],...
  'GLM_contrasts_flag',false,[false true],...
  'GLM_iresp_flag',false,[false true],...
  'GLM_skipTRs',0,[0,1000],...
  'GLM_minlag',0,[0,10],...
  'GLM_maxlag',4,[0,30],...
  'GLM_tksmooth',2,[0,1000],...
  'GLM_stim_fnames',[],[],...
  'GLM_num_conds',[],[],...
  'GLM_smoothsteps',0,[0,1000],...
  'GLM_sphsmoothsteps',0,[0,1000],...
  'GLM_pthresh',0.05,[0,1],...
...
  'GLM_ROI_VisitID',[],[],...
  'GLM_ROI_outstem',[],[],...
  'GLM_ROI_dir','label',[],...
  'GLM_ROI_stem',[],[],...
  'GLM_ROI_retfit_dir','retfit',[],...
  'GLM_ROI_retfit_flag',false,[false true],...
  'GLM_ROI_retfit_stem','retfit',[],...
  'GLM_ROI_retfit_thresh',0,[],...
  'GLM_ROI_norm_flag',true,[false true],...
  'GLM_ROI_iresp_flag',0,[0 1 2],...
  'GLM_ROI_iresp_t0',1,[1,Inf],...
  'GLM_ROI_iresp_t1',Inf,[1,Inf],...
  'GLM_ROI_iresp_baseline_flag',false,[false true],...
  'GLM_ROI_TR',1,[0.1,100],...
  'GLM_ROI_weights_flag',true,[false true],...
  'GLM_ROI_weights_dir',[],[],...
  'GLM_ROI_combine_ROIs_flag',false,[false true],...
  'GLM_ROI_names',{'v1','v2','v3'},[],...
  'GLM_ROI_cond_info',[],[],...
  'GLM_ROI_fname_conds',[],[],...
  'GLM_ROI_r_vec',7,[],...
  'GLM_ROI_th_vec',[45,135,225,315],[],...
  'GLM_ROI_ecc_width',10,[0,100],...
  'GLM_ROI_theta_width',90,[0,360],...
  'GLM_ROI_r_max',12.5,[0,Inf],...
  'GLM_ROI_rf_sizes',[0.66,1.03,1.88],[],...
  'GLM_ROI_rf_slopes',[0.06,0.10,0.15],[0,10],...
  'GLM_ROI_surround_flag',false,[false true],...
  'GLM_ROI_surround_rf_fact',1.5,[],...
  'GLM_ROI_surround_amp_fact',0.6,[],...
...
  'snums_valid',[],[],...
  'infix','corr_resBOLD',[],...
  'mc_flag',true,[false true],...
  'mc_inter_flag',true,[false true],...
  'regFS_flag',false,[false true],...
  'fstats_type',2,[0:2],...
  'cxfstatsflag',true,[false true],...
  'paint_flag',1,[0,2],...
  'projdist',1,[0,10],...
  'projfrac',0.5,[0,2],...
  'projfrac_flag',false,[false true],...
  'projdist_avg',[],[],...
  'projfrac_avg',[],[],...
  'resamp_flag',0,[0,1],...
  'force_repaint_flag',false,[false true],...
  'surfname','white',[],...
  'mask_midbrain_flag',false,[false true],...
  'sphere_flag',false,[false true],...
...
  'required_containers',{'proc_bold','fsurf'},[],...
  'modality','MRI',[],...
  'QC_raw',true,[false true],... % only applies if manual raw QC exists
  'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists
  'QC_recon',true,[false true],...
};
parms = mmil_args2parms(varargin,parms_filter);

% excl_tags are fields that should not be passed to MMIL_Analyze_BOLD_Exam
excl_tags = {'StudyInfo' 'RootDirs' 'batchname' 'qcflag'...
  'required_containers' 'modality' 'QC_raw' 'QC_BOLD' 'QC_recon'...
  'GLM_ROI_VisitID'};
tags = setdiff(fieldnames(parms),excl_tags);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_ProjID');
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

if strcmp(parms.infix,'none'), parms.infix = []; end;

% split GLM_fnames
StudyInfo = BOLD_MMIL_Check_StudyInfo(StudyInfo);

% cond_info file could be in ProjInfo or ContainerPath
if ~isempty(parms.GLM_ROI_fname_conds)
  if mmil_isrelative(parms.GLM_ROI_fname_conds)
    parms.GLM_ROI_fname_conds = [RootDirs.home '/ProjInfo/' ProjID '/'...
      parms.GLM_ROI_fname_conds];
  end;
  if ~exist(parms.GLM_ROI_fname_conds,'file')
    fprintf('%s: WARNING: missing cond_info file %s',...
      mfilename,parms.GLM_ROI_fname_conds);
  end;
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
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;
for i=1:length(StudyInfo)
  tmpInfo = StudyInfo(i);
  VisitID = tmpInfo.VisitID;
  GLM_ROI_VisitID = tmpInfo.GLM_ROI_VisitID;
  ContainerPath = fullfile(RootDirs.proc_bold,StudyInfo(i).proc_bold);
  FSContainerPath = fullfile(RootDirs.fsurf,StudyInfo(i).fsurf);
  tmp_args = {ContainerPath,FSContainerPath};
  tmp_parms = parms;

  for t=1:length(tags)
    tmp = mmil_getfield(tmpInfo,tags{t},[]);
    if ~isempty(tmp), tmp_parms = setfield(tmp_parms,tags{t},tmp); end;
  end;
  
  if isempty(tmp_parms.FP_snums), tmp_parms.FP_flag = 0; end;
  if isempty(tmp_parms.FE_snums), tmp_parms.FE_flag = 0; end;
  if isempty(tmp_parms.GLM_snums) |...
     ~isfield(tmp_parms,'GLM_stim_fnames') |...
      isempty(tmp_parms.GLM_stim_fnames)
    tmp_parms.GLM_flag = 0;
    tmp_parms.GLM_ROI_flag = 0;
  end;

  if isfield(tmpInfo,'BOLDScanNums'),
    tmp_parms.snums_valid = tmpInfo.BOLDScanNums;
  end;

  % create jobs only for subjects with FP, FE, or GLM snums
  if ~tmp_parms.FP_flag && ~tmp_parms.FE_flag &&...
     ~tmp_parms.GLM_flag && ~tmp_parms.GLM_ROI_flag
    continue;
  end;
  
  % set full path for directory of ROIs to use for GLM ROI analysis
  if tmp_parms.GLM_ROI_flag && ~isempty(tmp_parms.GLM_snums)
    if parms.GLM_ROI_retfit_flag
      tmp_parms.GLM_ROI_dir = [];
      if ~isempty(GLM_ROI_VisitID)
        tmp_path = MMIL_Get_Container(RootDirs,GLM_ROI_VisitID,'proc_bold');
      else
        tmp_path = ContainerPath;
      end;
      if isempty(tmp_path) || ~exist(tmp_path,'dir')
        fprintf('%s: WARNING: Container for GLM_ROI_retfit_dir %s not found\n',...
          mfilename,tmp_path);
        tmp_parms.GLM_ROI_retfit_flag = 0;
      else
        tmp_parms.GLM_ROI_retfit_dir = ...
          [tmp_path '/' tmp_parms.GLM_ROI_retfit_dir];
        if ~exist(tmp_parms.GLM_ROI_retfit_dir,'file')
          fprintf('%s: WARNING: GLM_ROI_retfit_dir %s not found\n',...
            mfilename,tmp_parms.GLM_ROI_retfit_dir);
          tmp_parms.GLM_ROI_retfit_flag = 0;
        end;
      end;
      tmp_val = mmil_getfield(StudyInfo(i),'RF_r_max',[]);
      if ~isempty(tmp_val)
        tmp_parms.GLM_ROI_r_max = tmp_val;
      end;
    else
      if ~isempty(GLM_ROI_VisitID)
        tmp_path = MMIL_Get_Container(RootDirs,GLM_ROI_VisitID);
      else    
        tmp_path = FSContainerPath;
      end;
      if isempty(tmp_path) || ~exist(tmp_path,'dir')
        fprintf('%s: WARNING: Container for GLM_ROI_dir %s not found\n',...
          mfilename,tmp_path);
        tmp_parms.GLM_ROI_dir = [];
      else
        tmp_parms.GLM_ROI_dir = [tmp_path '/' tmp_parms.GLM_ROI_dir];
        if ~exist(tmp_parms.GLM_ROI_dir,'file')
          fprintf('%s: WARNING: GLM_ROI_dir %s not found\n',...
            mfilename,tmp_parms.GLM_ROI_dir);
          tmp_parms.GLM_ROI_dir = [];
        end;
      end;
    end;
  else
    tmp_parms.GLM_ROI_dir = [];
    tmp_parms.GLM_ROI_retfit_flag = 0;
  end;    
  
  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem);
  jobfname = [batchdir '/' jobID '.m'];
  mmil_write_script(jobfname,'MMIL_Analyze_BOLD_Exam',tmp_args,tags,tmp_parms);

  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
  j=j+1;
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);

