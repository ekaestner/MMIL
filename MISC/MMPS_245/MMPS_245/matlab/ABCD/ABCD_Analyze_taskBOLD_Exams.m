function ABCD_Analyze_taskBOLD_Exams(ProjID,varargin)
%function ABCD_Analyze_taskBOLD_Exams(ProjID,[options])
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
%    SubjID, VisitID
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: proc_bold, fsurf (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied, MMIL_ProjInfo.csv is not required
%    {default = []}
%  'force_eprime_flag': [0|1] whether to create eprime info file with forceflag = 1
%    {default = 0}
%  'eprime_rootdir': root input directory containing ProjID dirs
%     {default = '/space/syn05/1/data/MMILDB'}
%  'eprime_indir' root input directory containing eprime files
%     if empty, will be in {rootdir}/{ProjID}/aux_incoming
%     {default = []}
%  'qcflag': [0|1] whether to exclude subjects with StudyInfo.QC=0
%    {default = 1}
%  'batchname': name of output batchdir
%    {default = 'MMIL_Analyze_taskBOLD_Exams'}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional Parameters applicable to data selection:
%  'tasknames': cell array of task names
%    {default = {'MID','SST','nBack'}}
%  'infix': input file will be sprintf('BOLD%d_%s.mgz',snum,infix)
%    to specify no infix (i.e. raw data), use 'none'
%    {default = 'corr_resBOLD'}
%  'mc_flag': [0|1] whether within-scan motion correction was done
%    Will use motion.1D files as nuisance regressors
%    {default = 1}
%  'mc_inter_flag': [0|1] whether between-scan motion correction was done
%    Allows for a single reference scan to be used for registration to FS
%    {default = 1}
%
% Optional Parameters specific to GLM analysis:
%  'analysis_prefix': string at beginning of analysis output subdirectories
%    {default = 'taskBOLD'}
%  'analysis_outfix': string attached to analysis output subdirectories
%    {default = 'analysis'}
%  'stim_times_flag': [0|1|2] use stimulus time files
%     0 = .1D, 1 = .txt and 2 = _block.txt 
%     {default = 1}
%  'stim_times_model': name of stimulus model
%    model names in 3dDeconvolve for -stim_times
%    allowed: 'SPMG','GAM','TENT', 'BLOCK'
%    {default = 'SPMG'}
%  'stim_times_nbasis': max number of basis functions for HRF
%    {default = 10}
%  'concat_flag': [0|1|2] analyze concatenated across scans
%    0: analyze each scan individually
%    1: analyze concatenated scans
%    2: analyze individually and concatenated
%    {default = 1}
%  'pthresh': probability threshold applied to GLM f-stats for each
%     condition
%    {default = 0.05}
%  'contrasts_flag': [0|1] calculate glt contrasts between each condition
%    {default = 0}
%  'iresp_flag': [0|1] output impulse response functions for each condition
%    {default = 0}
%
% Optional Parameters applicable to surface-based analysis:
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
%  'aseg_flag': [0|1] extract ROI averages for subcortical segmentation
%    {default = 0}
%  'aseg_erode_flag': [0|1] "erode" aseg ROIs by smoothing and thresholding
%    after resampling to BOLD resolution
%    {default = 0}
%  'aseg_erode_nvoxels': smoothing kernel sigma for aseg erosion (# of voxels)
%    {default = 1}
%  'fname_aseg': name of aseg file
%    if empty, will use aseg.mgz in subjdir/subj/mri
%    {default = []}
%  'aparc_flag': [0|1] extract ROI averages for cortical surface parcellation
%    {default = 0}
%  'fnames_aparc': cell array of annotation files
%    if empty, will use ?h.aparc.annot files in subjdir/subj/label
%     otherwise, must use ?h.<name>.annot naming convention
%    {default = []}
%  'fparc_flag': extract time series for fparc cortical surface ROIs
%    ignored if fnames_fparc is empty
%    {default = 1}
%  'fnames_fparc': cell array of annotation files in fsaverage space
%    will be resampled to individual subject space before use
%    {default = []}
%  'render_flag': [0|1] whether to generate images of surface maps
%     {default = 0}
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
% Created:  12/02/16 by Don Hagler
% Prev Mod: 10/25/17 by Dani Cornejo 
% Last Mod: 10/31/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchname','ABCD_Analyze_taskBOLD_Exams',[],...
  'fname_info',[],[],...
  'force_eprime_flag',false,[false true],...
  'eprime_rootdir','/space/syn05/1/data/MMILDB',[],...
  'eprime_indir',[],[],...
  'qcflag',true,[false true],...
  'forceflag',false,[false true],...
...
  'tasknames',{'MID','SST','nBack'},{'MID','SST','nBack'},...
  'infix','corr_resBOLD',[],...
  'mc_flag',true,[false true],...
  'mc_inter_flag',true,[false true],...
...
  'analysis_prefix','taskBOLD',[],...
  'analysis_outfix','analysis',[],...
  'stim_times_flag',1,[0:2],...
  'stim_times_model','SPMG',{'SPMG','TENT','GAM','BLOCK'},...
  'stim_times_nbasis',10,[1,100],...
  'concat_flag',1,[0:2],...
  'pthresh',0.05,[0,1],...
  'contrasts_flag',false,[false true],...
  'iresp_flag',false,[false true],...
...
  'regFS_flag',false,[false true],...
  'paint_flag',1,[0,2],...
  'aseg_flag',false,[false true],...
  'aseg_erode_flag',false,[false true],...
  'aseg_erode_nvoxels',1,[1:100],...
  'fname_aseg',[],[],...
  'aparc_flag',false,[false true],...
  'fnames_aparc',[],[],...
  'fparc_flag',true,[false true],...
  'fnames_fparc',[],[],...
  'render_flag',false,[false true],...
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
  'snums',[],[],...
  'snums_valid',[],[],...
  'skipTRs',0,[0,1000],...
  'minlag',0,[0,10],...
  'maxlag',14,[0,30],...
... % parameters for rendering output
  'curvflag',true,[false true],...
  'curvfact',0.2,[0,1],...
  'zoom',1,[0.01,100],...
  'tif_flag',true,[false true],...
  'tif_dpi',300,[10,1000],...
  'fmin',0.1,[0,100],...
  'fmid',1,[0.001,1000],...
  'fmax',2,[0.001,1000],...
  'view_surfname','inflated',{'inflated','white','pial'},...
  'viewlist',{'lat','med'},{'left','right','lat','med','pos',...
                            'ant','sup','ven','infr','infl'},...
...
  'verbose',2,[0:2],...
  'stim_fnames',[],[],...
  'smoothsteps',0,[0,1000],...
  'sphsmoothsteps',0,[0,1000],...
...
  'required_containers',{'proc_bold','fsurf'},[],...
  'modality','MRI',[],...
  'QC_raw',true,[false true],... % only applies if manual raw QC exists
  'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists
  'QC_recon',true,[false true],...
...
  'skipTRs_GE',5,[],...
  'skipTRs_Philips',8,[],...
  'skipTRs_Siemens',8,[],...
...
  'info_tags',{'VisitIDs','SubjIDs','StudyInfo','RootDirs','ignore_VisitInfo_flag',...
               'user','numvec_tags'},[],...
  'analy_tags',{'outstem','snums','analysis_outfix',...
                'stim_times_flag','stim_times_model','stim_times_nbasis',...
                'concat_flag',...
                'contrasts_flag','iresp_flag','skipTRs',...
                'minlag','maxlag','stim_fnames','smoothsteps',...
                'sphsmoothsteps','pthresh','snums_valid','infix','mc_flag','mc_inter_flag',...
                'regFS_flag','paint_flag',...
                'aseg_flag','aseg_erode_flag','aseg_erode_nvoxels',...
                'fname_aseg','aparc_flag','fnames_aparc',...
                'fparc_flag','fnames_fparc',...
                'render_flag',...
                'projdist','projfrac','projfrac_flag','projdist_avg',...
                'projfrac_avg','resamp_flag','force_repaint_flag','surfname',...
                'mask_midbrain_flag','sphere_flag','forceflag','fnamestem','out_ext',...
                'curvflag','curvfact','zoom','tif_flag','tif_dpi',...
                'fmin','fmid','fmax','view_surfname','viewlist'},[],...
};

parms = mmil_args2parms(varargin,parms_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = mmil_parms2args(parms,parms.info_tags);
[ProjInfo,StudyInfo,RootDirs] = MMIL_Quick_Check_ProjID(ProjID,args{:});

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

if ~iscell(parms.tasknames)
  parms.tasknames = {parms.tasknames};
end;
ntypes = length(parms.tasknames);

% create spreadsheet with names of eprime files
tparms = [];
%tparms.tasknames = parms.tasknames;
tparms.forceflag = parms.force_eprime_flag;
tparms.rootdir = parms.eprime_rootdir;
tparms.indir = parms.eprime_indir;
tparms.scanner_flag = 1;
args = mmil_parms2args(tparms);
fname_eprime_info = ABCD_Check_Eprime(ProjID,args{:});
all_eprime_info = mmil_csv2struct(fname_eprime_info);

all_VisitIDs = {StudyInfo.VisitID};

j = 1;
for k=1:ntypes
  taskname = parms.tasknames{k};
  parms.outstem = sprintf('%s_%s',parms.analysis_prefix,taskname);

  % reduce eprime_info to this task
  ind_task = find(strcmp({all_eprime_info.SeriesType},...
                  sprintf('fMRI_%s_task',taskname)));
  eprime_info = all_eprime_info(ind_task);

  % exclude entries with empty VisitID
  ind_valid = find(~cellfun(@isempty,{eprime_info.VisitID}));
  if length(ind_valid) < length(eprime_info)
    eprime_info = eprime_info(ind_valid);
  end;

  ind_folder = find([eprime_info.eprime_folder_found]);
  fprintf('%s: e-prime folder found for %d of %d series\n',...
    mfilename,length(ind_folder),length(eprime_info));

  ind_file = find([eprime_info.eprime_file_found]);
  fprintf('%s: e-prime file found for %d series\n',...
    mfilename,length(ind_file));
  
  ind_naming = find([eprime_info.eprime_file_naming]);
  fprintf('%s: correct naming of e-prime file for %d series\n',...
    mfilename,length(ind_naming));

  % reduce eprime_info to those with a file
  eprime_info = eprime_info(ind_file);
  eprime_VisitIDs = {eprime_info.VisitID};

  VisitIDs = intersect(eprime_VisitIDs,all_VisitIDs);
  nvisits = length(VisitIDs);
  fprintf('%s: attempting to create jobs for %d visits...\n',...
    mfilename,nvisits);
  
  for i=1:nvisits
    VisitID = VisitIDs{i};

    % find containers for this VisitID
    ContainerPath = MMIL_Get_Container(RootDirs,VisitID,'proc_bold');
%    RawContainerPath = MMIL_Get_Container(RootDirs,VisitID,'raw');
    FSContainerPath = MMIL_Get_Container(RootDirs,VisitID,'fsurf');

    if isempty(ContainerPath)
      fprintf('%s: WARNING: missing proc_bold container for %s\n',...
        mfilename,VisitID);
      continue;
    end;
%    if isempty(RawContainerPath)
%      fprintf('%s: WARNING: missing raw container for %s\n',...
%        mfilename,VisitID);
%      continue;
%    end;
    if isempty(FSContainerPath)
      fprintf('%s: WARNING: missing fsurf container for %s\n',...
        mfilename,VisitID);
      continue;
    end;

    % find eprime entries for this VisitID
    ind_eprime = find(strcmp(VisitID,eprime_VisitIDs));
    eprime_SeUIDs = {eprime_info(ind_eprime).SeriesInstanceUID};
    eprime_files = {eprime_info(ind_eprime).eprime_file_name};
    fname_eprime = unique(eprime_files);
    if length(fname_eprime)==1
      fname_eprime = fname_eprime{1};
    end;
    
    targs = {ContainerPath,FSContainerPath,fname_eprime};

    % load ContainerInfo from raw and proc containers
%    [RawContainerInfo,errcode] = MMIL_Load_ContainerInfo(RawContainerPath);
%    if errcode, continue; end;
    [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
    if errcode, continue; end;
    
    % set skipTRs according to Scanner Manufacturer
    switch ContainerInfo.Manufacturer
      case 'GE MEDICAL SYSTEMS'
        parms.skipTRs = parms.skipTRs_GE;
      case 'Philips Medical Systems'
        parms.skipTRs = parms.skipTRs_Philips;
      case 'SIEMENS'
        parms.skipTRs = parms.skipTRs_Siemens;
      otherwise
        fprintf('%s: WARNING: unsupported manufacturer "%s", setting skipTRs = 0\n',...
          mfilename,ContainerInfo.Manufacturer);
        parms.skipTRs = 0;
    end;
    
    % match scans based on SeUID
%    SeUIDs = {RawContainerInfo.SeriesInfo.SeriesInstanceUID};
    SeUIDs = {ContainerInfo.SeriesInfo.SeriesInstanceUID};
    [tmp,ind_eprime,ind_series] = intersect(eprime_SeUIDs,SeUIDs);
    if isempty(ind_series)
      fprintf('%s: WARNING: mismatch in SeriesInstanceUIDs for VisitID %s... skipping\n',...
        mfilename,VisitID);
      continue;
    end;
    
    % find corresponding BOLD snums
    ind_BOLD = [ContainerInfo.ScanInfo.BOLD.SeriesIndex];
    [tmp,ind_s,ind_B] = intersect(ind_series,ind_BOLD);
    parms.snums = ind_B;
    
    jstem = regexprep(VisitID,'\^','_');
    jstem = [jstem  '_' taskname];
    jstem = jstem(1:min(20,length(jstem)));
    jobID = sprintf('job_%03d_%s',j,jstem);
    jobfname = [batchdir '/' jobID '.m'];
    mmil_write_script(jobfname,'ABCD_Analyze_taskBOLD_Exam',...
      targs,parms.analy_tags,parms);

    fid = fopen(scriptlistfname,'a');
    fprintf(fid,'%s\n',jobID);
    fclose(fid);
    j=j+1;
  end;
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs3 %s\n',parms.batchname);

