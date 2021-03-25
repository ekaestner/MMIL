function [fnames_mgz,fnames_csv,fnames_log] = ABCD_Concat_taskBOLD_SurfStats(ProjID,varargin)
%function [fnames_mgz,fnames_csv,fnames_log] = ABCD_Concat_taskBOLD_SurfStats(ProjID,[options])
%
% Usage:
%  [fnames_mgz,fnames_csv,fnames_log] = ABCD_Concat_taskBOLD_SurfStats(ProjID,'key1', value1,...)
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv')
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%      proc_bold, fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
%     {default = 1}
%  'max_motion': exclude subjects with mean relative motion greater than this
%     set to Inf to include all subjects
%     {default = 0.5}
%
% Optional Parameters that specify which analysis results to compile:
%   'tasknames': cell array of task names
%     {default = {'MID','SST','nBack'}}
%   'analysis_outfix': string attached to analysis output subdirectories
%     {default = 'analysis'}
%   'instem': taskBOLD analysis file stem
%     {default = 'taskBOLD'}
%   'infix': string inside BOLD file names (e.g. 'corr_resBOLD')
%     {default = 'corr_resBOLD'}
%   'ico_order': icosahedral order (0-7)
%     {default = 7}
%   'ico_presmooth': number of smoothing steps before resampling to ico
%     NOTE: FWHM ~ 1.25 * sqrt(N)
%     {default = 0}
%   'ico_smooth': number of smoothing steps after resampling to ico
%     NOTE: FWHM ~ 1.25 * sqrt(N)
%     {default = 64}
%
% Optional Parameters:
%   'outdir': output directory
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'taskBOLD_SurfStats'}
%   'outstem': output file stem
%     relative to outdir unless full path given
%     {default = 'taskBOLD'}
%   'meanflag': [0|1] calculate mean instead of concatenating
%     {default = 0}
%   'options': option string to use any of the mri_concat command line options
%     {default = []}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fnames_mgz: cell array with element for each task
%     containing cell arrays with size = [nhemi,nconds]
%   fnames_csv: cell array of csv file for each task
%   fnames_log: cell array with element for each task
%     containing cell arrays with size = [nhemi,nconds]
%     that contain names of files with lists of input file names
%
% Created:  04/07/17 by Don Hagler
% Prev Mod: 09/08/17 by Don Hagler
% Last Mod: 10/23/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
fnames_mgz = []; fnames_csv = []; fnames_log = [];

parms = check_input(ProjID,varargin);
if parms.nsubs==0, return; end;

% concatenate files separately for each task
for t=1:parms.ntasks
  taskname = parms.tasknames{t};
  % set condition names for this task
  parms = set_condnames(parms,t);
  % concatenate files across subjects  
  [fnames_mgz_tmp,fnames_log_tmp] = concat_files(parms,t);
  fnames_mgz{t} = fnames_mgz_tmp;
  fnames_log{t} = fnames_log_tmp;
  fnames_csv{t} = write_csv(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'qcflag',true,[false true],...
    'max_motion',0.5,[],...
... % specify which results to compile
    'tasknames',{'MID','SST','nBack'},{'MID','SST','nBack'},...
    'analysis_outfix','analysis',[],...
    'indir','taskBOLD_analysis',[],...
    'instem','taskBOLD',[],...
    'infix','corr_resBOLD',[],...
    'ico_order',7,[0 7],...
    'ico_presmooth',0,[0,1000],...
    'ico_smooth',64,[0,1000],...
...
    'outdir','taskBOLD_SurfStats',[],...
    'outstem','taskBOLD',[],...
    'meanflag',false,[false true],...
    'options',[],[],...
    'verbose',1,[0:2],...
    'forceflag',false,[false true],...
...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'intype','mgz',{'mgh','mgz'},...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
...
    'eprime_tags',{'rootdir','indir','fname_info','tasknames',...
                   'outdir','outstem','forceflag'},[],...
    'projinfo_tags',{'VisitIDs','SubjIDs','StudyInfo','RootDirs','ignore_VisitInfo_flag',...
                 'user','numvec_tags'},[],...
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'concat_tags',{'meanflag','options','forceflag'},[],...
... % hidden
    'required_containers',{'proc_bold','fsurf'},[],...
    'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists %% NOTE: this is not yet implemented
    'QC_recon',true,[false true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~iscell(parms.tasknames)
    parms.tasknames = {parms.tasknames};
  end;
  parms.ntasks = length(parms.tasknames);

  % load project and study info
  args = mmil_parms2args(parms,parms.projinfo_tags);
  [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Quick_Check_ProjID(ProjID,args{:});
  if isempty(parms.StudyInfo), error('empty StudyInfo'); end;
  parms.nsubs = length(parms.StudyInfo);
  parms.nhemi = length(parms.hemilist);

  % load eprime summary (used to get snums for each task)
  parms = check_eprime(parms);

  % exclude visits without motion files
  parms = check_motion_files(parms);
  if parms.nsubs==0
    fprintf('%s: ERROR: no subjects with valid motion files\n',mfilename);
    return;
  end;

  % exclude visits with missing files
  for t=1:parms.ntasks
    parms = set_condnames(parms,t);
    parms = check_cond_files(parms,t);
  end;
  if parms.nsubs==0
    fprintf('%s: ERROR: no subjects with valid files\n',mfilename);
    return;
  end;
  
  % get container-related info
  if parms.continfo_flag
    parms.StudyInfo = MMIL_Get_ContInfo(parms.StudyInfo,parms.RootDirs);
  end;

  % get subject-related info
  if parms.subjinfo_flag
    parms.StudyInfo = ...
      MMIL_Get_SubjInfo(parms.StudyInfo,parms.RootDirs,parms.ProjID);
  end;

  if mmil_isrelative(parms.outstem)
    if mmil_isrelative(parms.outdir)
      parms.outdir = sprintf('%s/MetaData/%s/%s',getenv('HOME'),ProjID,parms.outdir);
    end;
    parms.outstem = [parms.outdir '/' parms.outstem];
  else
    parms.outdir = fileparts(parms.outstem);
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_eprime(parms)
  % create spreadsheet with eprime file info  (or load existing)
  tp = [];
  tp.tasknames = parms.tasknames;
  tp.forceflag = 0;
  args = mmil_parms2args(tp,parms.eprime_tags);
  parms.fname_eprime = ABCD_Check_Eprime(parms.ProjID,args{:});
  parms.all_eprime_info = mmil_csv2struct(parms.fname_eprime);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnames_mgz,fnames_log] = concat_files(parms,t)
  fnames_mgz = []; fnames_log = [];
  fnames_mgz = cell(parms.nhemi,parms.nconds);
  fnames_log = cell(parms.nhemi,parms.nconds);
  for c=1:parms.nconds
    condname = parms.condnames{c};
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      fnamelist = cell(parms.nsubs,1);
      for s=1:parms.nsubs
        fnamelist{s} = set_cond_fname(parms,s,h,t,c);
      end;
      fnames_mgz{h,c} = sprintf('%s_%s-%s.mgz',...
        parms.outstem,condname,hemi);
      fnames_log{h,c} = sprintf('%s_%s-%s.log',...
        parms.outstem,condname,hemi);
      if ~exist(fnames_mgz{h,c},'file') ||...
          ~exist(fnames_log{h,c},'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: concatenating data files into %s...\n',...
            mfilename,fnames_mgz{h,c});
        end;
        args = mmil_parms2args(parms,parms.concat_tags);
        mmil_concat(fnamelist,fnames_mgz{h,c},args{:});
        write_log(fnames_log{h,c},fnamelist);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_log(fname_log,fnamelist)
  fid = fopen(fname_log,'wt');
  if fid==-1
    error('failed to open %s for writing',fname_log);
  end;
  fprintf(fid,'%s\n',fnamelist{:});
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_csv = write_csv(parms)
  fname_csv = [];
  fname_csv = sprintf('%s_info.csv',parms.outstem);
  if ~exist(fname_csv,'file') || parms.forceflag
    if parms.verbose==2
      fprintf('%s: writing csv file...\n',mfilename);
    end;
    % initialize data cell array
    data = cell(parms.nsubs,0);
    col_labels = cell(1,0);
    % add SubjIDs as first column    
    data = cat(2,data,{parms.StudyInfo.SubjID}');
    col_labels = cat(2,col_labels,'SubjID');
    % add VisitIDs as second column
    data = cat(2,data,{parms.StudyInfo.VisitID}');
    col_labels = cat(2,col_labels,'VisitID');
    % compile TR, nreps, numTRs
    %% NOTE: TR, nreps, and numTRs have size = [nsubs,ntasks]
    %% NOTE: nreps and numTRs are redundant except for APE sequence
    [TR,nreps,numTRs] = compile_scan_info(parms);
    for t=1:parms.ntasks
      tname = parms.tasknames{t};
      % add TR
      data = cat(2,data,num2cell(TR(:,t)));
      col_labels = cat(2,col_labels,['TR_' tname]);
      % add nreps
      data = cat(2,data,num2cell(nreps(:,t)));
      col_labels = cat(2,col_labels,['nreps_' tname]);
      % add numTRs
      data = cat(2,data,num2cell(numTRs(:,t)));
      col_labels = cat(2,col_labels,['numTRs_' tname]);
      % add motion
      motion = compile_motion(parms,t);
      data = cat(2,data,num2cell(motion));
      col_labels = cat(2,col_labels,['motion_' tname]);
    end;
    % add info from StudyInfo
    for t=1:length(parms.info_tags)
      tag = parms.info_tags{t};
      if isfield(parms.StudyInfo,tag)
        info = {parms.StudyInfo.(tag)}';
        data = cat(2,data,info);
        col_labels = cat(2,col_labels,tag);
      end;
    end;
    % add column labels
    data = cat(1,col_labels,data);
    % write file
    mmil_write_csv(fname_csv,data);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_cond_files(parms,t)
  if parms.verbose==2
    fprintf('%s: checking cond files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    for h=1:parms.nhemi
      for c=1:parms.nconds
        fname = set_cond_fname(parms,s,h,t,c);
        if isempty(fname)
          if parms.verbose
            fprintf('%s: WARNING: no cond file found\n',mfilename);
          end;
          valid_flags(s) = 0;
          break;
        elseif ~exist(fname,'file')
          if parms.verbose
            fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
          end;
          valid_flags(s) = 0;
          break;
        end;
      end;
    end;
  end;
  ind_valid = find(valid_flags);
  if length(ind_valid)<parms.nsubs
    parms.StudyInfo = parms.StudyInfo(ind_valid);
    parms.nsubs = length(parms.StudyInfo);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_motion_files(parms)
  if parms.verbose==2
    fprintf('%s: checking motion files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,parms.ntasks);
  for s=1:parms.nsubs
    VisitID = parms.StudyInfo(s).VisitID;
    for t=1:parms.ntasks
      taskname = parms.tasknames{t};
      fname = set_motion_fname(parms,s,t);
      if isempty(fname)
        if parms.verbose
          fprintf('%s: WARNING: no motion file found for %s %s\n',...
            mfilename,VisitID,taskname);
        end;
        valid_flags(s,t) = 0;
      elseif ~exist(fname,'file')
        if parms.verbose
          fprintf('%s: WARNING: file %s not found for %s %s\n',...
            mfilename,fname,VisitID,taskname);
        end;
%        keyboard
        valid_flags(s,t) = 0;
      elseif parms.max_motion<Inf
        motion = load_subj_motion(parms,s,t);
        if motion > parms.max_motion
          fprintf('%s: WARNING: excluding visit %s %s with excessive motion...\n',...
            mfilename,VisitID,taskname);
          valid_flags(s,t) = 0;
        end;    
      end;
    end;
  end;
  ind_valid = find(all(valid_flags,2));
  if length(ind_valid)<parms.nsubs
    parms.StudyInfo = parms.StudyInfo(ind_valid);
    parms.nsubs = length(parms.StudyInfo);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_condnames(parms,t)
  switch upper(parms.tasknames{t})
    case 'MID'
      [stim_labels,conds_contrast] = abcd_set_contrasts_mid;
    case 'SST'
      [stim_labels,conds_contrast] = abcd_set_contrasts_sst;
    case 'NBACK'
      [stim_labels,conds_contrast] = abcd_set_contrasts_nback;
    otherwise
      error('prep_stims not implemented for %s',parms.tasknames{t});
  end;
  parms.condnames = {conds_contrast.name};
  parms.nconds = length(parms.condnames);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function snums = set_snums(parms,s,t)
  snums = [];
  VisitID = parms.StudyInfo(s).VisitID;
  % find containers for this VisitID
  ContainerPath = MMIL_Get_Container(parms.RootDirs,VisitID,'proc_bold');
%  RawContainerPath = MMIL_Get_Container(parms.RootDirs,VisitID,'raw');
  FSContainerPath = MMIL_Get_Container(parms.RootDirs,VisitID,'fsurf');
  if isempty(ContainerPath)
    fprintf('%s: WARNING: missing proc_bold container for %s\n',...
      mfilename,VisitID);
    return;
  end;
%  if isempty(RawContainerPath)
%    fprintf('%s: WARNING: missing raw container for %s\n',...
%      mfilename,VisitID);
%    return;
%  end;
  if isempty(FSContainerPath)
    fprintf('%s: WARNING: missing fsurf container for %s\n',...
      mfilename,VisitID);
    return;
  end;

  % reduce eprime_info to this task
  parms.taskname = parms.tasknames{t};
  all_tasknames = {parms.all_eprime_info.taskname};
  istask = strcmp(parms.taskname,all_tasknames);
  eprime_info = parms.all_eprime_info(istask);

  % find eprime entries for this VisitID
  eprime_VisitIDs = {eprime_info.VisitID};
  ind_eprime = find(strcmp(VisitID,eprime_VisitIDs));
  if isempty(ind_eprime)
    if parms.verbose
      fprintf('%s: WARNING: no series for VisitID %s... skipping\n',...
        mfilename,VisitID);
    end;
    return;
  end;
  eprime_SeUIDs = {eprime_info(ind_eprime).SeriesInstanceUID};
  eprime_files = {eprime_info(ind_eprime).eprime_file_name};
  if any(cellfun(@isempty,eprime_files))
    if parms.verbose
      fprintf('%s: WARNING: no eprime files for VisitID %s... skipping\n',...
        mfilename,VisitID);
    end;
    return;
  end;  
  fname_eprime = unique(eprime_files);
  if length(fname_eprime)==1
    fname_eprime = fname_eprime{1};
  end;

  % load ContainerInfo from raw and proc containers
%  [RawContainerInfo,errcode] = MMIL_Load_ContainerInfo(RawContainerPath);
%  if errcode, return; end;
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
  if errcode, return; end;

  % match scans based on SeUID
%  SeUIDs = {RawContainerInfo.SeriesInfo.SeriesInstanceUID};
  SeUIDs = {ContainerInfo.SeriesInfo.SeriesInstanceUID};
  [tmp,ind_eprime,ind_series] = intersect(eprime_SeUIDs,SeUIDs);
  if isempty(ind_series)
    if parms.verbose
      fprintf('%s: WARNING: mismatch in SeriesInstanceUIDs for VisitID %s... skipping\n',...
        mfilename,VisitID);
    end;
%    keyboard
    return;
  end;
  
  % find corresponding BOLD snums
  ind_BOLD = [ContainerInfo.ScanInfo.BOLD.SeriesIndex];
  [tmp,ind_s,ind_B] = intersect(ind_series,ind_BOLD);
  snums = ind_B;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [instem,indir] = set_input(parms,s,t)
  instem = []; indir = [];
  VisitID = parms.StudyInfo(s).VisitID;
  ContainerPath = MMIL_Get_Container(parms.RootDirs,VisitID,'proc_bold');
  tp = [];
  tp.snums = set_snums(parms,s,t);
  if isempty(tp.snums), return; end;
  tp.dirstem = [parms.instem '_' parms.tasknames{t}];
  tp.infix = parms.infix;
  tp.sphere_flag = 1;
  tp.smoothsteps = parms.ico_presmooth;
  tp.sphsmoothsteps = parms.ico_smooth;
  tp.analysis_outfix = parms.analysis_outfix;
  args = mmil_parms2args(tp);
  [instem,errcode] = BOLD_MMIL_Set_GLM_Stem(ContainerPath,args{:});
  indir = fileparts(instem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_cond_fname(parms,s,h,t,c)
  fname = [];
  instem = set_input(parms,s,t);
  [indir,instem] = fileparts(instem);
  if isempty(instem), return; end;
  contname = parms.condnames{c};
  fname = sprintf('%s/contrasts/%s_%s-%s.%s',...
    indir,instem,contname,parms.hemilist{h},parms.intype);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_info_fname(parms,s,t)
  fname = [];
  instem = set_input(parms,s,t);
  if isempty(instem), return; end;
  instem = regexprep(instem,'_3dDeconv.+','_info');
  fname = sprintf('%s.mat',instem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_motion_fname(parms,s,t)
  fname = [];
  instem = set_input(parms,s,t);
  if isempty(instem), return; end;
  instem = regexprep(instem,'_3dDeconv.+','_motion');
  fname = sprintf('%s.mat',instem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TR,nreps,numTRs] = load_subj_scan_info(parms,s,t)
  TR = nan; nreps = nan; numTRs = nan;
  fname = set_info_fname(parms,s,t);
  if isempty(fname), return; end;  
  info = load(fname);
  TR = mean(info.scan_info.TRs);
  nreps = sum(info.scan_info.nreps);
  numTRs = sum(info.scan_info.numTRs);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mean_motion = load_subj_motion(parms,s,t)
  mean_motion = nan;
  fname = set_motion_fname(parms,s,t);
  if isempty(fname), return; end;
  load(fname);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TR,nreps,numTRs] = compile_scan_info(parms)
  TR = nan(parms.nsubs,parms.ntasks);
  nreps = nan(parms.nsubs,parms.ntasks);
  numTRs = nan(parms.nsubs,parms.ntasks);
  for s=1:parms.nsubs
    for t=1:parms.ntasks
      [TR(s,t),nreps(s,t),numTRs(s,t)] = load_subj_scan_info(parms,s,t);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion = compile_motion(parms,t)
  motion = nan(parms.nsubs,1);
  for s=1:parms.nsubs
    motion(s) = load_subj_motion(parms,s,t);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

