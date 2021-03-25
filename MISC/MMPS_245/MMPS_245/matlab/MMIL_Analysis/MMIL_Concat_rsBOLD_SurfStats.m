function [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_rsBOLD_SurfStats(ProjID,varargin)
%function [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_rsBOLD_SurfStats(ProjID,[options])
%
% Usage:
%  [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_rsBOLD_SurfStats(ProjID,'key1', value1,...);
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
%   'indir': rsBOLD analysis subdirectory inside each BOLDPROC dir
%     {default = 'rsBOLD_analysis'}
%   'instem': rsBOLD analysis file stem
%     {default = 'rsBOLD'}
%   'roi_infix': additional string added to ROI output files
%     {default = []}
%   'corr_infix': string attached to rsBOLD correlation output files
%     {default = []}
%   'ico_order': icosahedral order (0-7)
%     {default = 7}
%   'ico_presmooth': number of smoothing steps applied on native surface
%     before resampling to ico
%     NOTE: FWHM ~ 1.25 * sqrt(N)
%     {default = 64}
%   'seed_roinames': name of seed ROI used to calculate correlation map
%     if empty, will load variance map instead of correlation map
%     may be cell array for multiple seeds
%     rsBOLD analysis must have been run with corr_roi_ico_maps_flag=1
%     {default = []}
%
% Optional Parameters:
%   'outdir': output directory
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'rsBOLD_SurfStats'}
%   'outstem': output file stem
%     relative to outdir unless full path given
%     {default = 'rsBOLD'}
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
%   fnames_mgz: cell array of output mgz file names
%     that contain concatenated surface stats
%     size = [nhemi,nseeds] or [nhemi,1] if seed_roinames is empty
%   fname_csv: output csv file name containing relevant info for each visit
%   fnames_log: cell array of output log file names
%     that contain lists of input file names
%     size = [nhemi,nseeds] or [nhemi,1] if seed_roinames is empty
%
% Created:  10/29/12 by Don Hagler
% Last Mod: 01/10/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
fnames_mgz = []; fname_csv = []; fnames_log = [];

parms = check_input(ProjID,varargin);
if parms.nsubs==0, return; end;

[fnames_mgz,fnames_log] = concat_files(parms);

fname_csv = write_csv(parms);

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
    'indir','rsBOLD_analysis',[],...
    'instem','rsBOLD',[],...
    'roi_infix',[],[],...
    'corr_infix',[],[],...
    'ico_order',7,[0 7],...
    'ico_presmooth',64,[0,1000],...
    'seed_roinames',[],[],...
...
    'outdir','rsBOLD_SurfStats',[],...
    'outstem','rsBOLD',[],...
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
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'concat_tags',{'meanflag','options','forceflag'},[],...
... % hidden
    'required_containers',{'proc_bold','fsurf'},[],...
    'QC_raw',true,[false true],... % only applies if manual raw QC exists
    'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists
    'QC_recon',true,[false true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(parms.StudyInfo), error('empty StudyInfo'); end;
  parms.nsubs = length(parms.StudyInfo);

  parms.instem = [parms.indir '/' parms.instem];
  parms.nhemi = length(parms.hemilist);
  if ~isempty(parms.roi_infix)
    parms.roi_instem = [parms.instem '_' parms.roi_infix];
  else
    parms.roi_instem = parms.instem;
  end;
  if ~isempty(parms.corr_infix)
    parms.corr_instem = [parms.roi_instem '_' parms.corr_infix];
  else
    parms.corr_instem = parms.roi_instem;
  end;

  parms = check_motion_files(parms);
  if parms.nsubs==0
    fprintf('%s: ERROR: no valid subjects\n',mfilename);
    return;
  end;

  if isempty(parms.seed_roinames)
    parms.nseeds = 0;
    parms = check_var_files(parms);
  else
    if ~iscell(parms.seed_roinames)
      parms.seed_roinames = {parms.seed_roinames};
    end;
    parms.nseeds = length(parms.seed_roinames);
    parms = check_corr_files(parms);
  end;
  if parms.nsubs==0
    fprintf('%s: ERROR: no valid subjects\n',mfilename);
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
      parms.outdir = [getenv('HOME') '/MetaData/' ProjID '/' parms.outdir];
    end;
    parms.outstem = [parms.outdir '/' parms.outstem];
  else
    parms.outdir = fileparts(parms.outstem);
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnames_mgz,fnames_log] = concat_files(parms)
  fnames_mgz = []; fnames_log = [];
  if parms.nseeds==0
    fnames_mgz = cell(parms.nhemi,1);
    fnames_log = cell(parms.nhemi,1);
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      fnamelist = cell(parms.nsubs,1);
      for s=1:parms.nsubs
        fnamelist{s} = set_var_fname(parms,s,h);
      end;
      fnames_mgz{h} = sprintf('%s-%s.mgz',parms.outstem,hemi);
      fnames_log{h} = sprintf('%s-%s.log',parms.outstem,hemi);
      if ~exist(fnames_mgz{h},'file') ||...
         ~exist(fnames_log{h},'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: concatenating data files into %s...\n',...
            mfilename,fnames_mgz{h});
        end;
        args = mmil_parms2args(parms,parms.concat_tags);
        mmil_concat(fnamelist,fnames_mgz{h},args{:});
        write_log(fnames_log{h},fnamelist);
      end;
    end;
  else
    fnames_mgz = cell(parms.nhemi,parms.nseeds);
    fnames_log = cell(parms.nhemi,parms.nseeds);
    for r=1:parms.nseeds
      seed_roiname = parms.seed_roinames{r};
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        fnamelist = cell(parms.nsubs,1);
        for s=1:parms.nsubs
          fnamelist{s} = set_corr_fname(parms,s,h,r);
        end;
        fnames_mgz{h,r} = sprintf('%s_%s-%s.mgz',...
          parms.outstem,seed_roiname,hemi);
        fnames_log{h,r} = sprintf('%s_%s-%s.log',...
          parms.outstem,seed_roiname,hemi);
        if ~exist(fnames_mgz{h,r},'file') ||...
           ~exist(fnames_log{h,r},'file') || parms.forceflag
          if parms.verbose==2
            fprintf('%s: concatenating data files into %s...\n',...
              mfilename,fnames_mgz{h,r});
          end;
          args = mmil_parms2args(parms,parms.concat_tags);
          mmil_concat(fnamelist,fnames_mgz{h,r},args{:});
          write_log(fnames_log{h,r},fnamelist);
        end;
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
    [TR,nreps,numTRs] = compile_scan_info(parms);
    % add TR
    data = cat(2,data,num2cell(TR));
    col_labels = cat(2,col_labels,'TR');
    % add nreps
    data = cat(2,data,num2cell(nreps));
    col_labels = cat(2,col_labels,'nreps');
    % add numTRs
    data = cat(2,data,num2cell(numTRs));
    col_labels = cat(2,col_labels,'numTRs');
    % add motion
    motion = compile_motion(parms);
    data = cat(2,data,num2cell(motion));
    col_labels = cat(2,col_labels,'motion');
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

function parms = check_var_files(parms)
  if parms.verbose==2
    fprintf('%s: checking var files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    for h=1:parms.nhemi
      fname = set_var_fname(parms,s,h);
      if ~exist(fname,'file')
        if parms.verbose
          fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
        end;
        valid_flags(s) = 0;
        break;
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

function parms = check_corr_files(parms)
  if parms.verbose==2
    fprintf('%s: checking corr files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    for h=1:parms.nhemi
      for r=1:parms.nseeds
        fname = set_corr_fname(parms,s,h,r);
        if ~exist(fname,'file')
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

function parms = check_motion_files(parms,results)
  if parms.verbose==2
    fprintf('%s: checking motion files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    fname = set_motion_fname(parms,s);
    if ~exist(fname,'file')
      if parms.verbose
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
      end;
      valid_flags(s) = 0;
    elseif parms.max_motion<Inf
      motion = load_subj_motion(parms,s);
      if motion > parms.max_motion
        fprintf('%s: WARNING: excluding visit %s with excessive motion...\n',...
          mfilename,parms.StudyInfo(s).VisitID);
        valid_flags(s) = 0;
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

function fname = set_corr_fname(parms,s,h,r)
  fname = sprintf('%s/%s/%s_%s_ico%d_corr-%s.%s',...
    parms.RootDirs.proc_bold,parms.StudyInfo(s).proc_bold,...
    parms.corr_instem,parms.seed_roinames{r},parms.ico_order,...
    parms.hemilist{h},parms.intype);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_var_fname(parms,s,h)
  fname = sprintf('%s/%s/%s-sm%d-ico%d-var-%s.%s',...
    parms.RootDirs.proc_bold,parms.StudyInfo(s).proc_bold,...
    parms.instem,parms.ico_presmooth,parms.ico_order,...
    parms.hemilist{h},parms.intype);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_info_fname(parms,s)
  fname = [parms.RootDirs.proc_bold '/' parms.StudyInfo(s).proc_bold '/' ...
           parms.instem '_info.mat'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_motion_fname(parms,s)
  fname = [parms.RootDirs.proc_bold '/' parms.StudyInfo(s).proc_bold '/' ...
           parms.instem '_motion.mat'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TR,nreps,numTRs] = load_subj_scan_info(parms,s)
  numTRs = nan;
  fname = set_info_fname(parms,s);
  load(fname);
  TR = mean(info.TRs);
  nreps = sum(info.nreps);
  numTRs = sum(info.numTRs);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mean_motion = load_subj_motion(parms,s)
  mean_motion = nan;
  fname = set_motion_fname(parms,s);
  load(fname);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TR,nreps,numTRs] = compile_scan_info(parms)
  TR = nan(parms.nsubs,1);
  nreps = nan(parms.nsubs,1);
  numTRs = nan(parms.nsubs,1);
  for s=1:parms.nsubs
    [TR(s),nreps(s),numTRs(s)] = load_subj_scan_info(parms,s);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion = compile_motion(parms)
  motion = nan(parms.nsubs,1);
  for s=1:parms.nsubs
    motion(s) = load_subj_motion(parms,s);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

