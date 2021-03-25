function [fnames_mgz,fname_csv] = MMIL_Concat_Long_SurfStats(ProjID,varargin)
%function [fnames_mgz,fname_csv] = MMIL_Concat_Long_SurfStats(ProjID,[options])
%
% Usage:
%  [fname_mgz,fname_csv] = MMIL_Concat_Long_SurfStats(ProjID,'key1', value1,...);
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
%  'QC_raw' : [0|1] use good raw QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_RawQC.mat
%     {default = 1}
%  'QC_recon' : [0|1] use recon QC if it exists for this project
%     /home/{user}/ProjInfo/{ProjID}/{ProjID}_FSReconQC.csv
%     {default = 1}
%  'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
%     {default = 1}
%
% Optional Parameters that specify which analysis results to compile:
%  'datatype': type of surface stats data
%    allowed data types: {'dv'}
%     'dv':      longitudinal volume change
%    {default = 'dv'}
%  'analysis_outdir': where MMIL_Analyze_Longitudinal_Exam puts output
%     files relative to Longitudinal Container directory
%    {default = 'analysis'}
%  'visitA': for dv, baseline visit number
%     {default = 1}
%  'visitB': for dv, followup visit number
%     {default = 2}
%  'projdist_list': vector of projection distances for T1surf and T1cont
%    {default = [-0.2,0.2]}
%  'smoothing': surface smoothing steps (on sphere)
%     slope of FWHM vs. sqrt(N) is ~1.13 for fsaverage (v3)
%     (FWHM = full-width-half-max smoothing kernel
%         N = number of smoothing steps)
%      with [176,705,2819], approx FWHM (mm) = 15,30,60
%      with [100,300,3000], approx FWHM (mm) = 11.3, 19.6, 61.9
%    {default = 176}
%  'mask_midbrain_flag', [0|1] whether mid brain and other
%     cortical regions marked "unknown" were masked out
%    {default = 0}
%
% Optional Parameters:
%   'outdir': output directory
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'MRI_SurfStats'}
%   'outstem': output file stem
%     relative to outdir unless full path given
%     {default = 'Long'}
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
% Created:  10/01/14 by Don Hagler
% Last Mod: 01/22/15 by Don Hagler
%

%% todo: loop over combinations of VisitA and VisitB
%     may be a vector of visit numbers used as baseline
%     may be a vector of visit numbers used as followup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
fnames_mgz = []; fname_csv = [];

[parms,fnames_mgz,fname_csv] = check_input(ProjID,varargin);
if parms.nsubs==0, return; end;

fnames_mgz = concat_files(parms);

fname_csv = write_csv(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,fnames_mgz,fname_csv] = check_input(ProjID,options)
  fnames_mgz = []; fname_csv = [];
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'QC_raw',true,[false true],...
    'QC_recon',true,[false true],...
    'qcflag',true,[false true],...
... % specify which results to compile
    'datatype','thick',{'thick','area','cortvol','T1surf','T1cont','dv'},...
    'analysis_outdir','analysis',[],...
    'long_analysis_outdir','analysis',[],...
    'projdist_list',[-0.2,0.2],[-5,5],...
    'smoothing',176,[0,Inf],...
    'mask_midbrain_flag',false,[false true],...
    'visitA',1,[],...
    'visitB',2,[],...
... % other
    'outdir','MRI_SurfStats',[],...
    'outstem','MRI',[],...
    'meanflag',false,[false true],...
    'options',[],[],...
    'verbose',1,[0:2],...
    'forceflag',false,[false true],...
... % undocumented:
    'hemilist',{'lh','rh'},{'lh' 'rh'},...
    'modality','MRI',[],...
    'intype','mgz',{'mgh','mgz'},...
    'outtype','mgz',{'mgh','mgz'},...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
...
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'concat_tags',{'meanflag','options','forceflag'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  switch parms.datatype
    case 'dv'
      parms.required_containers = {'long'};
    otherwise
      parms.required_containers = {'fsurf'};
  end;
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(parms.StudyInfo), error('empty StudyInfo'); end;
  parms.nsubs = length(parms.StudyInfo);
  if parms.nsubs==0
    fprintf('%s: ERROR: no valid subjects\n',mfilename);
    return;
  end;

  % select visits
  VisitNumbers = [parms.StudyInfo.VisitNumber];
  ind = find(VisitNumbers==parms.visitB);
  parms.StudyInfo = parms.StudyInfo(ind);
  parms.nsubs = length(parms.StudyInfo);
  if parms.nsubs==0
    fprintf('%s: ERROR: no subjects with visitB = %d\n',...
      mfilename,parms.visitB);
    return;
  end;

  parms.nhemi = length(parms.hemilist);
  parms.nsmooth = length(parms.smoothing);

  % check output files
  [runflag,fnames_mgz,fname_csv] = check_output_files(parms);
  if ~runflag
    parms.nsubs = 0;
    return;
  end;
  
  % check data files, exclude subjects with missing data
  parms = check_data_files(parms);
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

  if strcmp(parms.datatype,'T1cont')
    parms.projdist_list = parms.projdist_list(parms.projdist_list>0);
    if length(parms.projdist_list) == 0
      error('datatype T1cont requires at least one positive projdist value');
    end;
  end;

  if ~strcmp(parms.outstem,parms.datatype)
    parms.outstem = [parms.outstem '_' parms.datatype];
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

function [runflag,fnames_mgz,fname_csv] = check_output_files(parms)
  fnames_mgz = []; fname_csv = [];
  runflag = 0;
  fname_csv = set_fname_csv(parms);
  if ~exist(fname_csv,'file'), runflag = 1; return; end;
  fnames_mgz = cell(parms.nhemi,parms.nsmooth);
  for h=1:parms.nhemi
    for k=1:parms.nsmooth
      fname_out = set_fname_out(parms,h,k);
      if ~exist(fname_out,'file'), runflag = 1; return; end;
      fnames_mgz{h,k} = fname_out;
    end;
  end;
  fnames_mgz = squeeze(fnames_mgz);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnames_mgz = concat_files(parms)
  fnames_mgz = cell(parms.nhemi,parms.nsmooth);
  for h=1:parms.nhemi
    for k=1:parms.nsmooth
      fname_out = set_fname_out(parms,h,k);
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: concatenating data files into %s...\n',...
            mfilename,fname_out);
        end;
        fnamelist = cell(parms.nsubs,1);
        fname = set_datafile(parms,h,k);
        for s=1:parms.nsubs
          fnamelist{s} = sprintf('%s/%s',...
            set_datadir(parms,s),fname);
        end;
        args = mmil_parms2args(parms,parms.concat_tags);
        mmil_concat(fnamelist,fname_out,args{:});
      end;
      fnames_mgz{h,k} = fname_out;
    end;
  end;
  fnames_mgz = squeeze(fnames_mgz);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_csv = write_csv(parms)
  fname_csv = [];
  fname_csv = set_fname_csv(parms);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_data_files(parms)
  if parms.verbose==2
    fprintf('%s: checking data files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    for h=1:parms.nhemi
      for k=1:parms.nsmooth
        fname = sprintf('%s/%s',...
          set_datadir(parms,s),set_datafile(parms,h,k));
        if ~exist(fname,'file')
          if parms.verbose
            fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
          end;
          %% todo: only require at least one file for this setting?
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function datadir = set_datadir(parms,s)
  datadir = [];
  switch parms.datatype
    case 'dv'
      datadir = sprintf('%s/%s/%s/visit_%d_VS_%d',...
        parms.RootDirs.long,parms.StudyInfo(s).long,...
        parms.long_analysis_outdir,...
        parms.visitA,parms.visitB);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_datafile(parms,h,k)
  fname = [];
  switch parms.datatype
    case 'dv'      % longitudinal volume change
      fname = 'dv_Fine_paint';
      if parms.mask_midbrain_flag, fname = [fname '-mbmask']; end;
      fname = [fname '-sphere'];
  end;
  if parms.smoothing(k)~=0
    fname = sprintf('%s-sm%d',fname,parms.smoothing(k));
  end;
  fname = sprintf('%s-%s.%s',fname,parms.hemilist{h},parms.intype);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_fname_out(parms,h,k)
  outstem = parms.outstem;
  if parms.mask_midbrain_flag
    outstem = [outstem '_mbmask'];
  end;
  outstem = sprintf('%s_sm%d',outstem,parms.smoothing(k));
  fname_out = [outstem '-' parms.hemilist{h} '.' parms.outtype];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_csv = set_fname_csv(parms)
  fname_csv = sprintf('%s_info.csv',parms.outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


