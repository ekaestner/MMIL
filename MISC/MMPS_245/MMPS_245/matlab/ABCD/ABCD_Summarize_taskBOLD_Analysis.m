function ABCD_Summarize_taskBOLD_Analysis(ProjID,varargin)
%function ABCD_Summarize_taskBOLD_Analysis(ProjID,[options])
%
% Usage:
%  ABCD_Summarize_taskBOLD_Analysis(ProjID,'key1', value1,...)
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
%     {default = Inf}
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
%   'concat_flag': [0|1] summarize results from analysis of concatenated scans
%     0: summarize individual scan analysis
%     1: summarize concated scans analysis
%     {default = 0}
%   'max_nruns': maximum number of runs if concat_flag = 0
%     {default = 2}
%   'avg_flag': [0|1] average results across scans
%     ignored if concat_flag = 1 or if max_nruns < 2
%     {default = 1}
%   'aseg_flag': [0|1] compile aseg ROI results
%     {default = 1}
%   'aseg_erode_flag': [0|1] compile results from eroded aseg ROIs
%     {default = 1}
%   'aparc_flag': [0|1] compile aparc ROI results
%     {default = 1}
%   'fparc_flag': compile results for fparc cortical surface ROIs
%     {default = 1}
%   'fnames_fparc': cell array of annotation files in fsaverage space
%     used to determine fparc roinames and roicodes
%     {default = []}
%   'motion_flag': [0|1] compile head motion measures
%     {default = 1}
%
% Optional Parameters:
%   'outdir': output directory
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'ROI_Summaries'}
%   'outstem': output file stem
%     relative to outdir unless full path given
%     {default = 'taskBOLD'}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  08/29/17 by Don Hagler
% Prev Mod: 09/19/17 by Don Hagler
% Last Mod: 10/24/17 by Don Hagler
%

%% todo: compile dof values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(ProjID,varargin);
if parms.quitflag, return; end;

% set roi parameters
parms = set_roi_parms(parms);

% check files for each subject
parms = check_results(parms);
if parms.nsubs==0, return; end;

% create ROI summmary file
summarize_results(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'qcflag',true,[false true],...
    'max_motion',Inf,[],...
... % specify which results to compile
    'tasknames',{'MID','SST','nBack'},{'MID','SST','nBack'},...
    'analysis_outfix','analysis',[],...
    'indir','taskBOLD_analysis',[],...
    'instem','taskBOLD',[],...
    'infix','corr_resBOLD',[],...
    'concat_flag',0,[0:1],...
    'max_nruns',2,[1:10],...
    'avg_flag',true,[false true],...
    'aseg_flag',true,[false true],...
    'aseg_erode_flag',true,[false true],...
    'aparc_flag',true,[false true],...
    'fparc_flag',true,[false true],...
    'fnames_fparc',[],[],...
    'motion_flag',true,[false true],...
...
    'outdir','ROI_Summaries',[],...
    'outstem','taskBOLD',[],...
    'verbose',1,[0:2],...
    'forceflag',false,[false true],...
... % hidden
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'statnames',{'beta','tstat'},{'beta','tstat'},...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
    'fname_fscolorlut',[],[],...
    'aseg_roilist',[1:28,40:60],[],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'fparc_roilist',[21101:21300,22101:22300],[1,Inf],...
    'exclude_roilist',[1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59],[],...
... % hidden
    'required_containers',{'proc_bold','fsurf'},[],...
    'QC_BOLD',true,[false true],... % only applies if manual BOLDQC.csv file exists %% NOTE: this is not yet implemented
    'QC_recon',true,[false true],...
...
    'eprime_tags',{'rootdir','indir','fname_info','tasknames',...
                   'outdir','outstem','forceflag'},[],...
    'projinfo_tags',{'VisitIDs','SubjIDs','StudyInfo','RootDirs','ignore_VisitInfo_flag',...
                 'user','numvec_tags'},[],...
    'scan_info_tags',{'TR','numTRs'},{'TR','numTRs','nreps'},...
    'motion_tags',{'mean_motion'},{'mean_motion','mean_trans','mean_rot'},...
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
    'aseg_code_tags',{'aseg_aparc_flag','aseg_roilist','aparc_roilist',...
                      'exclude_roilist'},[],...
    'aparc_code_tags',{'aparc_roilist'},[],...
    'fparc_code_tags',{'fparc_roilist'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~iscell(parms.tasknames)
    parms.tasknames = {parms.tasknames};
  end;
  parms.ntasks = length(parms.tasknames);
  if ~iscell(parms.statnames)
    parms.statnames = {parms.statnames};
  end;
  parms.nstats = length(parms.statnames);
  
  % set runlist, depending on concat_flag and max_nruns
  if parms.concat_flag
    parms.runlist = 0;
    parms.avg_flag = 0;
  else
    parms.runlist = [1:parms.max_nruns];
    if parms.max_nruns < 2
      parms.avg_flag = 0;
    end;
  end;
  
  % check whether anything to do
  if ~parms.aseg_flag && ~parms.aparc_flag
    fprintf('%s: ERROR: aseg_flag = 0 and aparc_flag = 0 (nothing to do)\n',mfilename);
    parms.quitflag = 1;
    return;
  else
    parms.quitflag = 0;
  end;

  % check fnames_fparc
  if parms.fparc_flag && ~isempty(parms.fnames_fparc)
    if ~iscell(parms.fnames_fparc)
      parms.fnames_fparc = {parms.fnames_fparc};
    end;
    for f=1:length(parms.fnames_fparc)
      if ~exist(parms.fnames_fparc{f},'file')
        fprintf('%s: ERROR: fparc annot file %s not found\n',...
          mfilename,parms.fnames_fparc{f});
        parms.quitflag = 1;
        return;
      end;
    end;    
  end;

  % set output directory
  if mmil_isrelative(parms.outstem)
    if mmil_isrelative(parms.outdir)
      parms.outdir = sprintf('%s/MetaData/%s/%s',getenv('HOME'),ProjID,parms.outdir);
    end;
    parms.outstem = [parms.outdir '/' parms.outstem];
  else
    parms.outdir = fileparts(parms.outstem);
  end;
  mmil_mkdir(parms.outdir);
  
  % check output file
  fname_out = sprintf('%s.csv',parms.outstem);
  if exist(fname_out,'file') && ~parms.forceflag
    parms.quitflag = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_roi_parms(parms)
  % set roitypes
  parms.roitypes = {};
  if parms.aseg_flag, parms.roitypes{end+1} = 'aseg'; end;
  if parms.aparc_flag, parms.roitypes{end+1} = 'aparc'; end;
  if parms.fparc_flag, parms.roitypes{end+1} = 'fparc'; end;

  if parms.aseg_flag || parms.aparc_flag
    % read roi codes and names from FreeSurfer color lookup table
    if isempty(parms.fname_fscolorlut)
      MMPS_parms = getenv('MMPS_PARMS');
      parms.fname_fscolorlut = [MMPS_parms '/MMIL_FSColorLUT.txt'];
    end;
    if ~exist(parms.fname_fscolorlut,'file')
      error('FreeSurfer color lookup table file %s not found',...
        parms.fname_fscolorlut);
    end;
    [fs_roicodes,fs_roinames] = fs_colorlut(parms.fname_fscolorlut);
  end;

  % set parameters for aseg analysis
  if parms.aseg_flag
    args = mmil_parms2args(parms,parms.aseg_code_tags);
    parms.aseg_roicodes = mmil_aseg_roicodes(args{:});
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aseg_roicodes);
    parms.aseg_roicodes = fs_roicodes(ind_fs);
    parms.aseg_roinames = fs_roinames(ind_fs);
    parms.aseg_nrois = length(parms.aseg_roicodes);
  end;

  % set parameters for aparc analysis
  if parms.aparc_flag
    args = mmil_parms2args(parms,parms.aparc_code_tags);
    parms.aparc_roicodes = mmil_aparc_roicodes(args{:});
    [tmp,ind_fs]=intersect(fs_roicodes,parms.aparc_roicodes);
    parms.aparc_roicodes = fs_roicodes(ind_fs);
    parms.aparc_roinames = fs_roinames(ind_fs);
    parms.aparc_nrois = length(parms.aparc_roicodes);
  end;

  % set parameters for fparc analysis
  if parms.fparc_flag
    [fp_roinames,fp_roicodes] = get_fparc_roinames(parms);
    args = mmil_parms2args(parms,parms.fparc_code_tags);
    parms.fparc_roicodes = mmil_fparc_roicodes(args{:});
    [tmp,ind_fp]=intersect(fp_roicodes,parms.fparc_roicodes);
    parms.fparc_roicodes = fp_roicodes(ind_fp);
    parms.fparc_roinames = fp_roinames(ind_fp);
    parms.fparc_nrois = length(parms.fparc_roicodes);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roinames,roicodes] = get_fparc_roinames(parms)
  roinames = []; roicodes = [];
  for i=1:length(parms.fnames_fparc)
    fname = parms.fnames_fparc{i};
    
    % decide if fname_aparc is for lh or rh
    [tpath,tstem,text] = fileparts(fname);
    n = regexp(tstem,'(?<hemi>[lr]h)\.(?<name>.+)','names');
    if isempty(n)
      fprintf('%s: WARNING: annot file name %s has unexpected pattern...\n',...
        mfilename,fname);
      roicode_offset = 23100;
      hemi = 'xh';
    else
      hemi = n.hemi;
      switch hemi
        case 'lh'
          roicode_offset = 21100;
        case 'rh'
          roicode_offset = 22100;
      end;
    end;    
    [roinums,roilabels,ctab] = fs_read_annotation(fname);
    tmp_roicodes = roicode_offset + [1:length(roilabels)];
    tmp_roinames = cellfun(@(x) sprintf('ctx-%s-%s',hemi,x),roilabels,...
      'UniformOutput',false)';
    roicodes = cat(2,roicodes,tmp_roicodes);
    roinames = cat(2,roinames,tmp_roinames);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_results(parms)
  % load project and study info
  args = mmil_parms2args(parms,parms.projinfo_tags);
  [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Quick_Check_ProjID(parms.ProjID,args{:});
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
  parms = check_data_files(parms);
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

function summarize_results(parms)
  fname_out = sprintf('%s.csv',parms.outstem);
  if ~exist(fname_out,'file') || parms.forceflag
    % compile ROI averages
    results = compile_results(parms);
    % write summary files
    write_csv(fname_out,results,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_results(parms)
  if parms.verbose==2
    fprintf('%s: compiling analysis results...\n',mfilename);
  end;
  % initialize results
  results = init_results(parms);
  % compile results
  results.data = [];
  results.colnames = [];
  for t=1:parms.ntasks
    taskname = parms.tasknames{t};
    [condnames,nconds] = set_condnames(parms,t);
    % initialize data for averaging values across runs
    if parms.avg_flag
      avg_data = 0;
      avg_data_nruns = 0;
      avg_colnames = [];
    end;
    for k=1:length(parms.runlist)
      run = parms.runlist(k);
      run_data = [];
      run_colnames = [];
      for r=1:length(parms.roitypes)
        roitype = parms.roitypes{r};
        nrois = parms.([roitype '_nrois']);
        nvals = nrois*nconds*parms.nstats;
        roinames = parms.([roitype '_roinames']);
        % initialize data and column names
        tmp_data = nan(parms.nsubs,nvals);
        tmp_colnames = cell(1,nvals);
        % initialize column names for average
        if parms.avg_flag && k==1
          tmp_avg_colnames = cell(1,nvals);
        end;
        for m=1:parms.nstats
          statname = parms.statnames{m};
          for c=1:nconds
            condname = condnames{c};
            for i=1:nrois
              n = (m-1)*nconds*nrois + (c-1)*nrois + i;
              roiname = roinames{i};
              if parms.concat_flag || parms.max_nruns<2
                tmp_colnames{n} = sprintf('%s_%s_%s_%s',...
                                  taskname,condname,statname,roiname);
              else
                tmp_colnames{n} = sprintf('%s_run%d_%s_%s_%s',...
                                  taskname,run,condname,statname,roiname);
                if parms.avg_flag && k==1
                  tmp_avg_colnames{n} = sprintf('%s_%s_%s_%s',...
                                        taskname,condname,statname,roiname);
                end;
              end;
            end;
          end;
        end;        
        % compile results from each subject
        for s=1:parms.nsubs
          for m=1:parms.nstats
            statname = parms.statnames{m};
            roi_data = load_roi_data(parms,s,t,run,roitype,statname);
            if isempty(roi_data), continue; end;
            [tmp,ind_data,ind_sel] = ...
              intersect([roi_data.roicode],parms.([roitype '_roicodes']));
            for c=1:nconds
              for i=1:length(ind_sel)
                n = (m-1)*nconds*nrois + (c-1)*nrois + ind_sel(i);
                tmp_data(s,n) =...
                  roi_data(ind_data(i)).avg(:,c);
              end;
            end;
          end;
        end;
        % concatenate aseg and aparc into all
        run_data = cat(2,run_data,tmp_data);
        run_colnames = cat(2,run_colnames,tmp_colnames);
        if parms.avg_flag && k==1
          avg_colnames = cat(2,avg_colnames,tmp_avg_colnames);
        end;
      end;
      results.data = cat(2,results.data,run_data);
      results.colnames = cat(2,results.colnames,run_colnames);
      if parms.avg_flag
        avg_data_nruns = avg_data_nruns + 1.0*(~isnan(run_data));
        run_data(isnan(run_data)) = 0;
        avg_data = avg_data + run_data;
      end;        
      % compile TR, numTR and motion stats
      results = compile_scan_info(results,parms,t,run);
    end;    
    if parms.avg_flag
      results = combine_scan_info(results,parms,t);
      avg_data = avg_data ./ avg_data_nruns;
      results.data = cat(2,results.data,avg_data);
      results.colnames = cat(2,results.colnames,avg_colnames);
    end;
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results = [];
  results.StudyInfo = parms.StudyInfo;
  results.RootDirs = parms.RootDirs;
  results.SubjIDs = {parms.StudyInfo.SubjID};
  results.VisitIDs = {parms.StudyInfo.VisitID};
  results.STRUCT_VisitIDs = {parms.StudyInfo.STRUCT_VisitID};
  if ~isfield(parms.StudyInfo,'VisitNumber')
    results.VisitNumbers = zeros(parms.nsubs,1);
  else
    results.VisitNumbers = cell2mat({parms.StudyInfo.VisitNumber});
  end;
  % store relevant numbers
  results.nsubs = parms.nsubs;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = load_roi_data(parms,s,t,run,roitype,statname)
  roi_data = [];
  fname = set_data_fname(parms,s,t,run,roitype,statname);
  if isempty(fname) || ~exist(fname,'file'), return; end;
  load(fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compile_scan_info(results,parms,t,run)
  prefix = parms.tasknames{t};
  if run>0 && parms.max_nruns>1
    prefix = sprintf('%s_run%d',prefix,run);
  end;
  results.([prefix '_TR']) = nan(parms.nsubs,1);
  results.([prefix '_nreps']) = nan(parms.nsubs,1);
  results.([prefix '_numTRs']) = nan(parms.nsubs,1);
  if parms.motion_flag
    results.([prefix '_mean_motion']) = nan(parms.nsubs,1);
    results.([prefix '_mean_trans']) = nan(parms.nsubs,1);
    results.([prefix '_mean_rot']) = nan(parms.nsubs,1);
  end;
  for s=1:parms.nsubs
    [TR,nreps,numTRs] = load_subj_scan_info(parms,s,t,run);
    results.([prefix '_TR'])(s) = TR;
    results.([prefix '_nreps'])(s) = nreps; %% NOTE: nreps and numTRs are identical except for APE sequence
    results.([prefix '_numTRs'])(s) = numTRs;
    if parms.motion_flag
      motion_stats = load_subj_motion(parms,s,t,run);
      if ~isempty(motion_stats)
        results.([prefix '_mean_motion'])(s) = motion_stats.mean_motion;
        results.([prefix '_mean_trans'])(s) = motion_stats.mean_trans;
        results.([prefix '_mean_rot'])(s) = motion_stats.mean_rot;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = combine_scan_info(results,parms,t);
  nruns = 0;
  TR = 0;
  nreps = 0;
  numTRs = 0;
  if parms.motion_flag
    mean_motion = 0;
    mean_trans = 0;
    mean_rot = 0;
  end;
  for k=1:length(parms.runlist)
    run = parms.runlist(k);
    prefix = parms.tasknames{t};
    prefix = sprintf('%s_run%d',prefix,run);
    nruns = nruns + 1.0*(results.([prefix '_numTRs']) > 0);
    TR = TR + replace_nans(results.([prefix '_TR']));
    nreps = nreps + replace_nans(results.([prefix '_nreps']));
    numTRs = numTRs +  replace_nans(results.([prefix '_numTRs']));
    if parms.motion_flag
      mean_motion = mean_motion + replace_nans(results.([prefix '_mean_motion']));
      mean_trans = mean_trans + replace_nans(results.([prefix '_mean_trans']));
      mean_rot = mean_rot + replace_nans(results.([prefix '_mean_rot']));
    end;
  end;    
  prefix = parms.tasknames{t};
  results.([prefix '_TR']) =  TR ./ nruns;
  results.([prefix '_nreps']) = nreps;
  results.([prefix '_numTRs']) = numTRs;
  if parms.motion_flag
    results.([prefix '_mean_motion']) = mean_motion ./ nruns;
    results.([prefix '_mean_trans']) = mean_trans ./ nruns;
    results.([prefix '_mean_rot']) = mean_rot ./ nruns;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = replace_nans(vals)
  vals(isnan(vals)) = 0;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_data_files(parms)
  if parms.verbose==2
    fprintf('%s: checking data files...\n',mfilename);
  end;
  valid_flags = zeros(parms.nsubs,parms.ntasks);
  for s=1:parms.nsubs
    VisitID = parms.StudyInfo(s).VisitID;
    for t=1:parms.ntasks
      taskname = parms.tasknames{t};
      for m=1:parms.nstats
        statname = parms.statnames{m};
        for r=1:length(parms.roitypes)
          for k=1:length(parms.runlist)
            run = parms.runlist(k);
            snums = set_snums(parms,s,t,run);
            if isempty(snums)
              if parms.verbose
                if run~=0
                  fprintf('%s: WARNING: no scan for %s %s run %d\n',...
                    mfilename,VisitID,taskname,run);
                else
                  fprintf('%s: WARNING: no scan for %s %s\n',...
                    mfilename,VisitID,taskname);
                end;
              end;
              continue;
            else
              if parms.verbose
                if run~=0
                  fprintf('%s: checking results for %s %s run %d (scan %d)...\n',...
                    mfilename,VisitID,taskname,run,snums);
                else
                  fprintf('%s: checking results for %s %s (scan %s)...\n',...
                    mfilename,VisitID,taskname,strtrim(sprintf('%d ',snums)));
                end;
              end;
            end;
            fname = set_data_fname(parms,s,t,run,parms.roitypes{r},statname);
            if isempty(fname)
              if parms.verbose
                fprintf('%s: WARNING: no data file found for %s, task %s, roitype %s, stat %s\n',...
                  mfilename,VisitID,parms.tasknames{t},parms.roitypes{r},statname);
              end;
              continue;
            end;
            if ~exist(fname,'file')
              if parms.verbose
                fprintf('%s: WARNING: data file %s not found\n',mfilename,fname);
              end;
              continue;
            end;
          end;
          valid_flags(s,t) = 1;
        end;
      end;
    end;
  end;
  % require at least one task to be valid
  ind_valid = find(any(valid_flags,2)); 
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
  valid_flags = zeros(parms.nsubs,parms.ntasks);
  for s=1:parms.nsubs
    VisitID = parms.StudyInfo(s).VisitID;
    for t=1:parms.ntasks
      taskname = parms.tasknames{t};
      for k=1:length(parms.runlist)
        run = parms.runlist(k);
        snums = set_snums(parms,s,t,run);
        if isempty(snums)
          if parms.verbose
            if run~=0
              fprintf('%s: WARNING: no scan for %s %s run %d\n',...
                mfilename,VisitID,taskname,run);
            else
              fprintf('%s: WARNING: no scan for %s %s\n',...
                mfilename,VisitID,taskname);
            end;
          end;
          continue;
        end;
        fname = set_motion_fname(parms,s,t,run);
        if isempty(fname)
          if parms.verbose
            fprintf('%s: WARNING: no motion file found for %s %s\n',...
              mfilename,VisitID,taskname);
          end;
          continue;
        end;
        if ~exist(fname,'file')
          if parms.verbose
            fprintf('%s: WARNING: motion file %s not found for %s %s\n',...
              mfilename,fname,VisitID,taskname);
          end;
          continue;
        end;
        if parms.max_motion<Inf
          motion_stats = load_subj_motion(parms,s,t,run);
          if motion_stats.mean_motion > parms.max_motion
            fprintf('%s: WARNING: excluding visit %s %s with excessive motion...\n',...
              mfilename,VisitID,taskname);
            continue;
          end;    
        end;
        valid_flags(s,t) = 1;
      end;
    end;
  end;
  % require at least one task to be valid
  ind_valid = find(any(valid_flags,2)); 
  if length(ind_valid)<parms.nsubs
    parms.StudyInfo = parms.StudyInfo(ind_valid);
    parms.nsubs = length(parms.StudyInfo);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [condnames,nconds] = set_condnames(parms,t)
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
  condnames = cat(2,stim_labels,{conds_contrast.name});
  nconds = length(condnames);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function snums = set_snums(parms,s,t,run)
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
    return;
  end;
  
  % find corresponding BOLD snums
  ind_BOLD = [ContainerInfo.ScanInfo.BOLD.SeriesIndex];
  [tmp,ind_s,ind_B] = intersect(ind_series,ind_BOLD);
  snums = ind_B;
  if (run > 0) && ~isempty(snums)
    if run > length(snums)
      snums = [];
    else
      snums = snums(run);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [instem,indir] = set_input(parms,s,t,run)
  instem = []; indir = [];
  VisitID = parms.StudyInfo(s).VisitID;
  ContainerPath = MMIL_Get_Container(parms.RootDirs,VisitID,'proc_bold');
  tp = [];
  tp.snums = set_snums(parms,s,t,run);
  if isempty(tp.snums), return; end;
  tp.dirstem = [parms.instem '_' parms.tasknames{t}];
  tp.infix = parms.infix;
  tp.analysis_outfix = parms.analysis_outfix;
  args = mmil_parms2args(tp);
  [instem,errcode] = BOLD_MMIL_Set_GLM_Stem(ContainerPath,args{:});
  indir = fileparts(instem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_data_fname(parms,s,t,run,roitype,statname)
  fname = [];
  instem = set_input(parms,s,t,run);
  if isempty(instem), return; end;
  [indir,instem] = fileparts(instem);
  if strcmp(roitype,'aseg') && parms.aseg_erode_flag
    infix = 'aseg_erode';
  else
    infix = roitype;
  end;
  fname = sprintf('%s/%s_%s_%s_roi_data.mat',...
    indir,instem,statname,infix);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_info_fname(parms,s,t,run)
  fname = [];
  instem = set_input(parms,s,t,run);
  if isempty(instem), return; end;
  instem = regexprep(instem,'_3dDeconv','_info');
  fname = sprintf('%s.mat',instem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_motion_fname(parms,s,t,run)
  fname = [];
  instem = set_input(parms,s,t,run);
  if isempty(instem), return; end;
  instem = regexprep(instem,'_3dDeconv','_motion');
  fname = sprintf('%s.mat',instem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TR,nreps,numTRs] = load_subj_scan_info(parms,s,t,run)
  TR = nan; nreps = nan; numTRs = nan;
  fname = set_info_fname(parms,s,t,run);
  if isempty(fname) || ~exist(fname,'file'), return; end;  
  info = load(fname);
  TR = info.scan_info.TRs;
  nreps = info.scan_info.nreps;
  numTRs = info.scan_info.numTRs;
  if run==0
    TR = mean(TR);
    nreps = sum(nreps);
    numTRs = sum(numTRs);
  elseif length(TR)>1
    TR = TR(run);
    nreps = nreps(run);
    numTRs = numTRs(run);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_stats = load_subj_motion(parms,s,t,run)
  motion_stats = [];
  mean_motion = nan; mean_trans = nan; mean_rot = nan;
  fname = set_motion_fname(parms,s,t,run);
  if isempty(fname) || ~exist(fname,'file'), return; end;
  motion_stats = load(fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_csv(fname,results,parms)
  if parms.verbose==2
    fprintf('%s: writing results to csv %s...\n',mfilename,fname);
  end;
  % initialize data cell array
  data = cell(results.nsubs,0);
  col_labels = cell(1,0);
  % add SubjIDs as first column    
  data = cat(2,data,results.SubjIDs');
  col_labels = cat(2,col_labels,'SubjID');
  % add VisitIDs as second column
  data = cat(2,data,results.VisitIDs');
  col_labels = cat(2,col_labels,'VisitID');
  % add info from StudyInfo
  if parms.continfo_flag || parms.subjinfo_flag
    for t=1:length(parms.info_tags)
      tag = parms.info_tags{t};
      if isfield(results.StudyInfo,tag)
        info = {results.StudyInfo.(tag)}';
        data = cat(2,data,info);
        col_labels = cat(2,col_labels,tag);
      end;
    end;
  end;
  % add scan info and motion stats
  for t=1:length(parms.tasknames)
    taskname = parms.tasknames{t};
    runlist = parms.runlist;
    if parms.avg_flag
      runlist = [runlist 0];
    end;
    for k=1:length(runlist)
      run = runlist(k);
      prefix = parms.tasknames{t};
      if run>0
        prefix = sprintf('%s_run%d',prefix,run);
      end;
      tags = {};
      for i=1:length(parms.scan_info_tags)
        tags{end+1} = [prefix '_' parms.scan_info_tags{i}];
      end;
      if parms.motion_flag
        for i=1:length(parms.motion_tags)
          tags{end+1} = [prefix '_' parms.motion_tags{i}];
        end;
      end;
      for i=1:length(tags)
        tag = tags{i};
        data = cat(2,data,num2cell(results.(tag)));
        col_labels = cat(2,col_labels,tag);
      end;
    end;
  end;
  % add ROI data
  tmp_data = results.data;
  data = cat(2,data,num2cell(tmp_data));
  col_labels = cat(2,col_labels,mmil_rowvec(results.colnames));
  % add column labels
  data = cat(1,col_labels,data);
  % write file
  mmil_write_csv(fname,data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

