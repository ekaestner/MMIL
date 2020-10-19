function MEG_MMIL_Sensor_GroupAvg(ProjID,varargin)
%function MEG_MMIL_Sensor_GroupAvg(ProjID,varargin)
%
% Purpose: Calculate cross-subject average of sensor waveforms
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%    Must have 'MEG_VisitID' to indicate name of the orig_meg data directory
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include this field: proc_meg
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'outdir': output directory
%    full path or relative to /home/{user}/MetaData/{ProjID}
%    if ProjID is empty, will be relative to current working directory
%    {default = 'Sensor_GroupAvg'}
%  'outstem': output file stem (relative to outdir)
%    {default = 'Sensor'}
%  'proc_prefix': prefix of input processed MEG files
%    {default = 'proc'}
%  'proc_infix': infix of input processed MEG files
%     e.g. [], 'combcond', 'subnull'
%     {default = []}
%  'conditions': vector of condition numbers; if empty, use all conditions
%    {default = []}
%  'condition_names': cell array of condition names
%     If empty, use names like condition1, condition2, etc.
%    {default = []}
%  'usegrad_flag': [0|1] use gradiometer data, if available
%    {default = 1}
%  'usemag_flag': [0|1] use magnetometer data, if available
%    {default = 0}
%  'useEEG_flag': [0|1] use EEG data, if available
%    {default = 0}
%  'channames': cell array of channel names (sensors) to average
%    if supplied, ignores usegrad_flag, usemag_flag, and useEEG_flag
%    {default = []}
%  'grad_scalefact': scaling factor applied to gradiometer data
%     {default = 10^13}
%  'mag_scalefact': scaling factor applied to magnetometer data
%     {default = 10^15}
%  'EEG_scalefact': scaling factor applied to EEG data
%     {default = 10^6}
%  'power_flag': [0|1] calculate sensor power, otherwise abs
%    {default = 0}
%  'norm_flag': [0|1|2] normalize sensor waveforms
%     0: no normalization
%     1: normalize within subject waveforms to max of absolute values
%     2: normalize group average waveforms
%     3: normalize within subject and group average waveforms
%    {defaut = 1}
%  'time': vector of time points (msec) used to resample waveforms
%    if empty, will use t0, t1, and samp_dur to create time vector
%    {default = []}
%  't0': starting time point (msec)
%    {default = -100}
%  't1': ending time point (msec)
%    {default = 300}
%  'samp_dur': sample duration (msec) for common time points
%    ignored if time vector supplied
%    {default = 1}
%  'baseline_flag': [0|1] baseline correction of input data
%    {default = 1}
%  'baseline_t0': start time of baseline period (msec)
%    {default = -Inf} (start at beginning of prestimulus period)
%  'baseline_t1'  : end time of baseline period (msec)
%    {default = 0} (end at trigger onset)
%  'qcflag': [0|1] exclude subjects with StudyInfo.QC = 0
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default: 0}
%
% Created:  09/17/13 by Don Hagler
% Last Mod: 11/26/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: make plots with ts_plot_wforms?
%% todo: find peaks and plot responses?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
[RootDirs,StudyInfo,parms] = check_input(ProjID,varargin);

fname_results = [parms.outstem '_GroupAvg_results.mat'];
if ~exist(fname_results,'file') || parms.forceflag
  results = [];
  results.time = parms.time;
  results.avgS = single(0);
  results.stdS = single(0);
  results.contrasts = [];
  results.nconditions = [];
  results.subjdata = [];
  for s=1:length(StudyInfo)
    fprintf('%s: getting sensor results for SubjID %s, MEG_VisitID %s...\n',...
      mfilename,StudyInfo(s).SubjID,StudyInfo(s).MEG_VisitID);
    % get sensor waveforms for this session
    subjdata = load_subjdata(parms,RootDirs,StudyInfo(s));
    if isempty(subjdata), continue; end;
    % check conditions, initialize results fields if necessary
    results = check_conditions(parms,results,subjdata);
    % combine into struct array
    results.subjdata = [results.subjdata subjdata];
    % add to running totals
    results.avgS = results.avgS + subjdata.wforms;
    results.stdS = results.stdS + subjdata.wforms.^2;
  end;
  % calculate average, stdev, and stderr
  results.N = length(results.subjdata);
  if results.N>2
    results.stdS = sqrt((results.N*results.stdS - results.avgS.^2)./...
                    (results.N*(results.N-1)));
    results.stderrS = results.stdS/sqrt(results.N);
  else
    results.stdS = [];
    results.stderrS = [];
  end;
  if results.N>1
    results.avgS = results.avgS/results.N;
  end;
  if results.N==0, return; end;
  % save results to mat file
  save(fname_results,'results');
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RootDirs,StudyInfo,parms] = check_input(ProjID,options)
  parms_filter = {...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'outdir','Sensor_GroupAvg',[],...
    'outstem','Sensor',[],...
    'proc_prefix','proc',[],...
    'proc_infix',[],[],...
    'conditions',[],[],...
    'condition_names',[],[],...
    'usegrad_flag',true,[false true],...
    'usemag_flag',false,[false true],...
    'useEEG_flag',false,[false true],...
    'channames',[],[],...
    'grad_scalefact',10^13,[-Inf Inf],...
    'mag_scalefact',10^15,[-Inf Inf],...
    'EEG_scalefact',10^6,[-Inf Inf],...
    'power_flag',false,[false true],...
    'norm_flag',1,[0:3],...
    'time',[],[],...
    't0',-100,[],...
    't1',300,[],...
    'samp_dur',1,[],...
    'baseline_flag',true,[false true],...
    'baseline_t0',-Inf,[-Inf Inf],...
    'baseline_t1',0,[-Inf Inf],...
    'qcflag',true,[false true],...
    'forceflag',false,[false true],...  
  ...
    'badchans',[],[],...
    'badchanfile',[],[],...
    'required_containers',{'proc_meg'},[],...
  ...
    'chans_tags',{'badchans','badchanfile',...
                  'usegrad_flag','usemag_flag','useEEG_flag','channames'},[],...
    'baseline_tags',{'time','baseline_t0','baseline_t1'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if ~isempty(ProjInfo)
    % For arg names present in both varargin and ProjInfo
    % the varargin values will appear in merged_args
    ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
    merged_args = mmil_merge_args(options,ProjInfo_args);
    % check that parameters fit allowed range, use defaults if not supplied
    parms = mmil_args2parms(merged_args,parms_filter);
  end;

  nsubs  = length(StudyInfo);
  if ~nsubs
    error('no valid subjects in StudyInfo');
  end;
  if isfield(StudyInfo,'GroupAvg')
    for i=1:nsubs
      if isempty(StudyInfo(i).GroupAvg)
        StudyInfo(i).GroupAvg = 0;
      end;
    end;
    ind_GroupAvg = find([StudyInfo.GroupAvg]);
    StudyInfo = StudyInfo(ind_GroupAvg);
    nsubs  = length(StudyInfo);
    if ~nsubs
      error('no valid subjects with StudyInfo.GroupAvg=1');
    end;
  end;
  if isempty(parms.time)
    parms.time = [parms.t0:parms.samp_dur:parms.t1]/1000;
  end;
  % change outstem to be full path
  if mmil_isrelative(parms.outdir)
    if isempty(ProjID)
      parms.outdir = [pwd '/' parms.outdir];
    else
      parms.outdir = [RootDirs.home '/MetaData/' ProjID '/' parms.outdir];
    end;
  end;
  mmil_mkdir(parms.outdir);
  parms.outstem = [parms.outdir '/' parms.outstem];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = load_subjdata(parms,RootDirs,StudyInfo)
  subjdata = [];
  % load ContainerInfo
  ContainerPath = [RootDirs.proc_meg '/' StudyInfo.proc_meg];
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
  if errcode, return; end;
  % set input file name
  if ~isempty(parms.proc_infix)
    infix = ['_' parms.proc_infix];
  else
    nfiles = length(ContainerInfo.input_data_files);
    if nfiles>1
      infix = [];
    else
      infix = '_1';
    end;
  end;
  fname_data = sprintf('%s/matfiles/%s_avg_data%s.mat',...
    ContainerPath,parms.proc_prefix,infix);
  if ~exist(fname_data,'file')
    error('processed MEG file %s not found\n',fname_data);
  end;
  % load data
  fprintf('%s: loading data from %s...\n',mfilename,fname_data);
  load(fname_data);
  % set subjdata fields
  subjdata = StudyInfo;
  subjdata.ContainerPath = ContainerPath;
  subjdata.fname_data = fname_data;
  % get waveforms
  subjdata = get_wform(parms,subjdata,avg_data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = get_wform(parms,subjdata,avg_data)
  % set badchanfile if present
  if ~isempty(parms.badchanfile) && mmil_isrelative(parms.badchanfile)
    parms.badchanfile = [subjdata.ContainerPath '/' parms.badchanfile];
    if ~exist(parms.badchanfile)
      parms.badchanfile = [];
    end;
  end;
  % select valid channels
  args = mmil_parms2args(parms,parms.chans_tags);
  [ind_chans,ind_grad,ind_mag,ind_EEG] = ts_set_chans(avg_data,args{:});
  if isempty(ind_chans), error('no valid channels'); end;
  nchans = length(ind_chans);
  if isempty(parms.conditions)
    parms.conditions = [1:length(avg_data.averages)];
  end;
  nconds = length(parms.conditions);
  subj_time = avg_data.averages(1).time;
  ntpoints = length(parms.time);
  % set scalefacts
  scalefacts = ones(1,1,nchans);
  if ~isempty(ind_grad)
    [tmp,ind] = intersect(ind_chans,ind_grad);
    scalefacts(ind) = parms.grad_scalefact;
  end;
  if ~isempty(ind_mag)
    [tmp,ind] = intersect(ind_chans,ind_mag);
    scalefacts(ind) = parms.mag_scalefact;
  end;
  if ~isempty(ind_EEG)
    [tmp,ind] = intersect(ind_chans,ind_EEG);
    scalefacts(ind) = parms.EEG_scalefact;
  end;
  % get waveforms for selected channels and conditions from avg_data
  wforms = zeros(ntpoints,nconds,nchans);
  for c=1:nconds
    cond = parms.conditions(c);
    tmp_wforms = avg_data.averages(cond).data(ind_chans,:)';
    % resample to common time samples with spline
    if length(parms.time)~=length(subj_time) || any(parms.time~=subj_time)
      tmp_wforms = spline(subj_time,tmp_wforms',parms.time)';
    end;
    wforms(:,c,:) = reshape(tmp_wforms,[ntpoints,1,nchans]);
  end;
  % apply grad_scalefact, mag_scalefact, EEG_scalefact
  wforms = bsxfun(@times,wforms,scalefacts);    
  % calculate power or absolute values
  if parms.power_flag
    wforms = wforms.^2;
  else
    wforms = abs(wforms);
  end;
  % average across channels
  wforms = mean(wforms,3);
  % subtract baseline
  if parms.baseline_flag
    args = mmil_parms2args(parms,parms.baseline_tags);
    wforms = ts_baseline_wforms(wforms,args{:});
  end;
  % normalize waveforms
  if ismember(parms.norm_flag,[1,3])
    max_val = max(abs(wforms(:)));
    wforms = wforms/max_val;
  end;
  subjdata.wforms = wforms;
  subjdata.nconditions = nconds;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = check_conditions(parms,results,subjdata)
  if isempty(results.nconditions)
    results.nconditions = subjdata.nconditions;
  elseif results.nconditions ~= subjdata.nconditions
    error('number of condition levels for %s does not match\n',...
      subjdata.MEG_VisitID);
  end;
  if ~isempty(parms.condition_names)
    if numel(parms.condition_names)~=results.nconditions
      error('number of condition_names (%d) does not match nconditions (%d)',...
        length(parms.condition_names),results.nconditions);
    else
      results.condition_names = ...
        reshape(parms.condition_names,[1,results.nconditions]);
    end;
  else
    results.condition_names = cell(1,results.nconditions);
    for c=1:results.nconditions
      results.condition_names{c} = sprintf('condition %d',c);
    end;
    parms.condition_names = results.condition_names;
  end;
  if isempty(results.contrasts)
    results.contrasts = mmil_getfield(subjdata,...
      'contrasts',1:results.nconditions);
  end;
return;


