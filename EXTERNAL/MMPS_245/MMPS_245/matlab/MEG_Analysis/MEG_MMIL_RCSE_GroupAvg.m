function MEG_MMIL_RCSE_GroupAvg(ProjID,varargin)
%function MEG_MMIL_RCSE_GroupAvg(ProjID,varargin)
%
% Purpose: Calculate cross-subject average of RCSE source estimates
%   also performs peak detection and creates plots
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional project-specific input:
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
%
% Optional Input from ProjInfo, command line, or StudyInfo:
%  'RCSE_prefix': RCSE prefix
%    {default = 'RCSE'}
%  'RCSE_infix': extra string attached to RCSE output prefix
%    {default = []}
%  'RCSE_infix_flag': whether infix was automatically generated
%    ignored if RCSE_infix is not empty
%    {default = 0}
%  'RCSE_fstem_conds' - stem of csv file containing condition information
%    {default = 'cond_info'}
%
% Optional Parameters:
%  'diff_prefix': prefix of RCSE results to be subtracted
%    If empty, no subtraction will be performed
%    {default = []}
%  'diff_infix': infix for RCSE results to be subtracted
%    {default = []}
%  'fit_wforms_flag': [0|1] find peak amplitude and latency through
%     nonlinear fitting of multiple components (see rc_fit_wforms)
%     otherwise, find minima/maxima (see mmil_peakdet)
%  'use_components': vector of component numbers 1 to 7
%     to model if fit_wforms_flag = 1
%     {default = [1 2 3]}
%  'hiC_flag': [0|1] use source results for highest condition only
%     if 0, all subjects must have same number of condition levels
%    {default = 0}
%  'outdir': output directory
%    full path or relative to /home/{user}/MetaData/{ProjID}
%    if ProjID is empty, will be relative to current working directory
%    {default = 'RCSE_GroupAvg'}
%  'outstem': output file stem (relative to outdir)
%    {default = 'RCSE'}
%  'qcflag': [0|1] exclude subjects with StudyInfo.QC = 0
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  03/05/11 by Don Hagler
% Last Mod: 04/07/14 by Don Hagler
%

%% todo: move waveform fitting to separate function that uses aggregated data

%% todo: create a function to accept matrix of waveforms and do
%%   peak detection, averaging, plotting etc.
%%   without needing to know anything about StudyInfo, RootDirs, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

[RootDirs,StudyInfo,parms] = check_input(ProjID,varargin);

fname_results = [parms.outstem '_GroupAvg_results.mat'];
if ~exist(fname_results,'file') || parms.forceflag
  results = [];
  results.time = parms.time;
  results.avgS = single(0);
  results.stdS = single(0);
  results.nareas = parms.nareas;
  results.contrasts = [];
  results.nconditions = [];
  results.subjdata = [];
  for s=1:length(StudyInfo)
    fprintf('%s: getting RCSE results for SubjID %s, MEG_VisitID %s...\n',...
      mfilename,StudyInfo(s).SubjID,StudyInfo(s).MEG_VisitID);
    % get RCSE results for this session
    subjresults = load_results(parms,RootDirs,StudyInfo(s));
    if isempty(subjresults), continue; end;

    % extract waveforms, find peaks, calculate overall fit and err variance
    subjdata = get_subjdata(parms,subjresults,StudyInfo(s));

    if isempty(results.nconditions)
      results.nconditions = subjdata.nconditions;
    elseif results.nconditions ~= subjdata.nconditions
      error('number of condition levels for %s does not match - try hiC_flag = 1\n',...
        StudyInfo(s).MEG_VisitID);
    end;
    if isempty(results.contrasts)
      results.contrasts = mmil_getfield(subjresults,...
        'contrasts',1:results.nconditions);
    end;

    % vary linewidth for varying condition levels
    if results.nconditions > 1
      parms.linewidths = parms.min_linewidth + ...
                          ([1:results.nconditions]-1)*(parms.max_linewidth - parms.min_linewidth)/...
                          (results.nconditions-1);
    else
      parms.linewidths = parms.linewidth;
    end;
    % combine into struct array
    results.subjdata = [results.subjdata subjdata];
    % add to running totals
    results.avgS = results.avgS + subjdata.wform;
    results.stdS = results.stdS + subjdata.wform.^2;
  end;

  if isempty(parms.condition_values)
    parms.condition_values = unique(results.contrasts);
  end;
  if length(parms.condition_values) ~= results.nconditions
    error('number of condition_values (%d) does not match number of conditions (%d)\n',...
      length(parms.condition_values),results.nconditions);
  end;
  [tmp,parms.ind_max_cond] = max(parms.condition_values);

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
  results = get_group_responses(parms,results);
  if parms.plotflag
    % plot average waveforms
    plot_avg_waveforms(parms,results);
    % plot average and multi-subject responses
    plot_all_responses(parms,results);
  end;
  % write peaks to csv file
  write_peaks(parms,results.subjdata);
  % write best fit and err to csv file
  write_fitvar(parms,results.subjdata);
  % write average peaks and auc to csv file
  write_analysis(parms,results);
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
  ...
    'RCSE_prefix','RCSE',[],...
    'RCSE_infix',[],[],...
    'RCSE_infix_flag',false,[false true],...
    'RCSE_fstem_conds','cond_info',[],...
  ...
    'diff_prefix',[],[],...
    'diff_infix',[],[],...
    'fit_wforms_flag',false,[false true],...
    'use_components',[1 2 3],[1,7],...
    'hiC_flag',false,[false true],...
    'outdir','RCSE_GroupAvg',[],...
    'outstem','RCSE',[],...
    'qcflag',true,[false true],...
    'forceflag',false,[false true],...  
  ... for waveform analysis
    'smooth_sigma',0,[0,1000],...
    'auc_range',[0,350],[],...
    'auc_nbins',1,[1,100],...
    'auc_baseline_flag',true,[false true],...
    'auc_baseline_range',[-100,0],[],...
    'peak_range',[50,150],[],...
    'peak_pol',-1,[-1,1,1],...
    'peak_mindiff',0.5,[0,Inf],...
    'normflag',1,[0 1 2],...
    'powerflag',2,[0 1 2],...
    'firstflag',false,[false true],...
    'sfreq',1000,[],...
    't0',-100,[],...
    't1',350,[],...
    'time',[],[],...
  ... % for plotting analysis results
    'plotflag',true,[false true],...
    'plot_subj_flag',false,[false true],...
    'outdir',pwd,[],...
    'outstem',[],[],...
    'condition_values',[],[],...
    'condition_label','Stimulus Condition',[],...
    'eps_flag',false,[false true],...
    'visible_flag',false,[false true],...
    'linewidth',1,[],...
    'min_linewidth',1,[],...
    'max_linewidth',2.5,[],...
    'axes_linewidth',1,[],...
    'fontname','Arial',[],...
    'fontsize',12,[],...
    'ylim_auc',[0,1.3],[],...
    'ylim_peak',[0,1.3],[],...
    'ylim_latency',[50,150],[],...
    'ylim_wform',[-25,10],[],...
    'xlim',[0,1.05],[],...
    'units_wform','nA M',[],...
    'logx_flag',false,[false true],...
    'stderrflag',true,[false true],...
    'area_names',{'V1','V2','V3'},[],...
    'area_colors',{'b','g','r'},[],...
    'legend_loc','SouthEast',[],...
    'label_flag',true,[false true],...
  ... % avg waveform plotting
    'err_samp_dur',15,[],...
    'norm_wform_flag',false,[false true],...
    'errbar_flag',1,[0,1,2],... % 0=no error bars, 1=error bars, 2=shaded interval
    'fill_alpha',0.1,[0 1],...
  ... % used if fit_wforms_flag=1:
    'latency_bounds',[50,300],[],...
    'amplitude_bounds',[0.1,25],[],...
    'rise_tc_bounds',[2,6],[],...
    'fall_tc_bounds',[4,10],[],...
    'search_type','fmincon',[],...
    'rand_init_flag',1,[],...
    'delay_sf',0,[],...
    'niters',500,[],...
    'stepsize',0.05,[],...
    'linewidth_data',1,[],...
    'linewidth_fit',2,[],...
    'use_areas',[],[],...
    'contrast_latency_flag',0,[],...
    'visible_flag',0,[],...
    'polarity',[1 -1 -1 1 -1 1 -1],[],...
  ...% initial wform fit parameters
  ...%   (chosen to fit group-contrained solution for high contrast stimuli)
  ...%              P1     N1a    N1b    P1     N2     P2     N3
    'latency',  [   50     75     90    145    188    214    264;         % V1
                    65     92    107    167    200    220    299;         % V2
                    60     88    103    150    180    213    256 ],[],... % V3
    'amplitude',[  0.1      6      4      2      2      2      2;
                   0.2      3      3      2    0.1      2    0.3;
                   0.2      3      3      1    0.1      2    0.8 ],[],...
    'rise_tc',  [    3      3      3      3      3      3      3;
                     3      3      3      3      3      3      3;
                     3      3      3      3      3      3      3 ],[],...
    'fall_tc',  [    3      6      6      6      6      6      6;
                     3      6      6      6      6      6      6;
                     3      6      6      6      6      6      6],[],...
  ...
    'tc_bounds_1',[2 4],[],...
    'latency_bounds_2_3',[70 150],[],...
  ...
    'required_containers',{'proc_meg'},[],...
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
  if ~isfield(StudyInfo,'RCSE_fstem_conds') ||...
     isempty(StudyInfo.RCSE_fstem_conds)
    for i=1:nsubs
      StudyInfo(i).RCSE_fstem_conds = parms.RCSE_fstem_conds;
    end;
  end;
  if ~isfield(StudyInfo,'RCSE_infix') ||...
     isempty(StudyInfo.RCSE_infix)
    for i=1:nsubs
      StudyInfo(i).RCSE_infix = parms.RCSE_infix;
    end;
  end;
  parms.nareas = length(parms.area_names);
  if parms.norm_wform_flag
    parms.units_wform = 'normalized to V1 peak';
  end;
  if isempty(parms.time)
    samp_dur = 1000/parms.sfreq;
    parms.time = [parms.t0:samp_dur:parms.t1]/1000;
  end;
  parms.peak_latency_bounds = parms.latency_bounds;
  if strcmp(parms.ylim_auc,'none'), parms.ylim_auc = []; end;
  if strcmp(parms.ylim_peak,'none'), parms.ylim_peak = []; end;
  if strcmp(parms.ylim_latency,'none'), parms.ylim_latency = []; end;
  if strcmp(parms.ylim_wform,'none'), parms.ylim_wform = []; end;
  if parms.fit_wforms_flag
    parms.ncomponents = length(parms.polarity);
    tc_bounds = parms.fall_tc_bounds;
    parms.fall_tc_bounds = zeros(parms.nareas,parms.ncomponents,2);
    for a=1:parms.nareas
      for c=1:parms.ncomponents
        if c==1
          parms.fall_tc_bounds(a,c,:) = parms.tc_bounds_1;
        else
          parms.fall_tc_bounds(a,c,:) = tc_bounds;
        end;
      end;
    end;
    parms.latency_bounds = zeros(parms.nareas,parms.ncomponents,2);
    for a=1:parms.nareas
      for c=1:parms.ncomponents
        if ismember(c,[2 3])
          parms.latency_bounds(a,c,:) = parms.latency_bounds_2_3;
        else
          parms.latency_bounds(a,c,:) = parms.peak_latency_bounds;
        end;
      end;
    end;
    if isempty(parms.use_components)
      parms.use_components = [1:parms.ncomponents];
    end;
    parms.use_components = intersect(parms.use_components,[1:parms.ncomponents]);
    if length(parms.use_components) ~= parms.ncomponents
      parms.polarity = parms.polarity(parms.use_components);  
      parms.latency = parms.latency(:,parms.use_components);
      parms.amplitude = parms.amplitude(:,parms.use_components);
      parms.rise_tc = parms.rise_tc(:,parms.use_components);
      parms.fall_tc = parms.fall_tc(:,parms.use_components);
      parms.fall_tc_bounds = parms.fall_tc_bounds(:,parms.use_components,:);
      parms.latency_bounds = parms.latency_bounds(:,parms.use_components,:);
      parms.ncomponents = length(parms.use_components);
    end;
  end;
  parms.fontargs = {'FontSize',parms.fontsize,'FontName',parms.fontname};
  % set time bins for area under curve
  if parms.auc_nbins>1
    parms.auc_bins = zeros(parms.auc_nbins+1,2);
    range_dur = range(parms.auc_range);
    bin_dur = range_dur / parms.auc_nbins;
    for n=1:parms.auc_nbins
      parms.auc_bins(n,1) = parms.auc_range(1) + (n-1)*bin_dur;
      parms.auc_bins(n,2) = parms.auc_bins(n,1) + bin_dur;
    end;
    parms.auc_bins(parms.auc_nbins+1,:) = parms.auc_range;
  else
    parms.auc_bins = parms.auc_range;
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

function results = load_results(parms,RootDirs,StudyInfo)
  results = [];
  indir = [RootDirs.proc_meg '/' StudyInfo.proc_meg];
  infix = StudyInfo.RCSE_infix;
  if parms.RCSE_infix_flag && isempty(infix)
    infix = rc_RCSE_set_infix(parms);
  end;
  if ~isfield(StudyInfo,'RCSE_prefix') || isempty(StudyInfo.RCSE_prefix)
    StudyInfo.RCSE_prefix = parms.RCSE_prefix;
  end;
  if ~isempty(infix)
    StudyInfo.RCSE_prefix = [StudyInfo.RCSE_prefix '_' infix];
  end;
  fname = [indir '/matfiles/' StudyInfo.RCSE_prefix '_results.mat'];
  if ~exist(fname,'file')
    error('file %s not found\n',fname);
  end;
  load(fname);

  if ~isempty(parms.diff_prefix)
    fname = [indir '/matfiles/'...
             parms.diff_prefix '_' parms.diff_infix '_results.mat'];
    if ~exist(fname,'file')
      error('file %s not found\n',fname);
    end;
    tmp = load(fname);
    results.S = results.S - tmp.results.S;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = get_subjdata(parms,results,StudyInfo)
  subjdata = StudyInfo;
  subjdata.SubjID = StudyInfo.SubjID;
  subjdata.MEG_VisitID = StudyInfo.MEG_VisitID;
  subjdata.RCSE_prefix = mmil_getfield(StudyInfo,'RCSE_prefix',[]);
  subjdata.RCSE_fstem_conds = mmil_getfield(StudyInfo,'RCSE_fstem_conds',[]);
  % get source waveforms
  subjdata = get_wform(parms,results,subjdata);
  % get best fit variance and residual error
  subjdata = get_fitvar(parms,results,subjdata);
  % find peaks
  subjdata = get_peaks(parms,results,subjdata);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = get_wform(parms,results,subjdata)
  S = results.S;
  subjdata.nareas = parms.nareas;
  subjdata.nconditions = results.ncontrasts;
  if parms.hiC_flag && subjdata.nconditions>1
    % if multiple condition levels (or durations), use last one only
    S = S(:,:,end);
    subjdata.nconditions = 1;
  end;
  % resample to common time samples with spline
  S_res = zeros(length(parms.time),parms.nareas,subjdata.nconditions);
  for c=1:subjdata.nconditions
    for a=1:subjdata.nareas
      S_res(:,a,c) = spline(results.time,S(:,a,c),parms.time);
    end;
  end;
  % normalize to last condition value for V1
  if parms.norm_wform_flag
    S_res = S_res / max(abs(S_res(:,1,end)));
  end;
  subjdata.wform = S_res;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = get_fitvar(parms,results,subjdata)
  if ~isempty(results.Ygrad)
    var_Ygrad = var(results.Ygrad,0,2);
    var_Yfitgrad = var(results.Yfitgrad,0,2);
    var_Egrad = var(results.Egrad,0,2);
    max_var_Ygrad = max(var_Ygrad);
    var_Ygrad = var_Ygrad/max_var_Ygrad;
    var_Yfitgrad = var_Yfitgrad/max_var_Ygrad;
    var_Egrad = var_Egrad/max_var_Ygrad;
    subjdata.var_grad = spline(results.time,var_Ygrad,parms.time);
    subjdata.fit_grad = spline(results.time,var_Yfitgrad,parms.time);
    subjdata.err_grad = spline(results.time,var_Egrad,parms.time);
    [subjdata.fit_grad_max,ind] = max(var_Yfitgrad);
    subjdata.err_grad_max = var_Egrad(ind);
  else
    subjdata.var_grad = [];
    subjdata.fit_grad = [];
    subjdata.err_grad = [];
    subjdata.fit_grad_max = [];
    subjdata.err_grad_max = [];
  end;

  if ~isempty(results.Ymag)
    var_Ymag = var(results.Ymag,0,2);
    var_Yfitmag = var(results.Yfitmag,0,2);
    var_Emag = var(results.Emag,0,2);
    max_var_Ymag = max(var_Ymag);
    var_Ymag = var_Ymag/max_var_Ymag;
    var_Yfitmag = var_Yfitmag/max_var_Ymag;
    var_Emag = var_Emag/max_var_Ymag;
    subjdata.var_mag = spline(results.time,var_Ymag,parms.time);
    subjdata.fit_mag = spline(results.time,var_Yfitmag,parms.time);
    subjdata.err_mag = spline(results.time,var_Emag,parms.time);
    [subjdata.fit_mag_max,ind] = max(var_Yfitmag);
    subjdata.err_mag_max = var_Emag(ind);
  else
    subjdata.var_mag = [];
    subjdata.fit_mag = [];
    subjdata.err_mag = [];
    subjdata.fit_mag_max = [];
    subjdata.err_mag_max = [];
  end;

  if ~isempty(results.Yeeg)
    var_Yeeg = var(results.Yeeg,0,2);
    var_Yfiteeg = var(results.Yfiteeg,0,2);
    var_Eeeg = var(results.Eeeg,0,2);
    max_var_Yeeg = max(var_Yeeg);
    var_Yeeg = var_Yeeg/max_var_Yeeg;
    var_Yfiteeg = var_Yfiteeg/max_var_Yeeg;
    var_Eeeg = var_Eeeg/max_var_Yeeg;
    subjdata.var_eeg = spline(results.time,var_Yeeg,parms.time);
    subjdata.fit_eeg = spline(results.time,var_Yfiteeg,parms.time);
    subjdata.err_eeg = spline(results.time,var_Eeeg,parms.time);
    [subjdata.fit_eeg_max,ind] = max(var_Yfiteeg);
    subjdata.err_eeg_max = var_Eeeg(ind);
  else
    subjdata.var_eeg = [];
    subjdata.fit_eeg = [];
    subjdata.err_eeg = [];
    subjdata.fit_eeg_max = [];
    subjdata.err_eeg_max = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = get_peaks(parms,results,subjdata)
  % fit waveform components
  subjdata.wfit = [];
  if parms.fit_wforms_flag
    tmp_parms = parms;
    tmp_parms.outstem = sprintf('%s_%s_%s',...
      parms.outstem,subjdata.SubjID,subjdata.MEG_VisitID);
    tmp_parms.time = 1000*parms.time;
    tmp_parms.xlim = 1000*[parms.time(1),parms.time(end)];
    args = MMIL_Args(tmp_parms,'rc_fit_wforms');
    [subjdata.wfit.results,subjdata.wfit.wforms] = ...
      rc_fit_wforms(subjdata.wform,args{:});
    % create minmia and maxima structs
    clear minima maxima;
    for c=1:subjdata.nconditions
      for a=1:subjdata.nareas
        [tmp_minima,tmp_maxima] = get_minmax_wfit(parms,subjdata.wfit.results(a,c));
        minima(a,c) = tmp_minima;
        maxima(a,c) = tmp_maxima;
      end;
    end;
    minima = get_minmax_fit(minima,results);
    maxima = get_minmax_fit(maxima,results);
    subjdata.wfit.minima = minima;
    subjdata.wfit.maxima = maxima;
    wforms = subjdata.wfit.wforms;
  else
    wforms = subjdata.wform;
  end;

  % find peaks, calculate auc, make plots
  tmp_parms = parms;
  tmp_parms.outstem = sprintf('%s_%s_%s',...
    parms.outstem,subjdata.SubjID,subjdata.MEG_VisitID);
  tmp_parms.time = 1000*parms.time;
  tmp_parms.plotflag = parms.plot_subj_flag;
  args = MMIL_Args(tmp_parms,'rc_analyze_wforms');
  subjdata.analysis = rc_analyze_wforms(wforms,args{:});
  subjdata.analysis.minima = get_minmax_fit(subjdata.analysis.minima,results);
  subjdata.analysis.maxima = get_minmax_fit(subjdata.analysis.maxima,results);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [minima,maxima] = get_minmax_wfit(parms,wfit_results)
  minima = struct('latency',[],'amplitude',[],'deflection',[]);
  maxima = minima;
  j = 1; k = 1;
  for i=1:wfit_results.ncomponents
    latency = wfit_results.latency(i);
    amplitude = wfit_results.amplitude(i);
    rise_tc = wfit_results.rise_tc(i);
    fall_tc = wfit_results.fall_tc(i);
    if wfit_results.polarity(i)<0
      minima.latency(j) = latency;
      minima.amplitude(j) = -amplitude;
      minima.rise_tcs(j) = rise_tc;
      minima.fall_tcs(j) = fall_tc;
      j = j + 1;
    else
      maxima.latency(k) = latency;
      maxima.amplitude(k) = amplitude;
      maxima.rise_tcs(k) = rise_tc;
      maxima.fall_tcs(k) = fall_tc;
      k = k + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax = get_minmax_fit(minmax,results)
  for a=1:size(minmax,1)
    for c=1:size(minmax,2)
      for j=1:length(minmax(a,c).latency)
        % get goodness of fit
        [tmp,ind] = min(abs(minmax(a,c).latency(j)-results.time*1000));
        if ~isempty(results.Ygrad)
          var_Ygrad = var(results.Ygrad,0,2);
          var_Yfitgrad = var(results.Yfitgrad,0,2);
          minmax(a,c).fit_grad(j) = var_Yfitgrad(ind);
          minmax(a,c).var_grad(j) = var_Ygrad(ind);
        else
          minmax(a,c).fit_grad(j) = NaN;
          minmax(a,c).var_grad(j) = NaN;
        end;
        if ~isempty(results.Ymag)
          var_Ymag = var(results.Ymag,0,2);
          var_Yfitmag = var(results.Yfitmag,0,2);
          minmax(a,c).fit_mag(j) = var_Yfitmag(ind);
          minmax(a,c).var_mag(j) = var_Ymag(ind);
        else
          minmax(a,c).fit_mag(j) = NaN;
          minmax(a,c).var_mag(j) = NaN;
        end;
        if ~isempty(results.Yeeg)
          var_Yeeg = var(results.Yeeg,0,2);
          var_Yfiteeg = var(results.Yfiteeg,0,2);
          minmax(a,c).fit_eeg(j) = var_Yfiteeg(ind);
          minmax(a,c).var_eeg(j) = var_Yeeg(ind);
        else
          minmax(a,c).fit_eeg(j) = NaN;
          minmax(a,c).var_eeg(j) = NaN;
        end;
        % get mstart stdev
        if isfield(results,'S_stdv')
          minmax(a,c).mstart_stdv(j) = mean(mmil_rowvec(results.S_stdv(ind,:,:)));
        else
          minmax(a,c).mstart_stdv(j) = NaN;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = get_group_responses(parms,results)
  types = {'amplitude','deflection','latency','auc'};
  for t=1:length(types)
    type = types{t};
    results.(type) = [];
    if strcmp(type,'auc') && parms.auc_nbins>1
      binflag = 1;
      val_size = [results.nareas,results.nconditions,parms.auc_nbins+1,results.N];
      mean_size = [results.nareas,results.nconditions,parms.auc_nbins+1];
    else
      binflag = 0;
      val_size = [results.nareas,results.nconditions,results.N];
      mean_size = [results.nareas,results.nconditions];
    end;
    % compile results for all subjects
    results.(type).values = nan(val_size);
    for n=1:results.N
      analysis = results.subjdata(n).analysis;
      if binflag
        results.(type).values(:,:,:,n) = analysis.(type);
      else
        results.(type).values(:,:,n) = analysis.(type);    
      end;
    end;
    % normalize by condition with max value
    if parms.normflag && ~strcmp(type,'latency')
      results.(type).orig_values = results.(type).values;
      for n=1:results.N
        for a=1:results.nareas
          if binflag
            vals = results.(type).values(a,:,:,n);
            tmp_size = size(vals);
            vals = reshape(vals,[tmp_size(2),tmp_size(3)]);
            if parms.normflag==1
              vals = vals./(ones(size(vals,1),1)*vals(parms.ind_max_cond,:));
            elseif parms.normflag==2
              vals = vals./(ones(size(vals,1),1)*max(abs(vals),[],1));
            end;
            results.(type).values(a,:,:,n) = vals;
          else
            vals = squeeze(results.(type).values(a,:,n));
            if parms.normflag==1
              vals = vals/vals(parms.ind_max_cond);
            elseif parms.normflag==2
              vals = vals/max(abs(vals));
            end;
            results.(type).values(a,:,n) = vals;
          end;
        end;
      end;
    end;
    % calculate mean, standard deviation, and standard error
    results.(type).mean = nan(mean_size);
    results.(type).std = nan(mean_size);
    results.(type).stderr = nan(mean_size);
    results.(type).N = nan(mean_size);
    for a=1:results.nareas
      if binflag
        vals = results.(type).values(a,:,:,:);
        results.(type).mean(a,:,:) = mean(vals,4);
        results.(type).std(a,:,:) = std(vals,0,4);
        results.(type).stderr(a,:,:) = ...
          results.(type).std(a,:,:)/sqrt(results.N);
        results.(type).N(a,:,:) = results.N;
      else
        vals = results.(type).values(a,:,:);
        if numel(size(vals))==3
          vals = reshape(vals,[size(vals,2),size(vals,3)]);
        else
          vals = reshape(vals,[size(vals,1),size(vals,2)]);
        end;        
        ind_good = find(max(isnan(vals),[],1)==0);
        nbad = results.N - length(ind_good);
        if nbad>0
          fprintf('%s: WARNING: missing %s peak for %d/%d subjects\n',...
            mfilename,parms.area_names{a},nbad,results.N);
        end;
        if isempty(ind_good), continue; end;
        vals = vals(:,ind_good);
        results.(type).mean(a,:) = mean(vals,2);
        results.(type).std(a,:) = std(vals,0,2);
        results.(type).stderr(a,:) = ...
          results.(type).std(a,:)/sqrt(length(ind_good));
        results.(type).N(a,:) = length(ind_good);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_avg_waveforms(parms,results)
  figure; clf; hold on;
  % plot waveform
  for c=1:results.nconditions
    for a=1:results.nareas
      plot(1000*parms.time,results.avgS(:,a,c),...
        'Color',parms.area_colors{a},'LineWidth',parms.linewidths(c));
    end;
  end;
  % plot error bars
  for c=1:results.nconditions
    for a=1:results.nareas
      if ~isempty(results.stdS) && parms.errbar_flag>0
        if parms.errbar_flag==1
          tmp_time = 1000*parms.time(1:parms.err_samp_dur:end);
          tmp_vals = results.avgS(1:parms.err_samp_dur:end,a,c);
          tmp_err  = results.stderrS(1:parms.err_samp_dur:end,a,c);
          errorbar(tmp_time,tmp_vals,tmp_err,'.',...
            'Color',parms.area_colors{a},'LineWidth',parms.linewidth);
        elseif parms.errbar_flag==2
          % shaded confidence interval
          tmp_time = 1000*[parms.time fliplr(parms.time)];
          tmp_vals = mmil_rowvec(results.avgS(:,a,c));
          tmp_err  = mmil_rowvec(results.stderrS(:,a,c));
          tmp_lo  = tmp_vals - tmp_err;
          tmp_hi  = tmp_vals + tmp_err;
          tmp_int = [tmp_hi fliplr(tmp_lo)];
          h = fill(tmp_time,tmp_int,parms.area_colors{a});            
          set(h,'FaceAlpha',parms.fill_alpha,'EdgeAlpha',parms.fill_alpha);
        end;
      end;
    end;
  end;
  set(gca,'LineWidth',parms.axes_linewidth);
  if parms.label_flag
    ylabl = sprintf('Source Amplitude (%s)',parms.units_wform);
    toplabl = 'Average Source Waveforms';
    plot_waveform_labels(parms,toplabl,ylabl);
  end;
  save_plot(parms,[parms.outstem '_groupavg_wforms']);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_waveform_labels(parms,toplabl,ylabl)
  set(gca,'XLim',1000*[parms.time(1),parms.time(end)]);
  if ~isempty(parms.ylim_wform), set(gca,'YLim',parms.ylim_wform); end;
  set(gca,parms.fontargs{:});
  xlabel('time (msec)',parms.fontargs{:});
  ylabel(ylabl,parms.fontargs{:});
  title(toplabl,parms.fontargs{:});
  legend(parms.area_names,parms.fontargs{:},'Location',parms.legend_loc);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_all_responses(parms,results)
  types = {'amplitude','deflection','latency','auc'};
  ylabls = {'Peak Amplitude','Peak Deflection',...
            'Peak Latency (msec)','Area Under Curve'};
  toplabls = {...
    ['Peak Amplitude vs. ' parms.condition_label]...
    ['Peak Deflection vs. ' parms.condition_label]...
    ['Peak Latency vs. ' parms.condition_label]...
    ['Area Under Curve vs. ' parms.condition_label]...
  };
  for t=1:length(types)
    type = types{t};
    if strcmp(type,'auc') && parms.auc_nbins>1
      for n=1:parms.auc_nbins+1
        if n>parms.auc_nbins
          outfix = [];
        else
          outfix = sprintf('bin%d',n);
        end;
        toplabl = sprintf('%s   time range %d (%0.0f to %0.0f msec)',...
          toplabls{t},n,parms.auc_bins(n,1),parms.auc_bins(n,2));
        parms.auc_bin = n;
        plot_multisubj_responses(parms,results,type,...
          toplabls{t},ylabls{t},outfix);
        plot_groupavg_responses(parms,results,type,...
          toplabls{t},ylabls{t},outfix);
      end;
    else
      parms.auc_bin = 0;
      plot_multisubj_responses(parms,results,type,toplabls{t},ylabls{t});
      plot_groupavg_responses(parms,results,type,toplabls{t},ylabls{t});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multisubj_responses(parms,results,type,toplabl,ylabl,outfix)
  if ~exist('outfix','var'), outfix = []; end;
  if strcmp(type,'latency'), parms.normflag = 0; end;  
  figure; clf; hold on;
  for n=1:results.N
    for a=1:results.nareas
      if strcmp(type,'auc') && parms.auc_nbins>1
        resps = squeeze(results.(type).values(a,:,parms.auc_bin,n));
      else
        resps = squeeze(results.(type).values(a,:,n));
      end;
      % plot amplitude vs. conditions
      if parms.linewidth>0
        plot(parms.condition_values,resps,'o-',...
          'Color',parms.area_colors{a},'LineWidth',parms.linewidth);
      else
        plot(parms.condition_values,resps,'o','Color',parms.area_colors{a});
      end;
    end;
  end;
  set(gca,'LineWidth',parms.axes_linewidth);
  if parms.label_flag
    plot_response_labels(parms,type,toplabl,ylabl)
  end;
  outstem = [parms.outstem '_multisubj_' type];
  if ~isempty(outfix), outstem = [outstem '_' outfix]; end;
  if parms.normflag==1
    outstem = [outstem '_norm'];
  elseif parms.normflag==2
    outstem = [outstem '_normmax'];
  end;
  save_plot(parms,outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_groupavg_responses(parms,results,type,toplabl,ylabl,outfix)
  if ~exist('outfix','var'), outfix = []; end;
  if results.N<2, return; end;
  if strcmp(type,'latency'), parms.normflag = 0; end;  
  figure; clf; hold on;
  for a=1:results.nareas
    if strcmp(type,'auc') && parms.auc_nbins>1
      mean_resps = squeeze(results.(type).mean(a,:,parms.auc_bin));
      if parms.stderrflag
        std_resps = squeeze(results.(type).stderr(a,:,parms.auc_bin));
      else
        std_resps = squeeze(results.(type).std(a,:,parms.auc_bin));
      end;
    else
      mean_resps = results.(type).mean(a,:);
      if parms.stderrflag
        std_resps = results.(type).stderr(a,:);
      else
        std_resps = results.(type).std(a,:);
      end;
    end;
    % plot amplitude vs. conditions
    if parms.linewidth>0
      errorbar(parms.condition_values,mean_resps,std_resps,...
        'o-','Color',parms.area_colors{a},'LineWidth',parms.linewidth);
    else
      errorbar(parms.condition_values,mean_resps,std_resps,'o','Color',parms.area_colors{a});
    end;
  end;
  set(gca,'LineWidth',parms.axes_linewidth);
  plot_response_labels(parms,type,toplabl,ylabl)
  outstem = [parms.outstem '_groupavg_' type];
  if ~isempty(outfix), outstem = [outstem '_' outfix]; end;
  if parms.normflag==1
    outstem = [outstem '_norm'];
  elseif parms.normflag==2
    outstem = [outstem '_normmax'];
  end;
  save_plot(parms,outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_response_labels(parms,type,toplabl,ylabl)
  switch type
    case 'latency'
      tmp_ylim = parms.ylim_latency;
    case {'amplitude','deflection'}
      tmp_ylim = parms.ylim_peak;
    case 'auc'
      tmp_ylim = parms.ylim_auc;
    otherwise
      error('invalid type: %s',type);
  end;
  if ~isempty(tmp_ylim), set(gca,'YLim',tmp_ylim); end;
  set(gca,'XLim',parms.xlim);
  set(gca,parms.fontargs{:});
  xlabel(parms.condition_label,parms.fontargs{:});
  if strcmp(type,'latency')
    ylabel(ylabl',parms.fontargs{:});
  elseif parms.normflag
    ylabl = ['Relative ' ylabl];
  else
    ylabl = [ylabl ' (nA M)'];
  end;
  ylabel(ylabl,parms.fontargs{:});
  title(toplabl,parms.fontargs{:});
  if parms.logx_flag
    tmp = get(gca);
    tmp_ticks = tmp.XTickLabel;
    for i=1:size(tmp_ticks,1)
      tmp_str = tmp_ticks(i,:);
      tmp_val = str2num(tmp_str);
      tmp_val = 10^tmp_val;
      tmp_str = num2str(tmp_val,'%0.2f');
      tmp_ticks(i,1:length(tmp_str)) = tmp_str;
    end;
    set(gca,'XTickLabel',tmp_ticks);
  end;
  legend(parms.area_names,parms.fontargs{:},'Location',parms.legend_loc);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_plot(parms,outstem)
  if ~parms.visible_flag, set(gcf,'Visible','off'); end;
  print(gcf,'-dtiff',[outstem '.tif']);
  if parms.eps_flag, mmil_printeps(gcf,[outstem '.eps']); end;
  if ~parms.visible_flag, close(gcf); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_peaks(parms,subjdata)
  fname_out = [parms.outstem '_peaks.csv'];
  fid = fopen(fname_out,'wt');
  if fid==-1
    fprintf('%s: WARNNG: failed to open file %s for writing',fname_out);
  else
    fprintf(fid,'"SubjID","MEG_VisitID","RCSE_fstem_conds","area","condition","type","latency","amplitude","deflection",');
    fprintf(fid,'"var_grad","fit_grad","var_mag","fit_mag","var_eeg","fit_eeg",');
    fprintf(fid,'"mstart_stdv"\n');
    for N=1:length(subjdata)
      SubjID = subjdata(N).SubjID;
      MEG_VisitID = subjdata(N).MEG_VisitID;
      RCSE_fstem_conds = subjdata(N).RCSE_fstem_conds;
      tmp_subjdata = subjdata(N);
      for c=1:tmp_subjdata.nconditions
        for a=1:tmp_subjdata.nareas
          extypes = {'minimum','maximum'};
          for j=1:length(extypes)
            extype = extypes{j};
            if strcmp(extype,'minimum')
              minmax = tmp_subjdata.analysis.minima(a,c);
            else
              minmax = tmp_subjdata.analysis.maxima(a,c);
            end;          
            % list minima/maxima
            npeaks = length(minmax.latency);
            for i=1:npeaks
              latency = minmax.latency(i);
              amplitude = minmax.amplitude(i);
              deflection = minmax.deflection(i);
              fit_grad = minmax.fit_grad(i);
              var_grad = minmax.var_grad(i);
              fit_mag = minmax.fit_mag(i);
              var_mag = minmax.var_mag(i);
              fit_eeg = minmax.fit_eeg(i);
              var_eeg = minmax.var_eeg(i);
              mstart_stdv = minmax.mstart_stdv(i);
              fprintf(fid,'"%s","%s","%s","%s",%d,"%s",%0.2f,%0.2f,%0.2f',...
                SubjID,MEG_VisitID,RCSE_fstem_conds,parms.area_names{a},c,...
                extype,latency,amplitude,deflection);
              fprintf(fid,',%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f',...
                var_grad,fit_grad,var_mag,fit_mag,var_eeg,fit_eeg);
              fprintf(fid,',%0.2f',mstart_stdv);
              fprintf(fid,'\n');
            end;
          end;
        end;
      end;
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_analysis(parms,results)
  types = {'amplitude','deflection','latency','auc'};
  stats = {'mean','std','stderr','N'};
  % write values for each area and condition
  fname_out = [parms.outstem '_analysis.csv'];
  fid = fopen(fname_out,'wt');
  if fid==-1
    fprintf('%s: WARNNG: failed to open file %s for writing',fname_out);
  else
    fprintf(fid,'"area","condition"');
    for t=1:length(types)
      for s=1:length(stats)
        fprintf(fid,',"%s %s"',types{t},stats{s});
      end;
    end;
    fprintf(fid,'\n');
    for a=1:results.nareas
      for c=1:results.nconditions
         fprintf(fid,'"%s",%0.2f',...
           parms.area_names{a},parms.condition_values(c));
        for t=1:length(types)
          for s=1:length(stats)
            val = results.(types{t}).(stats{s})(a,c);
            fprintf(fid,',');
            mmil_write_value(fid,val);
          end;
        end;
        fprintf(fid,'\n');
      end;
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_fitvar(parms,subjdata)
  fname_out = [parms.outstem '_fitvar.csv'];
  fid = fopen(fname_out,'wt');
  if fid==-1
    fprintf('%s: WARNNG: failed to open file %s for writing',fname_out);
  else
    fprintf(fid,'"SubjID","MEG_VisitID",RCSE_fstem_conds,"fit grad","err grad","fit mag","err mag","fit eeg","err eeg"\n');
    for N=1:length(subjdata)
      fprintf(fid,'"%s","%s","%s"',...
        subjdata(N).SubjID,subjdata(N).MEG_VisitID,subjdata(N).RCSE_fstem_conds);
      if ~isempty(subjdata(N).fit_grad)
        fprintf(fid,',%0.2f,%0.2f',...
          subjdata(N).fit_grad,subjdata(N).err_grad_max);
      else
        fprintf(fid,',,');
      end;
      if ~isempty(subjdata(N).fit_mag)
        fprintf(fid,',%0.2f,%0.2f',...
          subjdata(N).fit_mag,subjdata(N).err_mag_max);
      else
        fprintf(fid,',,');
      end;
      if ~isempty(subjdata(N).fit_eeg)
        fprintf(fid,',%0.2f,%0.2f',...
          subjdata(N).fit_eeg,subjdata(N).err_eeg_max);
      else
        fprintf(fid,',,');
      end;
      fprintf(fid,'\n');
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

