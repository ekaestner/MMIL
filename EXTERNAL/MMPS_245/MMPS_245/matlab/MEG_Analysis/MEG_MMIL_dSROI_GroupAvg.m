function MEG_MMIL_dSROI_GroupAvg(ProjID,varargin)
%function MEG_MMIL_dSROI_GroupAvg(ProjID,varargin)
%
% Purpose: Calculate cross-subject average of dSPM ROI waveforms
%   also performs peak detection and creates plots
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
%  'conditions': vector of condition numbers; if empty, use all conditions
%    {default = []}
%  'condition_names': cell array of condition names
%     If empty, use names like condition1, condition2, etc.
%    {default = []}
%  'peaks_flag': [0|1] find peak amplitude and latency, etc.
%     otherwise, just plot waveforms
%     {default = 0}
%  'fit_wforms_flag': [0|1] find peak amplitude and latency through
%     nonlinear fitting of multiple components (see rc_fit_wforms)
%     otherwise, find minima/maxima (see mmil_peakdet)
%    {default = 0}
%  'use_components': vector of component numbers 1 to 7
%     to model if fit_wforms_flag = 1
%    {default = [1 2 3]}
%  'outdir': output directory
%    full path or relative to /home/{user}/MetaData/{ProjID}
%    if ProjID is empty, will be relative to current working directory
%    {default = 'dSPM_ROI_GroupAvg'}
%  'outstem': output file stem (relative to outdir)
%    {default = 'dSPM_ROI'}
%  'qcflag': [0|1] exclude subjects with StudyInfo.QC = 0
%    {default = 1}
%  'peak_mindiff': minimum difference for peak detection
%     {default = 1}
%  'overlay_roi_flag' [0|1] plot with ROIs overlayed on same plot
%    {default: 0}
%  'roi_colors' cell array of colors for each roi
%    {default: {'b','g','r','y','c','m','k'}}
%  'offset': subtract this value from individual subject waveforms
%    {default: 0}
%  'forceflag': [0|1] overwrite existing output
%    {default: 0}
%
% Optional Input from ProjInfo, command line, or StudyInfo
%  'dSROI_outdir': output directory, relative to ContainerPath
%    {default: 'dSPM_ROI_analysis'}
%  'dSROI_outstem': dSPM ROI output file stem
%    {default = 'dSPM_ROI_results'}
%
% Created:  10/27/11 by Don Hagler
% Last Mod: 08/23/16 by Don Hagler
%

%% todo: option roi_names to select which ROIs to include
%% todo: use ts_plot_wforms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
[RootDirs,StudyInfo,parms] = check_input(ProjID,varargin);

fname_results = [parms.outstem '_GroupAvg_results.mat'];
if ~exist(fname_results,'file') || parms.forceflag
  results = [];
  results.time = parms.time;
  results.avgS = single(0);
  results.stdS = single(0);
  results.nrois = [];
  results.contrasts = [];
  results.nconditions = [];
  results.subjdata = [];
  for s=1:length(StudyInfo)
    fprintf('%s: getting dSPM ROI results for SubjID %s, MEG_VisitID %s...\n',...
      mfilename,StudyInfo(s).SubjID,StudyInfo(s).MEG_VisitID);
    % get dSPM ROI results for this session
    subjresults = load_results(parms,RootDirs,StudyInfo(s));
    if isempty(subjresults), continue; end;

    % extract waveforms, find peaks, calculate overall fit and err variance
    subjdata = get_subjdata(parms,subjresults,StudyInfo(s));

    if isempty(results.nrois)
      results.nrois = subjdata.nrois;
      results.roinames = subjdata.roinames;
      parms.roinames = subjdata.roinames;
    elseif results.nrois ~= subjdata.nrois
      fprintf('%s: ERROR: number of ROIs for %s does not match\n',...
        mfilename,StudyInfo(s).MEG_VisitID);
      return;
    end;
    if isempty(results.nconditions)
      results.nconditions = subjdata.nconditions;
      if isempty(parms.condition_values)
        parms.xlim = [0,results.nconditions+1];
      end;
    elseif results.nconditions ~= subjdata.nconditions
      error('number of condition levels for %s does not match\n',...
        StudyInfo(s).MEG_VisitID);
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
  if parms.peaks_flag
    results = get_group_responses(parms,results);
  end;
  if parms.plotflag
    fprintf('%s: plotting results...\n',mfilename);
    % plot average waveforms
    plot_avg_waveforms(parms,results);
    if parms.peaks_flag
      % plot average and multi-subject responses
      plot_all_responses(parms,results);
    end;
  end;
  if parms.peaks_flag
    % write peaks to csv file
    write_peaks(parms,results.subjdata);
    % write average peaks and auc to csv file
    write_analysis(parms,results);
  end;
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
    'dSROI_outdir','dSPM_ROI_analysis',[],...
    'dSROI_outstem','dSPM_ROI_results',[],...
  ...
    'peaks_flag',false,[false true],...
    'fit_wforms_flag',false,[false true],...
    'use_components',[1 2 3],[1,7],...
    'outdir','dSPM_ROI_GroupAvg',[],...
    'outstem','dSPM_ROI',[],...
    'qcflag',true,[false true],...
    'offset',0,[-Inf,Inf],...
    'forceflag',false,[false true],...  
  ... for waveform analysis
    'smooth_sigma',0,[0,1000],...
    'auc_range',[50,300],[],...
    'auc_nbins',1,[1,100],...
    'auc_baseline_flag',true,[false true],...
    'auc_baseline_range',[-100,0],[],...
    'peak_range',[50,300],[],...
    'peak_pol',-1,[-1,1,1],...
    'peak_mindiff',1,[0,Inf],...
    'normflag',0,[0 1 2],...
    'powerflag',2,[0 1 2],...
    'firstflag',false,[false true],...
    'sfreq',1000,[],...
    't0',-100,[],...
    't1',300,[],...
    'time',[],[],...
  ... % for plotting analysis results
    'plotflag',true,[false true],...
    'plotsubjflag',false,[false true],...
    'conditions',[],[],...
    'condition_values',[],[],...
    'condition_names',[],[],...
    'condition_label','Stimulus Condition',[],...
    'eps_flag',false,[false true],...
    'visible_flag',false,[false true],...
    'linewidth',1,[],...
    'min_linewidth',1,[],...
    'max_linewidth',2.5,[],...
    'fontname','Arial',[],...
    'fontsize',12,[],...
    'ylim_auc',[0,50],[],...
    'ylim_peak',[0,20],[],...
    'ylim_latency',[50,300],[],...
    'ylim_wform',[-0.1,20],[],...
    'xlim_wform',[-100,300],[],...
    'xlim',[0,100],[],...
    'units_wform','nA M',[],...
    'logx_flag',false,[false true],...
    'stderrflag',true,[false true],...
    'overlay_roi_flag',false,[false true],...
    'roi_colors',{'b','g','r','y','c','m','k'},[],...
    'legend_loc','SouthEastOutside',[],...
    'label_flag',true,[false true],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
  ... % avg waveform plotting
    'samp_dur',1,[],...
    'err_samp_dur',15,[],...
    'norm_wform_flag',false,[false true],...
    'norm_time0',50,[],...
    'norm_time1',250,[],...
    'errbar_flag',1,[0,1,2],... % 0=no error bars, 1=error bars, 2=shaded interval
    'fill_alpha',0.1,[0 1],...
  ... % used if fit_wforms_flag=1:
    'latency_bounds',[50,700],[],...
    'amplitude_bounds',[0.1,30],[],...
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
    'polarity',[1 1 1 1 1 1 1],[],...
  ...% initial wform fit parameters
    'latency',  [  100    200    300    400    500    600    700],[],...
    'amplitude',[   5      5      5      5      5      5      5],[],...
    'rise_tc',  [   3      3      3      3      3      3      3],[],...
    'fall_tc',  [   6      6      6      6      6      6      6],[],...
  ...
    'tc_bounds_1',[2 4],[],...
    'latency_bounds_2_3',[],[],...
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
  if parms.norm_wform_flag
    parms.units_wform = 'normalized to max amplitude';
  end;
  if isempty(parms.time)
    parms.time = [floor(parms.xlim_wform(1)):parms.samp_dur:ceil(parms.xlim_wform(2))]/1000;
  end;
  parms.peak_latency_bounds = parms.latency_bounds;
  if strcmp(parms.ylim_auc,'none'), parms.ylim_auc = []; end;
  if strcmp(parms.ylim_peak,'none'), parms.ylim_peak = []; end;
  if strcmp(parms.ylim_latency,'none'), parms.ylim_latency = []; end;
  if strcmp(parms.ylim_wform,'none'), parms.ylim_wform = []; end;
  if parms.fit_wforms_flag
    parms.ncomponents = length(parms.polarity);
    tc_bounds = parms.fall_tc_bounds;
    parms.fall_tc_bounds = zeros(parms.nrois,parms.ncomponents,2);
    for r=1:parms.nrois
      for c=1:parms.ncomponents
        if c==1
          parms.fall_tc_bounds(r,c,:) = parms.tc_bounds_1;
        else
          parms.fall_tc_bounds(r,c,:) = tc_bounds;
        end;
      end;
    end;
    parms.latency_bounds = zeros(parms.nrois,parms.ncomponents,2);
    for r=1:parms.nrois
      for c=1:parms.ncomponents
        if ismember(c,[2 3])
          parms.latency_bounds(r,c,:) = parms.latency_bounds_2_3;
        else
          parms.latency_bounds(r,c,:) = parms.peak_latency_bounds;
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
  indir = [RootDirs.proc_meg '/' StudyInfo.proc_meg];
  if ~isfield(StudyInfo,'dSROI_outdir') || isempty(StudyInfo.dSROI_outdir)
    StudyInfo.dSROI_outdir = parms.dSROI_outdir;
  end;
  if ~isfield(StudyInfo,'dSROI_outstem') || isempty(StudyInfo.dSROI_outstem)
    StudyInfo.dSROI_outstem = parms.dSROI_outstem;
  end;
  clear results;
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    fname = sprintf('%s/%s/%s-%s.mat',indir,...\
      StudyInfo.dSROI_outdir, StudyInfo.dSROI_outstem,hemi);
    if ~exist(fname,'file')
      error('file %s not found\n',fname);
    end;
    results(h) = load(fname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = get_subjdata(parms,results,StudyInfo)
  subjdata = StudyInfo;
  subjdata.SubjID = StudyInfo.SubjID;
  subjdata.MEG_VisitID = StudyInfo.MEG_VisitID;
  subjdata.dSROI_outdir = mmil_getfield(StudyInfo,'dSROI_outdir',[]);
  subjdata.dSROI_outstem = mmil_getfield(StudyInfo,'dSROI_outstem',[]);
  % get source waveforms
  subjdata = get_wform(parms,results,subjdata);
  % find peaks
  if parms.peaks_flag
    subjdata = get_peaks(parms,results,subjdata);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = get_wform(parms,results,subjdata)
  subjdata.wform = [];
  subjdata.nrois = 0;
  subjdata.nconditions = 0;
  subjdata.roinames = {};
  subjdata.roihemis = [];
  subjdata.time = results(1).time/1000;
  S_res = [];
  for h=1:length(parms.hemilist)
    nrois = length(results(h).roinames);
    subjdata.nrois = subjdata.nrois + nrois;
    subjdata.roinames = cat(2,subjdata.roinames,results(h).roinames);    
    subjdata.roihemis = cat(2,subjdata.roihemis,h*ones(1,nrois));

    S = permute(results(h).wforms,[3,2,1]);
    subjdata.nconditions = size(S,3);
    if isempty(parms.conditions)
      subjdata.conditions = [1:subjdata.nconditions];
    else
      % maintain sequence of parms.conditions
      [tmp,ind_conds] = intersect(parms.conditions,[1:subjdata.nconditions]);
      subjdata.conditions = parms.conditions(sort(ind_conds));
      if isempty(subjdata.conditions), error('no valid conditions'); end;
      subjdata.nconditions = length(subjdata.conditions);
      S = S(:,:,subjdata.conditions);
    end;
    S = S - parms.offset;
    % resample to common time samples with spline
    S_tmp = zeros(length(parms.time),nrois);
    for c=1:subjdata.nconditions
      %% todo: don't need to loop over ROIs if time is last dim
      %% todo: resample is faster, but need timeseries struct (see fs_resample_frames)
      for r=1:nrois
        S_tmp(:,r,c) = spline(subjdata.time,S(:,r,c),parms.time);
      end;
    end;
    S_res = cat(2,S_res,S_tmp);
  end;
  % normalize to max value across time range
  %% todo: normalize to stdev of baseline?
  if parms.norm_wform_flag
    t0 = parms.norm_time0;
    t1 = parms.norm_time1;
    % find sample closest in time to that desired for baseline and noise
    [tmp,t0] = min(abs(parms.time-parms.norm_time0/1000));
    [tmp,t1] = min(abs(parms.time-parms.norm_time1/1000));
    S_res = S_res / max(mmil_rowvec(S_res(t0:t1,:,:)));
    %% todo: use something like 95th percentile?
  end;
  subjdata.wform = S_res;
  for r=1:subjdata.nrois
    subjdata.roinames{r} = sprintf('%s-%s',subjdata.roinames{r},...
      parms.hemilist{subjdata.roihemis(r)});
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
      for r=1:subjdata.nrois
        [tmp_minima,tmp_maxima] =...
          get_minmax_wfit(parms,subjdata.wfit.results(r,c));
        minima(r,c) = tmp_minima;
        maxima(r,c) = tmp_maxima;
      end;
    end;
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
  tmp_parms.area_names = subjdata.roinames;
  tmp_parms.area_colors = cell(1,subjdata.nrois);
  k = 1;
  for r=1:subjdata.nrois
    if k>length(parms.roi_colors), k = 1; end;
    tmp_parms.area_colors{r} = parms.roi_colors{k};
    k = k + 1;
  end;  
  tmp_parms.plotflag = parms.plotsubjflag;
  args = MMIL_Args(tmp_parms,'rc_analyze_wforms');  
  subjdata.analysis = rc_analyze_wforms(wforms,args{:});
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

function results = get_group_responses(parms,results)
  types = {'amplitude','deflection','latency','auc'};
  for t=1:length(types)
    type = types{t};
    results.(type) = [];
    if strcmp(type,'auc') && parms.auc_nbins>1
      binflag = 1;
      val_size = [results.nrois,results.nconditions,parms.auc_nbins+1,results.N];
      mean_size = [results.nrois,results.nconditions,parms.auc_nbins+1];
    else
      binflag = 0;
      val_size = [results.nrois,results.nconditions,results.N];
      mean_size = [results.nrois,results.nconditions];
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
        for r=1:results.nrois
          if binflag
            vals = results.(type).values(r,:,:,n);
            tmp_size = size(vals);
            if length(tmp_size)==3
              vals = reshape(vals,[tmp_size(2),tmp_size(3)]);
            end;
            if parms.normflag==1
              vals = vals./(ones(size(vals,1),1)*vals(parms.ind_max_cond,:));
            elseif parms.normflag==2
              vals = vals./(ones(size(vals,1),1)*max(abs(vals),[],1));
            end;
            results.(type).values(r,:,:,n) = vals;
          else
            vals = squeeze(results.(type).values(r,:,n));
            if parms.normflag==1
              vals = vals/vals(parms.ind_max_cond);
            elseif parms.normflag==2
              vals = vals/max(abs(vals));
            end;
            results.(type).values(r,:,n) = vals;
          end;
        end;
      end;
    end;
    % calculate mean, standard deviation, and standard error
    results.(type).mean = nan(mean_size);
    results.(type).std = nan(mean_size);
    results.(type).stderr = nan(mean_size);
    results.(type).N = nan(mean_size);
    for r=1:results.nrois
      if binflag
        vals = results.(type).values(r,:,:,:);
        results.(type).mean(r,:,:) = mean(vals,4);
        results.(type).std(r,:,:) = std(vals,0,4);
        results.(type).stderr(r,:,:) = ...
          results.(type).std(r,:,:)/sqrt(results.N);
        results.(type).N(r,:,:) = results.N;
      else
        vals = results.(type).values(r,:,:);
        tmp_size = size(vals);
        if length(tmp_size)==3
          vals = reshape(vals,[tmp_size(2),tmp_size(3)]);
        end;
        ind_good = find(max(isnan(vals),[],1)==0);
        nbad = results.N - length(ind_good);
        if nbad>0
          fprintf('%s: WARNING: missing %s peak for %d/%d subjects\n',...
            mfilename,results.roinames{r},nbad,results.N);
        end;
        if isempty(ind_good), continue; end;
        vals = vals(:,ind_good);
        results.(type).mean(r,:) = mean(vals,2);
        results.(type).std(r,:) = std(vals,0,2);
        results.(type).stderr(r,:) = ...
          results.(type).std(r,:)/sqrt(length(ind_good));
        results.(type).N(r,:) = length(ind_good);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_avg_waveforms(parms,results)
  if parms.overlay_roi_flag
    figure; clf; hold on;
    % plot waveform
    for c=1:results.nconditions
      k = 1;
      for r=1:results.nrois
        if k>length(parms.roi_colors), k=1; end;
        plot(1000*parms.time,results.avgS(:,r,c),...
          parms.roi_colors{k},'LineWidth',parms.linewidths(c));
        k = k + 1;
      end;
    end;
    % plot error bars
    if parms.errbar_flag>0 && ~isempty(results.stdS)
      for c=1:results.nconditions
        k = 1;
        for r=1:results.nrois
          if k>length(parms.roi_colors), k=1; end;
          if parms.errbar_flag==1
            tmp_time = 1000*parms.time(1:parms.err_samp_dur:end);
            tmp_vals = results.avgS(1:parms.err_samp_dur:end,r,c);
            tmp_err  = results.stderrS(1:parms.err_samp_dur:end,r,c);
            errorbar(tmp_time,tmp_vals,tmp_err,...
              [parms.roi_colors{k} '.'],'LineWidth',parms.linewidth);
          elseif parms.errbar_flag==2
            % shaded confidence interval      
            tmp_time = 1000*[parms.time fliplr(parms.time)];
            tmp_vals = mmil_rowvec(results.avgS(:,r,c));
            tmp_err  = mmil_rowvec(results.stderrS(:,r,c));
            tmp_lo  = tmp_vals - tmp_err;
            tmp_hi  = tmp_vals + tmp_err;
            tmp_intvl = [tmp_hi fliplr(tmp_lo)];
            h = fill(tmp_time,tmp_intvl,parms.roi_colors{k});            
            set(h,'FaceAlpha',parms.fill_alpha,'EdgeAlpha',parms.fill_alpha);
          end;
          k = k + 1;
        end;
      end;
    end;
    if parms.label_flag
      ylabl = sprintf('Source Amplitude (%s)',parms.units_wform);
      toplabl = 'Average Source Waveforms';
      plot_waveform_labels(parms,toplabl,ylabl);
    end;
    save_plot(parms,[parms.outstem '_avgwf']);
  else % plot separate figures for each roi
    for r=1:results.nrois
      figure; clf; hold on;
      % plot waveform
      k = 1;
      for c=1:results.nconditions
        if k>length(parms.roi_colors), k=1; end;
        plot(1000*parms.time,results.avgS(:,r,c),...
          parms.roi_colors{k},'LineWidth',parms.linewidths(c));
        k = k + 1;
      end;
      % plot error bars
      k = 1;
      for c=1:results.nconditions
        if ~isempty(results.stdS) && parms.errbar_flag>0
          if k>length(parms.roi_colors), k=1; end;
          if parms.errbar_flag==1
            tmp_vals = results.avgS(1:parms.err_samp_dur:end,r,c);
            tmp_err  = results.stderrS(1:parms.err_samp_dur:end,r,c);
            tmp_time = 1000*parms.time(1:parms.err_samp_dur:end);
            errorbar(tmp_time,tmp_vals,tmp_err,...
              [parms.roi_colors{k} '.'],'LineWidth',parms.linewidth);
          elseif parms.errbar_flag==2
            % shaded confidence interval      
            tmp_time = 1000*[parms.time fliplr(parms.time)];
            tmp_vals = mmil_rowvec(results.avgS(:,r,c));
            tmp_err  = mmil_rowvec(results.stderrS(:,r,c));
            tmp_lo  = tmp_vals - tmp_err;
            tmp_hi  = tmp_vals + tmp_err;
            tmp_intvl = [tmp_hi fliplr(tmp_lo)];
            h = fill(tmp_time,tmp_intvl,parms.roi_colors{k});            
            set(h,'FaceAlpha',parms.fill_alpha,'EdgeAlpha',parms.fill_alpha);
          end;
          k = k + 1;
        end;
      end;
      if parms.label_flag
        ylabl = sprintf('Source Amplitude (%s)',parms.units_wform);
        toplabl = 'Average Source Waveforms';
        plot_waveform_labels(parms,toplabl,ylabl);
      end;
      save_plot(parms,[parms.outstem '_' results.roinames{r} '_avgwf']);
    end; 
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_waveform_labels(parms,toplabl,ylabl)
  set(gca,'XLim',1000*[parms.time(1),parms.time(end)]);
  if ~isempty(parms.ylim_wform), set(gca,'YLim',parms.ylim_wform); end;
  set(gca,parms.fontargs{:});
  xlabel('time (msec)',parms.fontargs{:});
  ylabel(ylabl,parms.fontargs{:});
  title(toplabl,parms.fontargs{:});
  if parms.overlay_roi_flag
    legend(parms.roinames,parms.fontargs{:},'Location',parms.legend_loc);
  else
    legend(parms.condition_names,parms.fontargs{:},'Location',parms.legend_loc);
  end;
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
    k = 1;
    for r=1:results.nrois
      if k>length(parms.roi_colors), k = 1; end;
      if strcmp(type,'auc') && parms.auc_nbins>1
        resps = squeeze(results.(type).values(r,:,parms.auc_bin,n));
      else
        resps = squeeze(results.(type).values(r,:,n));
      end;
      % plot amplitude vs. conditions
      if parms.linewidth>0
        format_str = [parms.roi_colors{k} 'o-'];
        plot(parms.condition_values,resps,format_str,...
          'LineWidth',parms.linewidth);
      else
        format_str = [parms.roi_colors{k} 'o'];
        plot(parms.condition_values,resps,format_str);
      end;
      k = k + 1;
    end;
  end;
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
  k = 1;
  for r=1:results.nrois
    if k>length(parms.roi_colors), k = 1; end;
    if strcmp(type,'auc') && parms.auc_nbins>1
      mean_resps = squeeze(results.(type).mean(r,:,parms.auc_bin));
      if parms.stderrflag
        std_resps = squeeze(results.(type).stderr(r,:,parms.auc_bin));
      else
        std_resps = squeeze(results.(type).std(r,:,parms.auc_bin));
      end;
    else
      mean_resps = results.(type).mean(r,:);
      if parms.stderrflag
        std_resps = results.(type).stderr(r,:);
      else
        std_resps = results.(type).std(r,:);
      end;
    end;
    % plot amplitude vs. conditions
    if parms.linewidth>0
      format_str = [parms.roi_colors{k} 'o-'];
      errorbar(parms.condition_values,mean_resps,std_resps,...
        format_str,'LineWidth',parms.linewidth);
    else
      format_str = [parms.roi_colors{k} 'o'];
      errorbar(parms.condition_values,mean_resps,std_resps,format_str);
    end;
    k = k + 1;
  end;
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
  legend(parms.roinames,parms.fontargs{:},'Location',parms.legend_loc);
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
    fprintf(fid,'"SubjID""MEG_VisitID","area","condition","type","latency","amplitude","deflection"\n');
    for N=1:length(subjdata)
      SubjID = subjdata(N).SubjID;
      MEG_VisitID = subjdata(N).MEG_VisitID;
      tmp_subjdata = subjdata(N);
      for c=1:tmp_subjdata.nconditions
        for r=1:tmp_subjdata.nrois
          extypes = {'minimum','maximum'};
          for j=1:length(extypes)
            extype = extypes{j};
            if strcmp(extype,'minimum')
              minmax = tmp_subjdata.analysis.minima(r,c);
            else
              minmax = tmp_subjdata.analysis.maxima(r,c);
            end;          
            % list minima/maxima
            npeaks = length(minmax.latency);
            for i=1:npeaks
              latency = minmax.latency(i);
              amplitude = minmax.amplitude(i);
              deflection = minmax.deflection(i);
              fprintf(fid,'"%s","%s","%s",%d,"%s",%0.2f,%0.2f,%0.2f',...
                SubjID,MEG_VisitID,tmp_subjdata.roinames{r},c,...
                extype,latency,amplitude,deflection);
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
    for r=1:results.nrois
      for c=1:results.nconditions
         fprintf(fid,'"%s",%0.2f',...
           results.roinames{r},parms.condition_values(c));
        for t=1:length(types)
          for s=1:length(stats)
            val = results.(types{t}).(stats{s})(r,c);
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


