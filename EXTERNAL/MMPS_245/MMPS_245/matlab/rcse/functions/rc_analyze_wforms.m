function results = rc_analyze_wforms(wforms,varargin)
%function results = rc_analyze_wforms(wforms,[options])
%
% Purpose: analyze properties of source estimate waveforms
%   including peak amplitude and latency and area under curve
%
% Required Input:
%   wforms: 3D waveform vector
%     size of wforms should be [ntpoints,nareas,nconditions]
%
% Optional Parameters:
%   'smooth_sigma': Gaussian blurring kernel applied to waveforms
%     {default = 0}
%   'auc_range': time range (msec) within which to calculate area under curve
%     {default = [50,300]}
%   'auc_nbins': number of bins in which to subdivide auc_range
%     {default = 1}
%   'auc_baseline_flag': whether to subtract baseline area under curve
%     {default = 1}
%   'auc_baseline_range': time range (msec) to calculate baseline area under curve
%     {default = [-100,0]}
%   'peak_range': time range (msec) within which to select peaks
%     {default = [50,150]}
%   'peak_pol': peak polarity (-1 or 1)
%     {default = -1}
%   'peak_mindiff': minimum difference for peak detection
%     {default = 0.5}
%   'firstflag': use first peak in time range, otherwise use largest deflection
%     {default = 0}
%   'normflag': normalize peak amplitude and auc
%     0: raw values
%     1: normalized to value for max condition value
%     2: normalized to max value
%     {default = 1}
%   'powerflag': for area under curve
%     0: raw values      (signed magnitude)
%     1: absolute values (unsigned magnitude)
%     2: squared values  (power)
%     {default = 2}
%   'contrast_latency_flag': [0|1] require longer latency for lower condition values
%     {default = 0}
%   'latency_allowance': number of milliseconds to allow lower condition values
%      to have earlier latency (if contrast_latency_flag = 1)
%     {default = 5}
%   'sfreq': sampling frequency; used to calculate time vector
%     {default = 1000}
%   't0': start time of waveform (msec)
%     {default = -100}
%   't1': end time of waveform (msec)
%     {default = 300}
%   'time': time vector (msec)
%     if supplied, sfreq, t0, and t1 are ignored
%     length of time vector must match length of input wform
%     {default = []}
%
% Optional Parameters for Plotting
%   'plotflag': [0|1] create and save plots
%     {default = 1}
%   'outdir': output directory for plots
%     {default = pwd}
%   'outstem': output file stem
%     {default = []}
%   'condition_values': vector of x-axis values for each condition
%     number must equal nconditions; if empty, will use 1:nconditions
%     {default = []}
%   'condition_label': x-axis label
%     {default 'Stimulus Condition'}
%   'eps_flag': save plot in eps format (vector graphics)
%     {default = 0}
%   'visible_flag': [0|1] whether to display plots to screen
%     (otherwise save only)
%     {default = 0}
%   'linewidth': plot line width
%     {default = 1}
%   'min_linewidth': waveform line width for first condition
%     {default = 1}
%   'max_linewidth': waveform line width for last condition
%     {default = 2.5}
%   'fontname': font name for labels and axes
%     {default = 'Arial'}
%   'fontsize': font size for labels and axes
%     {default = 12}
%   'ylim_auc': y-axis limits for area under curve plots
%     {default = [0,1.3]}
%   'ylim_peak': y-axis limits for peak amplitude and deflection plots
%     {default = [0,1.3]}
%   'ylim_latency': y-axis limits for latency plots
%     {default = [50,150]}
%   'ylim_wform': y-axis limits for waveform plots
%     {default = [-25,10]}
%   'xlim': x-axis limits (condition values)
%     {default = [0,1.05]}
%   'xlim_wform': x-axis limits for waveform plots
%     if empty, use entire time range
%     {default = []}
%   'units_wform': unit of measurement for waveforms
%     {default = 'nA M'}
%   'logx_flag': [0|1] use log scale for x-axis
%     {default = 0}
%   'area_names': cell array of visual area names
%     {default = {'V1','V2','V3'}}
%   'area_colors': cell array of color codes for each area
%     {default = {'b','g','r'}
%   'legend_loc': location of area name legend
%     {default = 'SouthEast'}
%   'label_flag': [0|1] display labels on plots
%     {default = 1}
%   'mark_peaks_flag': [0|1] mark peaks on waveform plots
%     {default = 1}
%
% Created:  07/20/11 by Don Hagler
% Last Mod: 04/07/14 by Don Hagler
%

%% todo: use ts_wform_peaks, ts_wform_auc, etc.
%%       and ts_plot_wforms and ts_plot_resposnes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'smooth_sigma',0,[0,1000],...
  'auc_range',[50,300],[],...
  'auc_nbins',1,[1,100],...
  'auc_baseline_flag',true,[false true],...
  'auc_baseline_range',[-100,0],[],...
  'peak_range',[50,150],[],...
  'peak_pol',-1,[-1,1,1],...
  'peak_mindiff',0.5,[0,Inf],...
  'normflag',1,[0 1 2],...
  'powerflag',2,[0 1 2],...
  'firstflag',false,[false true],...
  'contrast_latency_flag',false,[false true],...
  'latency_allowance',5,[0,100],...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'time',[],[],...
... % for plotting
  'plotflag',true,[false true],...
  'outdir',pwd,[],...
  'outstem',[],[],...
  'condition_values',[],[],...
  'condition_label','Stimulus Condition',[],...
  'eps_flag',false,[false true],...
  'visible_flag',false,[false true],...
  'linewidth',1,[],...
  'min_linewidth',1,[],...
  'max_linewidth',2.5,[],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
  'ylim_auc',[0,1.3],[],...
  'ylim_peak',[0,1.3],[],...
  'ylim_latency',[50,150],[],...
  'ylim_wform',[-25,10],[],...
  'xlim',[0,1.05],[],...
  'xlim_wform',[],[],...
  'units_wform','nA M',[],...
  'logx_flag',false,[false true],...
  'area_names',{'V1','V2','V3'},[],...
  'area_colors',{'b','g','r'},[],...
  'legend_loc','SouthEast',[],...
  'label_flag',true,[false true],...
  'mark_peaks_flag',true,[false true],...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(wforms,varargin,parms_filter);
results = init_results(parms);
if parms.smooth_sigma
  wforms = smooth_waveforms(wforms,parms);
end;
results = get_auc(wforms,parms,results);
results = get_peaks(wforms,parms,results);
if parms.plotflag
  plot_all_responses(parms,results);
  plot_waveforms(wforms,parms,results);
  plot_cond_waveforms(wforms,parms,results);
end;
write_analysis(parms,results);
write_peaks(parms,results);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(wforms,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  parms.ntpoints = size(wforms,1);
  parms.nareas = size(wforms,2);
  parms.nconditions = size(wforms,3);
  parms.fontargs = {'FontSize',parms.fontsize,'FontName',parms.fontname};
  if length(parms.peak_range)~=2, error('peak range must have two elements'); end;
  if ~isempty(parms.outstem), parms.outstem = [parms.outstem '_']; end;
  mmil_mkdir(parms.outdir);
  if mmil_isrelative(parms.outstem)
    parms.outstem = [parms.outdir '/' parms.outstem];
  end;
  if strcmp(parms.ylim_auc,'none'), parms.ylim_auc = []; end;
  if strcmp(parms.ylim_peak,'none'), parms.ylim_peak = []; end;
  if strcmp(parms.ylim_latency,'none'), parms.ylim_latency = []; end;
  if strcmp(parms.ylim_wform,'none'), parms.ylim_wform = []; end;
  if isempty(parms.condition_values)
    parms.condition_values = 1:parms.nconditions;
  end;
  [tmp,parms.ind_max_cond] = max(parms.condition_values);
  if parms.logx_flag
    parms.condition_values = log10(parms.condition_values);
    parms.xlim = [min(parms.condition_values),0];
  end;
  if isempty(parms.time)
    sd = 1000/parms.sfreq;
    parms.time = [parms.t0:sd:parms.t1];
  end;
  if length(parms.time) ~= parms.ntpoints
    error('number of elements in time (%d) does not match wform (%d)',...
      length(parms.time),parms.ntpoints);
  end;
  if isempty(parms.xlim_wform)
    parms.xlim_wform = [parms.time(1),parms.time(end)];
  end;
  % vary linewidth for varying condition values
  if parms.nconditions > 1
    parms.linewidths = parms.min_linewidth + ...
      ([1:parms.nconditions]-1)*(parms.max_linewidth - parms.min_linewidth)/...
      (parms.nconditions-1);
  else
    parms.linewidths = parms.linewidth;
  end;
  if parms.auc_nbins>1
    parms.auc_bins = zeros(parms.auc_nbins+1,2); % extra bin for entire range
    range_dur = range(parms.auc_range);
    bin_dur = range_dur / parms.auc_nbins;
    for n=1:parms.auc_nbins
      parms.auc_bins(n,1) = parms.auc_range(1) + (n-1)*bin_dur;
      parms.auc_bins(n,2) = parms.auc_bins(n,1) + bin_dur;
    end;
    parms.auc_nbins = parms.auc_nbins + 1; % extra bin for entire range
    parms.auc_bins(parms.auc_nbins,:) = parms.auc_range;
  else
    parms.auc_bins = parms.auc_range;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results = [];
  fnames = {'nareas','nconditions','time'};
  for f=1:length(fnames)
    results.(fnames{f}) = parms.(fnames{f});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wforms = smooth_waveforms(wforms,parms)
  for a=1:parms.nareas
    for c=1:parms.nconditions
      wforms(:,a,c) = mmil_smooth(squeeze(wforms(:,a,c)),parms.smooth_sigma);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = get_auc(wforms,parms,results)
  auc = nan(results.nareas,results.nconditions,parms.auc_nbins);
  switch parms.powerflag
    case 1
      wforms = abs(wforms);
    case 2
      wforms = wforms.^2;
  end;
  for n=1:parms.auc_nbins
    [tmp,ind0] = min(abs(results.time - parms.auc_bins(n,1)));
    [tmp,ind1] = min(abs(results.time - parms.auc_bins(n,2)));
    auc_range = [ind0:ind1];
    n_auc_range = length(auc_range);
    if parms.auc_baseline_flag
      [tmp,b_ind0] = min(abs(results.time - parms.auc_baseline_range(1)));
      [tmp,b_ind1] = min(abs(results.time - parms.auc_baseline_range(2)));
      auc_b_range = [b_ind0:b_ind1];
    end;
    j = 0;
    for c=1:results.nconditions
      for a=1:results.nareas
        % area under curve (integrate over time range)
        if parms.auc_baseline_flag
          % average baseline across conditions
          auc_b = mean(mmil_rowvec(wforms(auc_b_range,a,:)));
        else
          auc_b = 0;
        end;
        auc(a,c,n) = mean(squeeze(wforms(auc_range,a,c))) - auc_b;
      end;
    end;
  end;
  results.auc = auc;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = get_peaks(wforms,parms,results)
  % find minima and maxima
  clear minima maxima;
  for c=1:results.nconditions
    for a=1:results.nareas
      % select time range of wform
      tmp_wform = wforms(:,a,c);
      [maxtab,mintab] = mmil_peakdet(tmp_wform,parms.peak_mindiff,results.time);
      tmp_minima = get_minmax(parms,mintab);
      tmp_maxima = get_minmax(parms,maxtab);
      tmp_minima = get_deflection(tmp_minima,tmp_maxima);
      tmp_maxima = get_deflection(tmp_maxima,tmp_minima);
      minima(a,c) = tmp_minima;
      maxima(a,c) = tmp_maxima;
    end;
  end;
  results.minima = minima;
  results.maxima = maxima;

  % get peak amplitude, deflection, and latency
  amplitude = nan(results.nareas,results.nconditions);
  deflection = amplitude;
  latency = amplitude;
  % peaks
  if parms.peak_pol>0
    extrema = results.maxima;
  else
    extrema = results.minima;
  end;

  if parms.contrast_latency_flag
    [tmp,cond_order] = sort(mmil_rowvec(parms.condition_values),2,'descend');
  else
    cond_order = [1:results.nconditions];
  end;

  for cond=1:results.nconditions
    c = cond_order(cond);
    for a=1:results.nareas
      % peak amplitude, deflection, and latency
      peaks = extrema(a,c); % for one waveform (area x contrast)
      npeaks = length(peaks.latency);
      % select one peak
      if npeaks
        % require that latency increases with decreasing condition value
        ind = 1:npeaks;
        if parms.contrast_latency_flag && cond>1
          min_latency = latency(a,cond_order(cond-1)) - parms.latency_allowance;
          if ~isnan(min_latency)
            ind = find(peaks.latency>=min_latency);
          end;
        end;
        npeaks = length(ind);
      end;
      if ~npeaks
        fprintf('%s: missing peak for condition %0.2f, area %s\n',...
          mfilename,parms.condition_values(c),parms.area_names{a});
        continue;
      end;
      % select first or largest peak
      if parms.firstflag || npeaks==1
        ind = ind(1);
      else
        [tmp,tmp_ind] = max(parms.peak_pol*peaks.amplitude(ind));
        ind = ind(tmp_ind);
      end;

      deflection(a,c) = peaks.deflection(ind);
      amplitude(a,c) = peaks.amplitude(ind);
      latency(a,c) = peaks.latency(ind);
    end;
  end;
  results.amplitude = amplitude;
  results.deflection = deflection;
  results.latency = latency;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax = get_minmax(parms,mtab)
  minmax = struct('latency',[],'amplitude',[],'deflection',[]);
  j = 1;
  for i=1:size(mtab,1)
    latency = mtab(i,1);
    amplitude = mtab(i,2);
    if latency<parms.peak_range(1) ||...
       latency>parms.peak_range(2)
      continue;
    end;
    minmax.latency(j) = latency;
    minmax.amplitude(j) = amplitude;
    j = j + 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax1 = get_deflection(minmax1,minmax2)
  % calculate deflection for minima or maxima
  for i=1:length(minmax1.latency)
    latency = minmax1.latency(i);
    % find preceding minimum
    tmp = latency - minmax2.latency;
    tmp(tmp<0)=Inf;
    [tmp,ind1] = min(tmp);
    if isempty(tmp) || isinf(tmp)
      amp1 = 0;
    else
      amp1 = minmax2.amplitude(ind1);
    end;

    % find subsequent maxima or minima
    tmp = minmax2.latency - latency;
    tmp(tmp<0)=Inf;
    [tmp,ind2] = min(tmp);
    if isempty(tmp) || isinf(tmp)
      amp2 = 0;
    else
      amp2 = minmax2.amplitude(ind2);
    end;

    minmax1.deflection(i) =...
      abs(minmax1.amplitude(i) - mean([amp1 amp2]));
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
    responses = results.(type);
    if strcmp(type,'auc') && parms.auc_nbins>1
      for n=1:parms.auc_nbins
        tmp_responses = responses(:,:,n);
        toplabl = sprintf('%s   time range %d (%0.0f to %0.0f msec)',...
          toplabls{t},n,parms.auc_bins(n,1),parms.auc_bins(n,2));
        if n==parms.auc_nbins
          outfix = []; % entire range
        else
          outfix = sprintf('bin%d',n);
        end;
        plot_responses(parms,tmp_responses,type,toplabls{t},ylabls{t},outfix);
      end;
    else
      plot_responses(parms,responses,type,toplabls{t},ylabls{t});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_responses(parms,responses,type,toplabl,ylabl,outfix)
  if ~exist('outfix','var'), outfix = []; end;
  if strcmp(type,'latency'), parms.normflag = 0; end;
  figure; clf; hold on;
  for a=1:parms.nareas
    resps = squeeze(responses(a,:));
    if parms.normflag==1
      resps = resps/resps(parms.ind_max_cond);
    elseif parms.normflag==2
      [maxval,ind] = max(abs(resps));
      maxval = resps(ind);
      resps = resps/maxval;
    end;
    % plot amplitude vs. conditions
    if parms.linewidth>0
      plot(parms.condition_values,resps,'o-',...
        'Color',parms.area_colors{a},'LineWidth',parms.linewidth);
    else
      plot(parms.condition_values,resps,'o','Color',parms.area_colors{a});
    end;
  end;
  if parms.label_flag
    plot_response_labels(parms,type,toplabl,ylabl)
  end;
  outstem = [parms.outstem type];
  if ~isempty(outfix), outstem = [outstem '_' outfix]; end;
  if parms.normflag==1
    outstem = [outstem '_norm'];
  elseif parms.normflag==2
    outstem = [outstem '_normmax'];
  end;
  save_plot(parms,outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function plot_waveforms(wforms,parms,results)
  figure; clf; hold on;
  % plot waveform
  for c=1:results.nconditions
    for a=1:results.nareas
      plot(results.time,wforms(:,a,c),...
        'Color',parms.area_colors{a},'LineWidth',parms.linewidths(c));
    end;
  end;
  if parms.mark_peaks_flag
    % mark peaks
    for c=1:results.nconditions
      for a=1:results.nareas
        % mark minima
        npeaks = length(results.minima(a,c).latency);
        for i=1:npeaks
          amplitude = results.minima(a,c).amplitude(i);
          latency = results.minima(a,c).latency(i);
          plot(latency,amplitude,'O',...
            'Color',parms.area_colors{a},'LineWidth',parms.linewidth);
        end;
        % mark maxima
        npeaks = length(results.maxima(a,c).latency);
        for i=1:npeaks
          amplitude = results.maxima(a,c).amplitude(i);
          latency = results.maxima(a,c).latency(i);
          plot(latency,amplitude,'X',...
            'Color',parms.area_colors{a},'LineWidth',parms.linewidth);
        end;
      end;
    end;
  end;
  if parms.label_flag
    ylabl = sprintf('Source Amplitude (%s)',parms.units_wform);
    toplabl = 'Source Estimate Waveforms';
    plot_waveform_labels(parms,toplabl,ylabl);
  end;
  outstem = [parms.outstem 'wforms'];
  save_plot(parms,outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_cond_waveforms(wforms,parms,results)
  % plot waveform
  for c=1:results.nconditions
    figure(c); clf; hold on;
    for a=1:results.nareas
      plot(results.time,wforms(:,a,c),...
        'Color',parms.area_colors{a},'LineWidth',parms.linewidths(c));
    end;
    if ~parms.visible_flag, set(gcf,'Visible','off'); end;
  end;
  for c=1:results.nconditions
    figure(c);
    if parms.mark_peaks_flag
      % mark peaks
      for a=1:results.nareas
        % mark minima
        npeaks = length(results.minima(a,c).latency);
        for i=1:npeaks
          amplitude = results.minima(a,c).amplitude(i);
          latency = results.minima(a,c).latency(i);
          plot(latency,amplitude,'O',...
            'Color',parms.area_colors{a},'LineWidth',parms.linewidth);
        end;
        % mark maxima
        npeaks = length(results.maxima(a,c).latency);
        for i=1:npeaks
          amplitude = results.maxima(a,c).amplitude(i);
          latency = results.maxima(a,c).latency(i);
          plot(latency,amplitude,'X',...
            'Color',parms.area_colors{a},'LineWidth',parms.linewidth);
        end;
      end;
    end;
    if parms.label_flag
      ylabl = sprintf('Source Amplitude (%s)',parms.units_wform);
      toplabl = sprintf('%s = %0.2f',...
        parms.condition_label,parms.condition_values(c));
      plot_waveform_labels(parms,toplabl,ylabl);
    end;
    outstem = sprintf('%scond%0.2f_wforms',...
      parms.outstem,parms.condition_values(c));
    save_plot(parms,outstem);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_waveform_labels(parms,toplabl,ylabl)
  set(gca,'XLim',parms.xlim_wform);
  if ~isempty(parms.ylim_wform), set(gca,'YLim',parms.ylim_wform); end;
  set(gca,parms.fontargs{:});
  xlabel('time (msec)',parms.fontargs{:});
  ylabel(ylabl,parms.fontargs{:});
  title(toplabl,parms.fontargs{:});
  legend(parms.area_names,parms.fontargs{:},'Location',parms.legend_loc);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_plot(parms,outstem)
  if ~parms.visible_flag, set(gcf,'Visible','off'); end;
  print(gcf,'-dtiff',[outstem '.tif']);
  if parms.eps_flag, mmil_printeps(gcf,[outstem '.eps']); end;
  if ~parms.visible_flag, close(gcf); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_analysis(parms,results)
  % write values for each area and condition
  fname_out = [parms.outstem 'analysis.csv'];
  fid = fopen(fname_out,'wt');
  if fid==-1
    fprintf('%s: WARNNG: failed to open file %s for writing',fname_out);
  else
    fprintf(fid,'"area","condition","latency","amplitude","deflection","area under curve"\n');
    for c=1:results.nconditions
      for a=1:results.nareas
        latency = results.latency(a,c);
        amplitude = results.amplitude(a,c);
        deflection = results.deflection(a,c);
        auc = results.auc(a,c);
        fprintf(fid,'"%s",%0.2f,%0.2f,%0.2f,%0.2f,%0.2f\n',...
          parms.area_names{a},parms.condition_values(c),...
          latency,amplitude,deflection,auc);
      end;
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_peaks(parms,results)
  % write all peaks
  fname_out = [parms.outstem 'peaks.csv'];
  fid = fopen(fname_out,'wt');
  if fid==-1
    fprintf('%s: WARNNG: failed to open file %s for writing',fname_out);
  else
    fprintf(fid,'"area","condition","type","latency","amplitude","deflection"\n');
    for c=1:results.nconditions
      for a=1:results.nareas
        extypes = {'minimum','maximum'};
        for j=1:length(extypes)
          extype = extypes{j};
          if strcmp(extype,'minimum')
            minmax = results.minima(a,c);
          else
            minmax = results.maxima(a,c);
          end;          
          % list minima/maxima
          npeaks = length(minmax.latency);
          for i=1:npeaks
            latency = minmax.latency(i);
            amplitude = minmax.amplitude(i);
            deflection = minmax.deflection(i);
            fprintf(fid,'"%s",%0.2f,"%s",%0.2f,%0.2f,%0.2f\n',...
              parms.area_names{a},parms.condition_values(c),...
              extype,latency,amplitude,deflection);
          end;
        end;
      end;
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

