function resamp_results = MEG_MMIL_Sensor_GroupAvg_Resamp(fnames,varargin)
%function resamp_results = MEG_MMIL_Sensor_GroupAvg_Resamp(fnames,[options])
%
% Purpose: Calculate mean and variation across subjects
%   of a set of sensor waveforms using bootstrap resampling
%
% Required Input:
%   fnames: cell array of Sensor_GroupAvg results mat files
%     to be averaged together
%
% Optional Input ('key',arg pairs):
%  'outdir': output directory
%    {default = [pwd '/Sensor_GroupAvg_Resamp']}
%  'outstem': output file stem
%    {default = 'Sensor'}
%  'forceflag': [0|1] overwrite existing output
%    {default: 0}
%
% Output:
%   resamp_result: struct array with resampling results
%     including means, sem, confidence intervals, etc.
%
% Created:  09/18/13 by Don Hagler
% Last Mod: 12/09/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: separate, general function for plotting resamp_results
%% todo: plot histograms of onset latency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resamp_results = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(fnames,varargin);

mmil_mkdir(parms.outdir);
fname_results = [parms.outdir '/' parms.outstem '_GroupAvg_results.mat'];
if ~exist(fname_results,'file') || parms.forceflag
  % load previously calculated GroupAvg results for each condition
  %   (which contain waveforms for each subject)
  fprintf('%s: loading Sensor GroupAvg results for %s...\n',...
    mfilename,parms.outstem);
  [results,parms] = load_results(parms);
  % set parameters for resampling and plotting
  parms = set_parms(parms);
  % resample waveforms
  fprintf('%s: resampling waveforms for %s...\n',mfilename,parms.outstem);
  resamp_results = resamp_wforms(parms,results);
  % calculate stats for across-condition average
  if parms.cond_avg_flag
    fprintf('%s: averaging across conditions for %s...\n',...
      mfilename,parms.outstem);
    resamp_results = calc_cond_avg(parms,resamp_results);
  end;
  % save results in csv files
  if parms.csv_flag
    fprintf('%s: summarize results in csv files for %s...\n',...
      mfilename,parms.outstem);
    summarize_results(parms,resamp_results);
  end;
  % save results to mat file
  save(fname_results,'resamp_results','-v7.3');
else
  if parms.plot_wforms_flag || parms.plot_resps_flag
    fprintf('%s: loading Sensor GroupAvg results for %s...\n',...
      mfilename,parms.outstem);
    [results,parms] = load_results(parms);
    % set parameters for resampling and plotting
    parms = set_parms(parms);
  end;
  load(fname_results);
end;

% plot waveforms
if parms.plot_wforms_flag
  fprintf('%s: plotting waveforms for %s...\n',mfilename,parms.outstem);
  parms = plot_wforms(parms,resamp_results);
end;
% plot responses
if parms.plot_resps_flag
  fprintf('%s: plotting responses for %s...\n',mfilename,parms.outstem);
  parms = plot_responses(parms,resamp_results);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fnames,options)
  parms_filter = {...
    'fnames',fnames,[],...
    'outdir',[pwd '/Sensor_GroupAvg_Resamp'],[],...
    'outstem','Sensor',[],...
    'forceflag',false,[false true],...  
    'conditions',[0.15,0.71,0.95],[],...
    'condnames',{'LC','MC','HC'},[],...
    'cond_varname','contrast',[],...
    'cond_label','Luminance Contrast',[],...
    'cond_lim',[0,1],[],...
    'cond_avg_flag',true,[false true],...
    'niters',2000,[],...
    'bias_corr_flag',true,[false true],...
    'alpha',0.05,[0.001,0.1],...
    'ylim_wform',[-1,10],[],...
    'ylim_latency',[0,200],[],...
    'ylim_peak_norm',[-1,2],[],...
    'ylim_peak_amp',[-0.5,10],[],...
    'ylim_auc_norm',[-1,2],[],...
    'ylim_auc',[-10,50],[],...
    'xlim_wform',[-100,300],[],...
    'plot_wforms_flag',true,[false,true],...
    'plot_resps_flag',true,[false,true],...
    'mark_peak_flag',true,[false,true],...
    'mark_onset_flag',true,[false,true],...
    'onset_marker_width',2,[1,1000],...
    'onset_kappa',4.5,[1,100],...
    'onset_baseline_range',[-100,0],[],...
    'cond_colors',{'k','b','r','m','g','y','c'},[],...
    'wforms_err_flag',true,[false true],...
    'fill_err_flag',true,[false true],...
    'normflag',false,[false true],...
    'peak_pol',1,[-1,1],...
    'peak_range',[0,200],[],...
    'peak_mindiff',0.1,[],...
    'resfact',1,[],...
    'units','a.u.',[],...
    'norm_units','a.u.',[],...
    'resps_fig_size',[],[],...
    'resps_legend_loc','NorthEastOutside',[],...
    'resps_linewidth',1,[],...
    'resps_zero_line_flag',true,[false,true],...
    'wforms_fig_size',[],[],...
    'wforms_legend_loc','NorthEastOutside',[],...
    'wforms_linewidth',1,[],...
    'wforms_zero_line_flag',false,[false,true],...
    'eps_flag',true,[false true],...
    'tif_dpi',300,[10,10000],...
    'fontname','Arial',[],...
    'fontsize',12,[],...
    'ranges_flag',false,[false true],...
    'onset_flag',true,[false true],...
    'peak_flag',true,[false true],...
    'auc_flag',true,[false true],...
    'ylabel_flag',true,[false true],...
    'title_flag',true,[false true],...
    'legend_flag',true,[false true],...
    'markersize',6,[],...
    'csv_flag',true,[false true],...
    'visible_flag',false,[false true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  if parms.normflag
    parms.units = parms.norm_units;
  end;
  parms.nfiles = length(parms.fnames);
  parms.f = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results,parms] = load_results(parms);
  % load results
  clear results;
  for f=1:parms.nfiles
    tmp = load(parms.fnames{f});
    results{f} = tmp.results;
    if f==1
      % get information from first results struct
      % NOTE: will have errors if other results have different N, etc.
      parms.time = tmp.results.time*1000;
      parms.nsubs = tmp.results.N;
      parms.ntpoints = length(parms.time);
      parms.nconds = tmp.results.nconditions;
    end;      
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_parms(parms)
  if ~parms.onset_flag
    parms.mark_onset_flag = false;
  end;
  if ~parms.peak_flag
    parms.mark_peak_flag = false;
  end;

  % parameters for resampling
  parms.resamp_parms = [];
  parms.resamp_parms.peak_flag = parms.peak_flag;
  parms.resamp_parms.onset_flag = parms.onset_flag;
  parms.resamp_parms.auc_flag = parms.auc_flag;
  parms.resamp_parms.niters = parms.niters;
  parms.resamp_parms.randseed_flag = 0;
  parms.resamp_parms.bias_corr_flag = parms.bias_corr_flag;
  parms.resamp_parms.alpha = parms.alpha;
  % for peaks
  parms.resamp_parms.peak_pol = parms.peak_pol;
  parms.resamp_parms.peak_range = parms.peak_range;
  parms.resamp_parms.peak_mindiff = parms.peak_mindiff;
  parms.resamp_parms.resfact = parms.resfact;
  % for onset
  parms.resamp_parms.onset_baseline_range = parms.onset_baseline_range;
  parms.resamp_parms.onset_baseline_collapse_flag = 2;
  parms.resamp_parms.onset_minimum = 40;
  parms.resamp_parms.onset_polarity = 1;
  parms.resamp_parms.onset_baseline_flag = 2;
  parms.resamp_parms.onset_method = 'quartiles';
  parms.resamp_parms.onset_kappa = parms.onset_kappa;
  parms.resamp_parms.time = parms.time;

  % parameters for plotting waveforms
  parms.plot_parms = [];
  if parms.ylabel_flag
    parms.plot_parms.ylabel = sprintf('Response Amplitude (%s)',...
        parms.units);
  end;
  parms.plot_parms.xlim = parms.xlim_wform;
  parms.plot_parms.legend_loc = parms.wforms_legend_loc;
  parms.plot_parms.outdir = parms.outdir;
  parms.plot_parms.baseline_flag = 1;
  parms.plot_parms.linewidth = parms.wforms_linewidth;
  parms.plot_parms.fill_err_flag = parms.fill_err_flag;
  parms.plot_parms.fill_alpha = 0.4;
  parms.plot_parms.errbar_interval = 5;
  parms.plot_parms.visible_flag = parms.visible_flag;
  parms.plot_parms.relative_err_flag = 0;
  parms.plot_parms.zero_line_flag = parms.wforms_zero_line_flag;
  parms.plot_parms.eps_flag = parms.eps_flag;
  parms.plot_parms.fig_size = parms.wforms_fig_size;
  parms.plot_parms.fontname = parms.fontname;
  parms.plot_parms.fontsize = parms.fontsize;
  parms.plot_parms.tif_dpi = parms.tif_dpi;
  parms.plot_parms.time = parms.time;

  % for plotting responses
  parms.rplot_parms = [];
  parms.rplot_parms.outdir = parms.outdir;
  parms.rplot_parms.colors = [];
  parms.rplot_parms.conditions = parms.conditions;
  parms.rplot_parms.xlabel = parms.cond_label;
  parms.rplot_parms.xlim = parms.cond_lim;
  parms.rplot_parms.legend_loc = parms.resps_legend_loc;
  parms.rplot_parms.linewidth = parms.resps_linewidth;
  parms.rplot_parms.visible_flag = parms.visible_flag;
  parms.rplot_parms.relative_err_flag = 0;
  parms.rplot_parms.zero_line_flag = parms.resps_zero_line_flag;
  parms.rplot_parms.eps_flag = parms.eps_flag;
  parms.rplot_parms.fig_size = parms.resps_fig_size;
  parms.rplot_parms.fontname = parms.fontname;
  parms.rplot_parms.fontsize = parms.fontsize;
  parms.rplot_parms.tif_dpi = parms.tif_dpi;
  parms.rplot_parms.markersize = parms.markersize;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resamp_results = resamp_wforms(parms,results)
  clear resamp_results;
  resamp_args = mmil_parms2args(parms.resamp_parms);
  % compile waveforms
  wforms = zeros(parms.ntpoints,parms.nconds,parms.nsubs);
  for f=1:parms.nfiles
    tmp_results = results{f};
    % compile data for all subjects
    for s=1:parms.nsubs
      tmp_wforms = tmp_results.subjdata(s).wforms;
      wforms(:,:,s) = wforms(:,:,s) + tmp_wforms;
    end;
  end;
  wforms = wforms / parms.nfiles;
  % normalize waveforms by max across time and conditions
  if parms.normflag
    max_vals = squeeze(max(squeeze(max(abs(wforms),[],1)),[],1));
    max_mat = repmat(reshape(max_vals,[1,1,parms.nsubs]),...
                     [parms.ntpoints,parms.nconds,1]);
    wforms = wforms ./ max_mat;
  end;
  % resample waveforms
  for c=1:parms.nconds
    tmp_wforms = squeeze(wforms(:,c,:));
    tmp_results = ts_bootstrap_wforms(tmp_wforms,resamp_args{:});
    tmp_results.cond = parms.condnames{c};
    resamp_results(c) = tmp_results;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resamp_results = calc_cond_avg(parms,resamp_results)
  cond_name = 'average';
  if parms.peak_flag
    % calculate average peak amplitude and latency
    peak = [];
    % calculate average amplitudes and latencies
    for c=1:parms.nconds
      if isempty(peak)
        peak = resamp_results(c).peak;
      else
        peak.amplitudes = peak.amplitudes + ...
          resamp_results(c).peak.amplitudes;
        peak.latencies = peak.latencies + ...
          resamp_results(c).peak.latencies;
        if parms.bias_corr_flag
          peak.amplitudes_jackknife = peak.amplitudes_jackknife + ...
            resamp_results(c).peak.amplitudes_jackknife;
          peak.latencies_jackknife = peak.latencies_jackknife + ...
            resamp_results(c).peak.latencies_jackknife;
        end;
        peak.amplitude_mean = peak.amplitude_mean + ...
          resamp_results(c).peak.amplitude_mean;
        peak.latency_mean = peak.latency_mean + ...
          resamp_results(c).peak.latency_mean;
      end;
    end;
    peak.amplitudes = peak.amplitudes / parms.nconds;
    peak.latencies = peak.latencies / parms.nconds;
    if parms.bias_corr_flag
      peak.amplitudes_jackknife = peak.amplitudes_jackknife / parms.nconds;
      peak.latencies_jackknife = peak.latencies_jackknife / parms.nconds;
    end;
    peak.amplitude_mean = peak.amplitude_mean / parms.nconds;
    peak.latency_mean = peak.latency_mean / parms.nconds;
    peak.amplitude = ts_nan_mean(peak.amplitudes);
    peak.latency = ts_nan_mean(peak.latencies);
    % calculate bootstrap sem
    peak.amplitude_sem = ts_nan_std(peak.amplitudes,0,2);
    peak.latency_sem = ts_nan_std(peak.latencies,0,2);
    % calculate confidence intervals
    if parms.bias_corr_flag
      peak.amplitude_ci = ts_calc_bs_ci_bias(peak.amplitudes,...
        peak.amplitudes_jackknife,peak.amplitude_mean,parms.alpha);
      peak.latency_ci = ts_calc_bs_ci_bias(peak.latencies,...
        peak.latencies_jackknife,peak.latency_mean,parms.alpha);
    else
      peak.amplitude_ci = ts_calc_bs_ci(peak.amplitudes,parms.alpha);
      peak.latency_ci = ts_calc_bs_ci(peak.latencies,parms.alpha);
    end;
    % add to resamp_results
    resamp_results(parms.nconds+1).peak = peak;
    resamp_results(parms.nconds+1).cond = cond_name;
  end;

  if parms.onset_flag
    % calculate average onset latency
    onset = [];
    % calculate average latencies
    for c=1:parms.nconds
      if isempty(onset)
        onset = resamp_results(c).onset;
      else
        onset.latencies = onset.latencies + ...
          resamp_results(c).onset.latencies;
        if parms.bias_corr_flag
          onset.latencies_jackknife = onset.latencies_jackknife + ...
            resamp_results(c).onset.latencies_jackknife;
        end;
        onset.latency_mean = onset.latency_mean + ...
          resamp_results(c).onset.latency_mean;
      end;
    end;
    onset.latencies = onset.latencies / parms.nconds;
    if parms.bias_corr_flag
      onset.latencies_jackknife = onset.latencies_jackknife / parms.nconds;
    end;
    onset.latency_mean = onset.latency_mean / parms.nconds;
    onset.latency = ts_nan_mean(onset.latencies);
    % calculate bootstrap sem
    onset.latency_sem = ts_nan_std(onset.latencies,0,2);
    % calculate confidence intervals
    if parms.bias_corr_flag
      onset.latency_ci = ts_calc_bs_ci_bias(onset.latencies,...
        onset.latencies_jackknife,onset.latency_mean,parms.alpha);
    else
      onset.latency_ci = ts_calc_bs_ci(onset.latencies,parms.alpha);
    end;
    % add to resamp_results
    resamp_results(parms.nconds+1).onset = onset;
    resamp_results(parms.nconds+1).cond = cond_name;
  end;

  if parms.auc_flag
    % calculate average auc value
    auc = [];
    % calculate average vals
    for c=1:parms.nconds
      if isempty(auc)
        auc = resamp_results(c).auc;
      else
        auc.vals = auc.vals + ...
          resamp_results(c).auc.vals;
        if parms.bias_corr_flag
          auc.vals_jackknife = auc.vals_jackknife + ...
            resamp_results(c).auc.vals_jackknife;
        end;
        auc.val_mean = auc.val_mean + ...
          resamp_results(c).auc.val_mean;
      end;
    end;
    auc.vals = auc.vals / parms.nconds;
    if parms.bias_corr_flag
      auc.vals_jackknife = auc.vals_jackknife / parms.nconds;
    end;
    auc.val_mean = auc.val_mean / parms.nconds;
    auc.mean = ts_nan_mean(auc.vals);
    % calculate bootstrap sem
    auc.sem = ts_nan_std(auc.vals,0,2);
    % calculate confidence intervals
    if parms.bias_corr_flag
      auc.ci = ts_calc_bs_ci_bias(auc.vals,...
        auc.vals_jackknife,auc.val_mean,parms.alpha);
    else
      auc.ci = ts_calc_bs_ci(auc.vals,parms.alpha);
    end;
    % add to resamp_results
    resamp_results(parms.nconds+1).auc = auc;
    resamp_results(parms.nconds+1).cond = cond_name;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_wforms(parms,resamp_results)
  % plot waveforms with all conds for each area
  fprintf('%s: plotting waveforms for %s...\n',mfilename,parms.outstem);
  % set plotting parameters
  plot_parms = parms.plot_parms;
  plot_parms.outstem = sprintf('%s_wforms',parms.outstem);
  if parms.title_flag
    plot_parms.title = ...
      sprintf('%s group response',regexprep(parms.outstem,'_',' '));
  end;
  if parms.legend_flag
    plot_parms.condnames = parms.condnames;
  end;
  plot_parms.ylim = parms.ylim_wform;
  plot_parms.colors = parms.cond_colors;
  % compile data
  wforms = zeros(parms.ntpoints,parms.nconds);
  if parms.wforms_err_flag
    plot_parms.wforms_err = zeros(parms.ntpoints,parms.nconds,2);
  else
    plot_parms.wforms_err = [];
  end;
  plot_parms.points = [];
  plot_parms.ranges = [];
  for c=1:parms.nconds
    tmp_results = resamp_results(c);
    wforms(:,c) = tmp_results.wforms_mean;
    if parms.wforms_err_flag
      plot_parms.wforms_err(:,c,:) = reshape(tmp_results.wforms_bs_ci,...
        [parms.ntpoints,1,2]);
    end;
    if parms.ranges_flag & ~isempty(tmp_results.ranges_sig)
      fprintf('%s: significant range for r=%d, c=%d\n',...
        mfilename,r,c);
      plot_parms.ranges(c).onset = ...
        plot_parms.time(tmp_results.ranges_sig.t_onset);
      plot_parms.ranges(c).offset = ...
        plot_parms.time(tmp_results.ranges_sig.t_offset);
    else
      plot_sig_ranges(c).onset = [];
      plot_sig_ranges(c).offset = [];
    end;
  end;
  % plot waveforms
  figure(parms.f); parms.f = parms.f + 1;
  plot_args = mmil_parms2args(plot_parms);
  ts_plot_wforms(wforms,plot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_responses(parms,resamp_results)
  % plot peak latency vs. condition value
  if parms.peak_flag
    parms = plot_latency(parms,resamp_results,'peak');
  end;
  % plot onset latency vs. condition value
  if parms.resamp_parms.onset_flag
    parms = plot_latency(parms,resamp_results,'onset');
  end;
  % plot peak amplitude vs. condition value
  if parms.peak_flag
    for normflag=0:1
      parms = plot_peak_amplitude(parms,resamp_results,normflag);
    end;
  end;
  % plot area under the curve vs. condition value
  if parms.resamp_parms.auc_flag
    for normflag=0:1
      parms = plot_area_under_curve(parms,resamp_results,normflag);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_latency(parms,resamp_results,type_str)
  fprintf('%s: plotting %s latency for %s...\n',...
    mfilename,type_str,parms.outstem);
  rplot_parms = parms.rplot_parms;
  rplot_parms.outstem = sprintf('%s_%s_latency',parms.outstem,type_str);
  if parms.title_flag
    rplot_parms.title = ...
      sprintf('%s latency vs. %s',type_str,parms.cond_varname);
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = 'Latency (msec)';
  end;
  if parms.legend_flag
    rplot_parms.roinames = [];
  end;
  rplot_parms.ylim = parms.ylim_latency;
  rplot_parms.normflag = 0;
  rplot_parms.responses_err = zeros(1,parms.nconds,2);
  responses = zeros(1,parms.nconds);
  for c=1:parms.nconds
    responses(c) = resamp_results(c).(type_str).latency;
    rplot_parms.responses_err(1,c,:) = resamp_results(c).(type_str).latency_ci;
  end;
  figure(parms.f); parms.f = parms.f + 1;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_peak_amplitude(parms,resamp_results,normflag)
  rplot_parms = parms.rplot_parms;
  if normflag
    fprintf('%s: plotting normalized peak amplitude for %s...\n',...
      mfilename,parms.outstem);
    rplot_parms.outstem = sprintf('%s_peak_amplitude_norm',parms.outstem);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('normalized peak amplitude vs. %s',parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = ...
        sprintf('Normalized Peak Amplitude (%s)',parms.norm_units);
    end;
    rplot_parms.ylim = parms.ylim_peak_norm;
    rplot_parms.normflag = 2;
  else
    fprintf('%s: plotting peak amplitude for %s...\n',...
      mfilename,parms.outstem);
    rplot_parms.outstem = sprintf('%s_peak_amplitude',parms.outstem);
    if parms.title_flag
      rplot_parms.title = sprintf('peak amplitude vs. %s',parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Peak Amplitude (%s)',...
        parms.units);
    end;
    rplot_parms.ylim = parms.ylim_peak_amp;
    rplot_parms.normflag = 0;
  end;
  if parms.legend_flag
    rplot_parms.roinames = [];
  end;
  rplot_parms.responses_err = zeros(1,parms.nconds,2);
  responses = zeros(1,parms.nconds);
  for c=1:parms.nconds
    responses(c) = resamp_results(c).peak.amplitude;
    rplot_parms.responses_err(1,c,:) = resamp_results(c).peak.amplitude_ci;
  end;
  figure(parms.f); parms.f = parms.f + 1;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_area_under_curve(parms,resamp_results,normflag)
  rplot_parms = parms.rplot_parms;
  if normflag
    fprintf('%s: plotting normalized area under curve for %s...\n',...
      mfilename,parms.outstem);
    rplot_parms.outstem = sprintf('%s_auc_norm',parms.outstem);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('normalized area under curve vs. %s',parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = ...
        sprintf('Normalized Area Under Curve (%s)',...
          parms.norm_units);
    end;
    rplot_parms.ylim = parms.ylim_auc_norm;
    rplot_parms.normflag = 2;
  else
    fprintf('%s: plotting area under curve for %s...\n',...
      mfilename,parms.outstem);
    rplot_parms.outstem = sprintf('%s_auc',parms.outstem);
    if parms.title_flag
      rplot_parms.title = sprintf('area under curve vs. %s',parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Area Under Curve (%s)',...
        parms.units);
    end;
    rplot_parms.ylim = parms.ylim_auc;
    rplot_parms.normflag = 0;
  end;
  if parms.legend_flag
    rplot_parms.roinames = [];
  end;
  rplot_parms.responses_err = zeros(1,parms.nconds,2);
  responses = zeros(1,parms.nconds);
  for c=1:parms.nconds
    responses(c) = resamp_results(c).auc.mean;
    rplot_parms.responses_err(1,c,:) = resamp_results(1,c).auc.ci;
  end;
  figure(parms.f); parms.f = parms.f + 1;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_results(parms,resamp_results)
  outstem = [parms.outdir '/' parms.outstem];
  ts_summarize_resamp_results(resamp_results,'fstem_out',outstem);
return;

