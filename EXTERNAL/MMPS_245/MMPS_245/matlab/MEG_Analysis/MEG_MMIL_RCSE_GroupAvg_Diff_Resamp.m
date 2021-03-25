function resamp_results = MEG_MMIL_RCSE_GroupAvg_Diff_Resamp(rootdir,varargin)
%function resamp_results = MEG_MMIL_RCSE_GroupAvg_Diff_Resamp(rootdir,[options])
%
% Purpose: Compare two sets of RCSE waveforms across subjects
%   using bootstrap resampling
%
% Required Input:
%   rootdir: root directory containing RCSE_GroupAvg files
%
% Optional Input ('key',arg pairs):
%   'outdir': output directory
%     {default = [pwd '/RCSE_GroupAvg_Diff_Resamp']}
%   'outstem': output file stem
%     {default = 'RCSE_diff'}
%   'prefix': RCSE prefix
%     {default = 'RCSE_prior'}
%   'infix_list': cell array of RCSE infixes
%     must have 2 entries
%     {default = {'uplow1','uplow2'}}
%   'outfix_list': cell array of outfixes used for legends
%     must have 2 entries
%     {default = {'upper','lower'}}
%
% Output:
%   resamp_result: struct array with resampling results
%     including means, sem, confidence intervals, etc.
%
% Created:  02/28/12 by Don Hagler
% Last Mod: 12/09/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: call ts_permute_wforms if ranges_flag = 1
% Purpose: Compare two sets of RCSE waveforms across subjects
%   using bootstrap resampling and permutation testing

%% todo: separate, general function for plotting resamp_results
%%       compatible with both this function and MEG_MMIL_dSROI_GroupAvg_Diff_Resamp
%% todo: add resamp_parms, plot_parms,and rplot_parms
%%       to parms_filter, use tags to distribute
%% todo: plot histograms of onset latency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
resamp_results = [];
parms = check_input(rootdir,varargin);

mmil_mkdir(parms.outdir);
fname_results = [parms.outdir '/' parms.outstem '_GroupAvg_results.mat'];
if ~exist(fname_results,'file') || parms.forceflag
  % load previously calculated GroupAvg results
  %   (which contain waveforms for each subject)
  fprintf('%s: loading RCSE GroupAvg results for %s...\n',mfilename,parms.outstem);
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
  % summarize results
  if parms.csv_flag
    fprintf('%s: summarize results in csv files for %s...\n',...
      mfilename,parms.outstem);
    summarize_results(parms,resamp_results);
  end;
  % save results to mat file
  save(fname_results,'resamp_results','-v7.3');
else
  if parms.plot_wforms_flag || parms.plot_resps_flag
    fprintf('%s: loading RCSE GroupAvg results for %s...\n',...
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

function parms = check_input(rootdir,options)
  parms_filter = {...
    'rootdir',rootdir,[],...
    'outdir',[pwd '/RCSE_GroupAvg_Diff_Resamp'],[],...
    'outstem','RCSE_diff',[],...
    'prefix','RCSE_prior',[],...
    'infix_list',{'uplow1','uplow2'},[],...
    'outfix_list',{'upper','lower'},[],...
    'area_names',{'V1','V2','V3'},[],...
    'condnames',{'LC','MC','HC'},[],...
    'conditions',[0.15,0.71,0.95],[],...
    'cond_varname','contrast',[],...
    'cond_label','Luminance Contrast',[],...
    'cond_lim',[0,1],[],...
    'cond_avg_flag',true,[false true],...
    'niters',2000,[],...
    'bias_corr_flag',true,[false true],...
    'alpha',0.05,[0.001,0.1],...
    'ylim_wform',[-20,10],[],...
    'ylim_latency',[0,200],[],...
    'ylim_peak_norm',[-1,2],[],...
    'ylim_peak_amp',[-20,10],[],...
    'ylim_auc_norm',[-1,2],[],...
    'ylim_auc',[-10,50],[],...
    'ylim_wform_diff',[-5,5],[],...
    'ylim_latency_diff',[-100,100],[],...
    'ylim_peak_norm_diff',[-5,5],[],...
    'ylim_peak_amp_diff',[-5,5],[],...
    'ylim_auc_norm_diff',[-5,5],[],...
    'ylim_auc_diff',[-30,30],[],...
    'xlim_wform',[-100,300],[],...
    'ranges_color','y',[],...
    'plot_wforms_flag',true,[false,true],...
    'plot_resps_flag',true,[false,true],...
    'mark_peak_flag',true,[false,true],...
    'mark_onset_flag',true,[false,true],...
    'onset_marker_width',2,[1,1000],...
    'onset_kappa',3,[1,100],...
    'onset_baseline_range',[-100,0],[],...
    'area_colors',{'b','g','r'},[],...
    'pair_colors',{'r','b'},[],...
    'cond_colors',{'k','b','r'},[],...
    'fill_err_flag',true,[false true],...
    'fill_alpha',0.4,[0,1],...
    'normflag',false,[false true],...
    'peak_pol',-1,[-1,1],...
    'peak_range',[0,200],[],...
    'peak_mindiff',0.25,[],...
    'resfact',1,[],...
    'units','nA M',[],...
    'norm_units','a.u.',[],...
    'wforms_err_flag',true,[false true],...
    'resps_fig_size',[],[],...
    'resps_legend_loc','NorthEast',[],...
    'resps_linewidth',1,[],...
    'resps_zero_line_flag',true,[false,true],...
    'wforms_fig_size',[],[],...
    'wforms_legend_loc','SouthEast',[],...
    'wforms_linewidth',1,[],...
    'wforms_zero_line_flag',false,[false,true],...
    'eps_flag',true,[false true],...
    'tif_dpi',300,[10,10000],...
    'fontname','Arial',[],...
    'fontsize',12,[],...
    'onset_flag',true,[false true],...
    'peak_flag',true,[false true],...
    'auc_flag',true,[false true],...
    'ylabel_flag',true,[false true],...
    'title_flag',true,[false true],...
    'legend_flag',true,[false true],...
    'markersize',6,[],...
    'axes_linewidth',1,[],...
    'csv_flag',true,[false true],...
    'visible_flag',false,[false true],...
    'forceflag',false,[false true],...
  ...
    'ranges_flag',false,[false true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  if parms.normflag
    parms.units = parms.norm_units;
  end;
  if ~exist(parms.rootdir,'dir')
    error('rootdir %s not found',parms.rootdir);
  end;
  if length(parms.infix_list)~=2
    error('infix_list must have two conditions for comparison');
  end;
  if length(parms.outfix_list)~=2
    error('length of outfix_list must match infix_list');
  end;
  parms.pval_range = linspace(1/parms.niters,1,parms.niters);
  parms.f = 1;
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
  %resamp_parms.ranges_flag = parms.ranges_flag;
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
  parms.resamp_parms.onset_polarity = -1;
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
  parms.plot_parms.ylim = parms.ylim_wform;
  parms.plot_parms.xlim = parms.xlim_wform;
  parms.plot_parms.legend_loc = parms.wforms_legend_loc;
  parms.plot_parms.outdir = parms.outdir;
  parms.plot_parms.baseline_flag = 0;
  parms.plot_parms.linewidth = parms.wforms_linewidth;
  parms.plot_parms.axes_linewidth = parms.axes_linewidth;
  parms.plot_parms.fill_err_flag = parms.fill_err_flag;
  parms.plot_parms.fill_alpha = parms.fill_alpha;
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
  parms.rplot_parms.colors = parms.area_colors;
  parms.rplot_parms.conditions = parms.conditions;
  parms.rplot_parms.xlabel = parms.cond_label;
  parms.rplot_parms.xlim = parms.cond_lim;
  parms.rplot_parms.legend_loc = parms.resps_legend_loc;
  parms.rplot_parms.linewidth = parms.resps_linewidth;
  parms.rplot_parms.axes_linewidth = parms.axes_linewidth;
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

function [results,parms] = load_results(parms);
  clear results;
  for p=1:2
    infix = parms.infix_list{p};
    if isempty(infix)
      fname = sprintf('%s/%s_GroupAvg_results.mat',...
        parms.rootdir,parms.prefix);
    else
      fname = sprintf('%s/%s_%s_GroupAvg_results.mat',...
        parms.rootdir,parms.prefix,parms.infix_list{p});
    end;
    tmp = load(fname);
    results(p) = tmp.results;
  end;
  parms.nsubs = results(1).N;
  parms.ntpoints = length(results(1).time);
  parms.nareas = results(1).nareas;
  parms.nconds = results(1).nconditions;
  parms.time = results(1).time*1000;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resamp_results = resamp_wforms(parms,results)
  clear resamp_results;
  resamp_args = mmil_parms2args(parms.resamp_parms);
  for a=1:parms.nareas
    % compile waveforms
    wforms = zeros(parms.ntpoints,parms.nconds,parms.nsubs,2);
    for s=1:parms.nsubs
      for p=1:2
        wforms(:,:,s,p) = squeeze(results(p).subjdata(s).wform(:,a,:));
      end;
    end;
    % normalize waveforms by max across time and conditions
    %  averaged across the pair of prefix+infix (to be unbiased)
    if parms.normflag
      max_vals = mean(squeeze(max(squeeze(max(abs(wforms),[],1)),[],1)),2);
      max_mat = repmat(reshape(max_vals,[1,1,parms.nsubs]),...
        [parms.ntpoints,parms.nconds,1,2]);
      wforms = wforms ./ max_mat;
    end;
    % resample waveforms
    for c=1:parms.nconds
      tmp_wforms = squeeze(wforms(:,c,:,:));
      tmp_results = ts_bootstrap_wforms(tmp_wforms,resamp_args{:});
      tmp_results.roi = parms.area_names{a};
      tmp_results.cond = parms.condnames{c};
      resamp_results(a,c) = tmp_results;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resamp_results = calc_cond_avg(parms,resamp_results)
  cond_name = 'average';
  for a=1:parms.nareas
    area_name = resamp_results(a,1).roi;
    if parms.peak_flag
      % calculate average peak amplitude and latency
      tags = {'peak1','peak2','peak'};
      for i=1:length(tags)
        tag = tags{i};
        peak = [];
        % calculate average amplitudes and latencies
        for c=1:parms.nconds
          if isempty(peak)
            peak = resamp_results(a,c).(tag);
          else
            peak.amplitudes = peak.amplitudes + ...
              resamp_results(a,c).(tag).amplitudes;
            peak.latencies = peak.latencies + ...
              resamp_results(a,c).(tag).latencies;
            if parms.bias_corr_flag
              peak.amplitudes_jackknife = peak.amplitudes_jackknife + ...
                resamp_results(a,c).(tag).amplitudes_jackknife;
              peak.latencies_jackknife = peak.latencies_jackknife + ...
                resamp_results(a,c).(tag).latencies_jackknife;
            end;
            peak.amplitude_mean = peak.amplitude_mean + ...
              resamp_results(a,c).(tag).amplitude_mean;
            peak.latency_mean = peak.latency_mean + ...
              resamp_results(a,c).(tag).latency_mean;
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
        % calculate null pvals if peak (not peak1 or peak2)
        if strcmp(tag,'peak')
          if parms.bias_corr_flag
            peak.amplitude_pval = ...
              ts_calc_bs_null_pval_bias(peak.amplitudes,...
                peak.amplitudes_jackknife,...
                peak.amplitude_mean,parms.pval_range);
            peak.latency_pval = ...
              ts_calc_bs_null_pval_bias(peak.latencies,...
                peak.latencies_jackknife,...
                peak.latency_mean,parms.pval_range);
          else
            peak.amplitude_pval = ...
              ts_calc_bs_null_pval(peak.amplitudes,parms.pval_range);
            peak.latency_pval = ...
              ts_calc_bs_null_pval(peak.latencies,parms.pval_range);
          end;
        end;
        % add to resamp_results
        resamp_results(a,parms.nconds+1).(tag) = peak;
        resamp_results(a,parms.nconds+1).roi = area_name;
        resamp_results(a,parms.nconds+1).cond = cond_name;
      end;
    end;

    if parms.onset_flag
      % calculate average onset latency
      tags = {'onset1','onset2','onset'};
      for i=1:length(tags)
        tag = tags{i};
        onset = [];
        % calculate average latencies
        for c=1:parms.nconds
          if isempty(onset)
            onset = resamp_results(a,c).(tag);
          else
            onset.latencies = onset.latencies + ...
              resamp_results(a,c).(tag).latencies;
            if parms.bias_corr_flag
              onset.latencies_jackknife = onset.latencies_jackknife + ...
                resamp_results(a,c).(tag).latencies_jackknife;
            end;
            onset.latency_mean = onset.latency_mean + ...
              resamp_results(a,c).(tag).latency_mean;
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
        % calculate null pvals if onset (not onset1 or onset2)
        if strcmp(tag,'onset')
          if parms.bias_corr_flag
            onset.latency_pval = ...
              ts_calc_bs_null_pval_bias(onset.latencies,...
                onset.latencies_jackknife,...
                onset.latency_mean,parms.pval_range);
          else
            onset.latency_pval = ...
              ts_calc_bs_null_pval(onset.latencies,parms.pval_range);
          end;
        end;
        % add to resamp_results
        resamp_results(a,parms.nconds+1).(tag) = onset;
        resamp_results(a,parms.nconds+1).roi = area_name;
        resamp_results(a,parms.nconds+1).cond = cond_name;
      end;
    end;

    if parms.auc_flag
      % calculate average auc value
      tags = {'auc1','auc2','auc'};
      for i=1:length(tags)
        tag = tags{i};
        auc = [];
        % calculate average vals
        for c=1:parms.nconds
          if isempty(auc)
            auc = resamp_results(a,c).(tag);
          else
            auc.vals = auc.vals + ...
              resamp_results(a,c).(tag).vals;
            if parms.bias_corr_flag
              auc.vals_jackknife = auc.vals_jackknife + ...
                resamp_results(a,c).(tag).vals_jackknife;
            end;
            auc.val_mean = auc.val_mean + ...
              resamp_results(a,c).(tag).val_mean;
          end;
        end;
        auc.vals = auc.vals / parms.nconds;
        if parms.bias_corr_flag
          auc.vals_jackknife = auc.vals_jackknife / parms.nconds;
        end;
        auc.val_mean = auc.val_mean / parms.nconds;
        auc.mean = ts_nan_mean(auc.vals);
        % calculate bootstrap sem
        auc.mean_sem = ts_nan_std(auc.vals,0,2);
        % calculate confidence intervals
        if parms.bias_corr_flag
          auc.ci = ts_calc_bs_ci_bias(auc.vals,...
            auc.vals_jackknife,auc.val_mean,parms.alpha);
        else
          auc.ci = ts_calc_bs_ci(auc.vals,parms.alpha);
        end;
        % calculate null pvals if auc (not auc1 or auc2)
        if strcmp(tag,'auc')
          if parms.bias_corr_flag
            auc.mean_pval = ...
              ts_calc_bs_null_pval_bias(auc.vals,...
                auc.vals_jackknife,...
                auc.val_mean,parms.pval_range);
          else
            auc.mean_pval = ...
              ts_calc_bs_null_pval(auc.vals,parms.pval_range);
          end;
        end;
        % add to resamp_results
        resamp_results(a,parms.nconds+1).(tag) = auc;
        resamp_results(a,parms.nconds+1).roi = area_name;
        resamp_results(a,parms.nconds+1).cond = cond_name;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_wforms(parms,resamp_results)
  % plot paired waveforms for each condition and area
  parms = plot_wforms_conds(parms,resamp_results);
  % plot difference waveforms for each area
  parms = plot_wforms_diffs(parms,resamp_results);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_wforms_conds(parms,resamp_results)
  % plot pairs of waveforms for each area and condition separately
  for a=1:parms.nareas
    area = parms.area_names{a};
    for c=1:parms.nconds
      plot_parms = parms.plot_parms;
      plot_parms.outstem = sprintf('%s_%s_cond%d_wforms',...
        parms.outstem,parms.area_names{a},c);
      if parms.title_flag
        plot_parms.title = sprintf('group %s response; %s = %0.2f',...
          parms.area_names{a},parms.cond_varname,parms.rplot_parms.conditions(c));
      end;
      wforms = ...
        cat(2,resamp_results(a,c).wforms1_mean,resamp_results(a,c).wforms2_mean);
      if parms.wforms_err_flag
        plot_parms.wforms_err = ...
          cat(2,resamp_results(a,c).wforms1_bs_ci,resamp_results(a,c).wforms2_bs_ci);
      end;
      if parms.legend_flag
        plot_parms.condnames = parms.outfix_list;
      end;
      if parms.ylabel_flag
        plot_parms.ylabel = sprintf('Response Amplitude (%s)',...
            parms.units);
      end;
      plot_parms.ylim = parms.ylim_wform;
      plot_ranges = [];
      if parms.mark_onset_flag
        plot_ranges(1).onset = resamp_results(a,c).onset1.latency;
        plot_ranges(1).offset = resamp_results(a,c).onset1.latency + parms.onset_marker_width;
        plot_ranges(2).onset = resamp_results(a,c).onset2.latency;
        plot_ranges(2).offset = resamp_results(a,c).onset2.latency + parms.onset_marker_width;
        plot_parms.ranges_color = [];
      elseif parms.ranges_flag & ~isempty(resamp_results(a,c).ranges_sig)
        plot_ranges(1).onset = plot_parms.time(resamp_results(a,c).ranges_sig.t_onset);
        plot_ranges(1).offset = plot_parms.time(resamp_results(a,c).ranges_sig.t_offset);
        plot_ranges(2).onset = [];
        plot_ranges(2).offset = [];
        plot_parms.ranges_color = parms.ranges_color;
      end;
      plot_points = [];
      if parms.mark_peak_flag
        plot_points(1).amplitude = resamp_results(a,c).peak1.amplitude;
        plot_points(1).latency = resamp_results(a,c).peak1.latency;
        plot_points(1).type = parms.resamp_parms.peak_pol;
        plot_points(2).amplitude = resamp_results(a,c).peak2.amplitude;
        plot_points(2).latency = resamp_results(a,c).peak2.latency;
        plot_points(2).type = parms.resamp_parms.peak_pol;
      end;
      plot_parms.ranges = plot_ranges;
      plot_parms.points = plot_points;
      plot_parms.colors = parms.pair_colors;
      plot_args = mmil_parms2args(plot_parms);
      figure(parms.f); parms.f = parms.f + 1;
      ts_plot_wforms(wforms,plot_args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_wforms_diffs(parms,resamp_results)
  for a=1:parms.nareas
    plot_parms = parms.plot_parms;
    wforms = [];
    if parms.wforms_err_flag
      plot_parms.wforms_err = zeros(parms.ntpoints,parms.nconds,2);
    else
      plot_parms.wforms_err = [];
    end;
    plot_sig_ranges = [];
    for c=1:parms.nconds
      tmp_results = resamp_results(a,c);
      if parms.ranges_flag & ~isempty(tmp_results.ranges_sig)
        plot_sig_ranges(c).onset = plot_parms.time(tmp_results.ranges_sig.t_onset);
        plot_sig_ranges(c).offset = plot_parms.time(tmp_results.ranges_sig.t_offset);
      else
        plot_sig_ranges(c).onset = [];
        plot_sig_ranges(c).offset = [];
      end;
      wforms = cat(2,wforms,tmp_results.wforms_mean);
      if parms.wforms_err_flag
        plot_parms.wforms_err(:,c,:) = ...
          reshape(tmp_results.wforms_bs_ci,[parms.ntpoints,1,2]);
      end;
    end;
    plot_parms.outstem = sprintf('%s_%s_diff_wforms',...
      parms.outstem,parms.area_names{a});
    if parms.title_flag
      plot_parms.title = sprintf('group %s response',parms.area_names{a});
    end;
    if parms.legend_flag
      plot_parms.condnames = parms.condnames;
    end;
    if parms.ylabel_flag
      plot_parms.ylabel = sprintf('Difference of Response Amplitude (%s)',...
          parms.units);
    end;
    plot_parms.ylim = parms.ylim_wform_diff;
    plot_parms.ranges = plot_sig_ranges; % use plot_sig_ranges to show sig ranges
    plot_parms.ranges_color = [];
    plot_parms.points = [];
    plot_parms.colors = parms.cond_colors;
    plot_args = mmil_parms2args(plot_parms);
    figure(parms.f); parms.f = parms.f + 1;
    ts_plot_wforms(wforms,plot_args{:});    
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_responses(parms,resamp_results)
  % plot peak latency vs. condition value
  if parms.peak_flag
    parms = plot_peak_latency(parms,resamp_results);
    parms = plot_peak_latency_diff(parms,resamp_results);
  end;
  % plot onset latency vs. condition value
  if parms.resamp_parms.onset_flag
    parms = plot_onset_latency(parms,resamp_results);
    parms = plot_onset_latency_diff(parms,resamp_results);
  end;
  % plot peak amplitude vs. condition value
  if parms.peak_flag
    for normflag=0:1
      parms = plot_peak_amplitude(parms,resamp_results,normflag);
      parms = plot_peak_amplitude_diff(parms,resamp_results,normflag);
    end;
  end;
  % plot area under the curve vs. condition value
  if parms.resamp_parms.auc_flag
    for normflag=0:1
      parms = plot_area_under_curve(parms,resamp_results,normflag);
      parms = plot_area_under_curve_diff(parms,resamp_results,normflag);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_peak_latency(parms,resamp_results)
  rplot_parms = parms.rplot_parms;
  rplot_parms.outstem = sprintf('%s_peak_latency',parms.outstem);
  if parms.title_flag
    rplot_parms.title = sprintf('peak latency vs. %s',parms.cond_varname);
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = 'Peak Latency (msec)';
  end;
  if parms.legend_flag
    rplot_parms.roinames = parms.outfix_list;
  end;
  rplot_parms.ylim = parms.ylim_latency;
  rplot_parms.normflag = 0;
  for a=1:parms.nareas
    area = parms.area_names{a};
    rplot_parms.outstem = sprintf('%s_%s_peak_latency',parms.outstem,area);
    if parms.title_flag
      rplot_parms.title = sprintf('%s peak latency vs. %s',...
        area,parms.cond_varname);
    end;
    responses = zeros(2,parms.nconds);
    rplot_parms.responses_err = zeros(2,parms.nconds,2);
    for c=1:parms.nconds
      responses(1,c) = resamp_results(a,c).peak1.latency;
      rplot_parms.responses_err(1,c,:) = resamp_results(a,c).peak1.latency_ci;
      responses(2,c) = resamp_results(a,c).peak2.latency;
      rplot_parms.responses_err(2,c,:) = resamp_results(a,c).peak2.latency_ci;
    end;
    rplot_parms.colors = parms.pair_colors;
    rplot_args = mmil_parms2args(rplot_parms);
    figure(parms.f); parms.f = parms.f + 1;
    ts_plot_responses(responses,rplot_args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_peak_latency_diff(parms,resamp_results)
  rplot_parms = parms.rplot_parms;
  rplot_parms.outstem = sprintf('%s_peak_latency_diff',parms.outstem);
  if parms.title_flag
    rplot_parms.title = ...
      sprintf('difference of peak latency vs. %s',parms.cond_varname);
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = 'Difference of Peak Latency (msec)';
  end;
  if parms.legend_flag
    rplot_parms.roinames = parms.area_names;
  end;
  rplot_parms.ylim = parms.ylim_latency_diff;
  rplot_parms.normflag = 0;
  responses = zeros(parms.nareas,parms.nconds);
  rplot_parms.responses_err = zeros(parms.nareas,parms.nconds,2);
  for a=1:parms.nareas
    for c=1:parms.nconds
      responses(a,c) = resamp_results(a,c).peak.latency;
      rplot_parms.responses_err(a,c,:) = resamp_results(a,c).peak.latency_ci;
    end;
  end;
  rplot_parms.colors = parms.area_colors;
  rplot_args = mmil_parms2args(rplot_parms);
  figure(parms.f); parms.f = parms.f + 1;
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_onset_latency(parms,resamp_results)
  rplot_parms = parms.rplot_parms;
  if parms.ylabel_flag
    rplot_parms.ylabel = 'Onset Latency (msec)';
  end;
  if parms.legend_flag
    rplot_parms.roinames = parms.outfix_list;
  end;
  rplot_parms.ylim = parms.ylim_latency;
  rplot_parms.normflag = 0;
  for a=1:parms.nareas
    area = parms.area_names{a};
    rplot_parms.outstem = sprintf('%s_%s_onset_latency',parms.outstem,area);
    if parms.title_flag
      rplot_parms.title = sprintf('%s onset latency vs. %s',...
        area,parms.cond_varname);
    end;
    responses = zeros(2,parms.nconds);
    rplot_parms.responses_err = zeros(2,parms.nconds,2);
    for c=1:parms.nconds
      responses(1,c) = resamp_results(a,c).onset1.latency;
      rplot_parms.responses_err(1,c,:) = resamp_results(a,c).onset1.latency_ci;
      responses(2,c) = resamp_results(a,c).onset2.latency;
      rplot_parms.responses_err(2,c,:) = resamp_results(a,c).onset2.latency_ci;
    end;
    rplot_parms.colors = parms.pair_colors;
    rplot_args = mmil_parms2args(rplot_parms);
    figure(parms.f); parms.f = parms.f + 1;
    ts_plot_responses(responses,rplot_args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_onset_latency_diff(parms,resamp_results)
  rplot_parms = parms.rplot_parms;
  rplot_parms.outstem = sprintf('%s_onset_latency_diff',parms.outstem);
  if parms.title_flag
    rplot_parms.title = ...
      sprintf('difference of onset latency vs. %s',parms.cond_varname);
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = 'Difference of Onset Latency (msec)';
  end;
  if parms.legend_flag
    rplot_parms.roinames = parms.area_names;
  end;
  rplot_parms.ylim = parms.ylim_latency_diff;
  rplot_parms.normflag = 0;
  responses = zeros(parms.nareas,parms.nconds);
  rplot_parms.responses_err = zeros(parms.nareas,parms.nconds,2);
  for a=1:parms.nareas
    for c=1:parms.nconds
      responses(a,c) = resamp_results(a,c).onset.latency;
      rplot_parms.responses_err(a,c,:) = resamp_results(a,c).onset.latency_ci;
    end;
  end;
  rplot_parms.colors = parms.area_colors;
  rplot_args = mmil_parms2args(rplot_parms);
  figure(parms.f); parms.f = parms.f + 1;
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_peak_amplitude(parms,resamp_results,normflag)
  rplot_parms = parms.rplot_parms;
  if normflag
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Normalized Peak Amplitude (%s)',...
        parms.norm_units);
    end;
    rplot_parms.ylim = parms.ylim_peak_norm;
    rplot_parms.normflag = 2;
  else
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Peak Amplitude (%s)',...
        parms.units);
    end;
    rplot_parms.ylim = parms.ylim_peak_amp;
    rplot_parms.normflag = 0;
  end;
  if parms.legend_flag
    rplot_parms.roinames = parms.outfix_list;
  end;
  for a=1:parms.nareas
    area = parms.area_names{a};
    if normflag
      rplot_parms.outstem = sprintf('%s_%s_peak_amplitude_norm',parms.outstem,area);
      if parms.title_flag
        rplot_parms.title = ...
          sprintf('%s normalized peak amplitude vs. %s',area,parms.cond_varname);
      end;
    else
      rplot_parms.outstem = sprintf('%s_%s_peak_amplitude',parms.outstem,area);
      if parms.title_flag
        rplot_parms.title = ...
          sprintf('%s peak amplitude vs. %s',area,parms.cond_varname);
      end;
    end;
    responses = zeros(2,parms.nconds);
    rplot_parms.responses_err = zeros(2,parms.nconds,2);
    for c=1:parms.nconds
      responses(1,c) = resamp_results(a,c).peak1.amplitude;
      rplot_parms.responses_err(1,c,:) = resamp_results(a,c).peak1.amplitude_ci;
      responses(2,c) = resamp_results(a,c).peak2.amplitude;
      rplot_parms.responses_err(2,c,:) = resamp_results(a,c).peak2.amplitude_ci;
    end;
    rplot_parms.colors = parms.pair_colors;
    rplot_args = mmil_parms2args(rplot_parms);
    figure(parms.f); parms.f = parms.f + 1;
    ts_plot_responses(responses,rplot_args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_peak_amplitude_diff(parms,resamp_results,normflag)
  rplot_parms = parms.rplot_parms;
  if normflag
    rplot_parms.outstem = sprintf('%s_peak_amplitude_diff_norm',parms.outstem);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('normalized difference of peak amplitude vs. %s',parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Normalized Difference of Peak Amplitude (%s)',...
        parms.norm_units);
    end;
    rplot_parms.ylim = parms.ylim_peak_norm;
    rplot_parms.normflag = 2;
  else
    rplot_parms.outstem = sprintf('%s_peak_amplitude_diff',parms.outstem);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('difference of peak amplitude vs. %s',parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Difference of Peak Amplitude (%s)',...
        parms.units);
    end;
    rplot_parms.ylim = parms.ylim_peak_amp;
    rplot_parms.normflag = 0;
  end;
  if parms.legend_flag
    rplot_parms.roinames = parms.area_names;
  end;
  responses = zeros(parms.nareas,parms.nconds);
  rplot_parms.responses_err = zeros(parms.nareas,parms.nconds,2);
  for a=1:parms.nareas
    for c=1:parms.nconds
      responses(a,c) = resamp_results(a,c).peak.amplitude;
      rplot_parms.responses_err(a,c,:) = resamp_results(a,c).peak.amplitude_ci;
    end;
  end;
  rplot_parms.colors = parms.area_colors;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_area_under_curve(parms,resamp_results,normflag)
  rplot_parms = parms.rplot_parms;
  if normflag
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Normalized Area Under Curve (%s)',...
        parms.norm_units);
    end;
    rplot_parms.ylim = parms.ylim_auc_norm;
    rplot_parms.normflag = 2;
  else
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Area Under Curve (%s)',...
        parms.units);
    end;
    rplot_parms.ylim = parms.ylim_auc;
    rplot_parms.normflag = 0;
  end;
  if parms.legend_flag
    rplot_parms.roinames = parms.outfix_list;
  end;
  for a=1:parms.nareas
    area = parms.area_names{a};
    if normflag
      rplot_parms.outstem = sprintf('%s_%s_auc_norm',parms.outstem,area);
      if parms.title_flag
        rplot_parms.title = ...
          sprintf('%s normalized area under curve vs. %s',area,parms.cond_varname);
      end;
    else
      rplot_parms.outstem = sprintf('%s_%s_auc',parms.outstem,area);
      if parms.title_flag
        rplot_parms.title = ...
          sprintf('%s area under curve vs. %s',area,parms.cond_varname);
      end;
    end;
    responses = zeros(2,parms.nconds);
    rplot_parms.responses_err = zeros(2,parms.nconds,2);
    for c=1:parms.nconds
      responses(1,c) = resamp_results(a,c).auc1.mean;
      rplot_parms.responses_err(1,c,:) = resamp_results(a,c).auc1.ci;
      responses(2,c) = resamp_results(a,c).auc2.mean;
      rplot_parms.responses_err(2,c,:) = resamp_results(a,c).auc2.ci;
    end;
    rplot_parms.colors = parms.pair_colors;
    rplot_args = mmil_parms2args(rplot_parms);
    figure(parms.f); parms.f = parms.f + 1;
    ts_plot_responses(responses,rplot_args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_area_under_curve_diff(parms,resamp_results,normflag)
  rplot_parms = parms.rplot_parms;
  if normflag
    rplot_parms.outstem = sprintf('%s_auc_diff_norm',parms.outstem);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('normalized difference of area under curve vs. %s',parms.cond_varname);
    end;
    if parms.ylabel_flag
     rplot_parms.ylabel = sprintf('Normalized Difference of Area Under Curve (%s)',...
        parms.norm_units);
    end;
    rplot_parms.ylim = parms.ylim_auc_norm_diff;
    rplot_parms.normflag = 2;
  else
    rplot_parms.outstem = sprintf('%s_auc_diff',parms.outstem);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('difference of area under curve vs. %s',parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Difference of Area Under Curve (%s)',...
        parms.units);
    end;
    rplot_parms.ylim = parms.ylim_auc_diff;
    rplot_parms.normflag = 0;
  end;
  if parms.legend_flag
    rplot_parms.roinames = parms.area_names;
  end;
  responses = zeros(parms.nareas,parms.nconds);
  rplot_parms.responses_err = zeros(parms.nareas,parms.nconds,2);
  for a=1:parms.nareas
    for c=1:parms.nconds
      responses(a,c) = resamp_results(a,c).auc.mean;
      rplot_parms.responses_err(a,c,:) = resamp_results(a,c).auc.ci;
    end;
  end;
  rplot_parms.colors = parms.area_colors;
  rplot_args = mmil_parms2args(rplot_parms);
  figure(parms.f); parms.f = parms.f + 1;
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_results(parms,resamp_results)
  outstem = [parms.outdir '/' parms.outstem];
  ts_summarize_resamp_results(resamp_results,'fstem_out',outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
