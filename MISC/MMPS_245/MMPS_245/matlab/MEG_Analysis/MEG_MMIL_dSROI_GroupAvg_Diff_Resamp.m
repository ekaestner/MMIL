function resamp_results = MEG_MMIL_dSROI_GroupAvg_Diff_Resamp(fname_sets,varargin)
%function resamp_results = MEG_MMIL_dSROI_GroupAvg_Diff_Resamp(fname_sets,[options])
%
% Purpose: Compare two sets of dSPM-ROI waveforms across subjects
%   using bootstrap resampling
%
% Required Input:
%   fname_sets: cell array with two elements, each containing cell array
%     of dSROI_GroupAvg results mat files to be averaged together
%
% Optional Input ('key',arg pairs):
%  'outdir': output directory
%    {default = [pwd '/dSROI_GroupAvg_Diff_Resamp']}
%  'outstem': output file stem
%    {default = 'dSROI_Diff'}
%  'outfix_list': cell array of outfixes used for legends
%    must have 2 entries, corresponding to fname_sets
%    {default = {'A','B'}}
%  'roinames': cell array of ROIs to use
%    if empty, use all ROIs found in dSROI results
%    {default = []}
%  'forceflag': [0|1] overwrite existing output
%    {default: 0}
%
% Output:
%   resamp_result: struct array with resampling results
%     including means, sem, confidence intervals, etc.
%
% Created:  03/02/12 by Don Hagler
% Last Mod: 12/09/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: add more options to help (e.g. fname_sets_hemi)
%% todo: separate, general function for plotting resamp_results
%%       compatible with both this function and MEG_MMIL_GroupRCSE_Diff_Resamp
%% todo: plot histograms of onset latency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resamp_results = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = check_input(fname_sets,varargin);

mmil_mkdir(parms.outdir);
fname_results = [parms.outdir '/' parms.outstem '_GroupAvg_results.mat'];
if ~exist(fname_results,'file') || parms.forceflag
  % load previously calculated GroupAvg results for each condition
  %   (which contain waveforms for each subject)
  fprintf('%s: loading dSPM-ROI GroupAvg results for %s...\n',...
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
    fprintf('%s: loading dSPM-ROI GroupAvg results for %s...\n',...
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

function parms = check_input(fname_sets,options)
  parms_filter = {...
    'fname_sets',fname_sets,[],...
    'outdir',[pwd '/dSROI_GroupAvg_Diff_Resamp'],[],...
    'outstem','dSROI_Diff',[],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'roinames',[],[],...
    'forceflag',false,[false true],...  
    'condnames',{'LC','MC','HC'},[],...
    'conditions',[0.15,0.71,0.95],[],...
    'cond_varname','contrast',[],...
    'cond_label','Luminance Contrast',[],...
    'cond_lim',[0,1],[],...
    'cond_avg_flag',true,[false true],...
    'fname_sets_hemi',[],[],...
    'outfix_list',{'A','B'},[],...
    'niters',2000,[],...
    'bias_corr_flag',true,[false true],...
    'alpha',0.05,[0.001,0.1],...
    'ylim_wform',[-1,10],[],...
    'ylim_latency',[0,200],[],...
    'ylim_peak_norm',[-1,2],[],...
    'ylim_peak_amp',[-0.5,10],[],...
    'ylim_auc_norm',[-1,2],[],...
    'ylim_auc',[-10,50],[],...
    'ylim_wform_diff',[-5,5],[],...
    'ylim_latency_diff',[-100,100],[],...
    'ylim_peak_norm_diff',[-2,2],[],...
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
    'onset_kappa',4.5,[1,100],...
    'onset_baseline_range',[-100,0],[],...
    'roi_colors',{'b','g','r','m','k','y','c'},[],...
    'pair_colors',{'r','b'},[],...
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
  if length(parms.fname_sets)~=2
    error('fname_sets must have two sets of conditions for comparison');
  end;
  if ~isempty(parms.fname_sets_hemi) && length(parms.fname_sets_hemi)~=2
    error('fname_sets_hemi must match fname_sets');
  end;
  if length(parms.outfix_list)~=2
    error('length of outfix_list must match fname_sets');
  end;
  for p=1:length(parms.fname_sets)
    if ~iscell(parms.fname_sets{p})
      parms.fname_sets{p} = {parms.fname_sets{p}};
    end;
    nfiles = length(parms.fname_sets{p});
    if ~isempty(parms.fname_sets_hemi)
      if length(parms.fname_sets_hemi{p})~=nfiles
        error('fname_sets_hemi numbers of elements must match elements fname_sets');
      end;
      %% todo: check that all hemi values are either 1 or 2
    end;
  end;
  parms.pval_range = linspace(1/parms.niters,1,parms.niters);
  parms.f = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results,parms] = load_results(parms);
  % load results
  clear results;
  for p=1:length(parms.fname_sets)
    fnames = parms.fname_sets{p};
    nfiles = length(fnames);
    for f=1:length(fnames)
      tmp = load(fnames{f});
      results{p}{f} = tmp.results;
      if p==1 && f==1
        % get information from first results struct
        % NOTE: will have errors if other results have different N, roinames, etc.
        parms.time = tmp.results.time*1000;
        parms.nsubs = tmp.results.N;
        parms.ntpoints = length(parms.time);
        parms.full_roinames = tmp.results.roinames;
        parms.all_roinames = unique(regexprep(parms.full_roinames,'-[lr]h',''));
        if isempty(parms.roinames)
          parms.roinames = parms.all_roinames;
        else
          [tmp_roinames,ind_rois] = intersect(parms.roinames,parms.all_roinames);
          parms.roinames = parms.roinames(sort(ind_rois));
        end;
        parms.nrois = length(parms.roinames);
        if ~parms.nrois, error('no valid roinames'); end;
        parms.nconds = tmp.results.nconditions;
      end;      
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
  parms.rplot_parms.colors = parms.roi_colors;
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
  for r=1:parms.nrois
    roiname = parms.roinames{r};
    roistr = regexprep(roiname,'_',' ');
    plot_parms.title = roistr;
    ind_r = find(strcmp(parms.all_roinames,roiname));
    
    % compile waveforms
    wforms = zeros(parms.ntpoints,parms.nconds,parms.nsubs,2);
    for p=1:2
      nfiles = length(results{p});
      for f=1:nfiles
        tmp_results = results{p}{f};
        if ~isempty(parms.fname_sets_hemi)
          hemilist = parms.hemilist(parms.fname_sets_hemi{p}(f));
        else
          hemilist = parms.hemilist;
        end;
        for h=1:length(hemilist)
          % get indicies of rois with this name
          hemi = hemilist{h};
          full_roiname = sprintf('%s-%s',roiname,hemi);
          ind_roi = find(strcmp(parms.full_roinames,full_roiname));
          % compile data for all subjects
          for s=1:parms.nsubs
            tmp_wforms = squeeze(tmp_results.subjdata(s).wform(:,ind_roi,:));
            wforms(:,:,s,p) = wforms(:,:,s,p) + tmp_wforms;
          end;
        end;
      end;
      wforms(:,:,:,p) = wforms(:,:,:,p) / (nfiles * length(hemilist));
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
      tmp_results.roi = parms.roinames{r};
      tmp_results.cond = parms.condnames{c};
      resamp_results(r,c) = tmp_results;
    end;
  end;
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: rois not areas
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
        auc.sem = ts_nan_sem(auc.vals,0,2);
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
  % plot difference waveforms with all conds for each area
  parms = plot_wforms_diff(parms,resamp_results);
  % plot pairs for each cond and each area
  parms = plot_wforms_pair(parms,resamp_results);  
  % plot all conds for each area and each set
  parms = plot_wforms_conds(parms,resamp_results);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_wforms_diff(parms,resamp_results)
  % plot difference waveforms with all conds for each area
  fprintf('%s: plotting difference waveforms for %s...\n',...
    mfilename,parms.outstem);
  for r=1:parms.nrois
    roiname = parms.roinames{r};
    roistr = regexprep(roiname,'_',' ');
    % set plotting parameters
    plot_parms = parms.plot_parms;
    plot_parms.outstem = sprintf('%s_%s_wforms_diff',parms.outstem,roiname);
    if parms.title_flag
      plot_parms.title = sprintf('group %s response',roistr);
    end;
    if parms.legend_flag
      plot_parms.condnames = parms.condnames;
    end;
    plot_parms.ylim = parms.ylim_wform_diff;
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
      tmp_results = resamp_results(r,c);
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
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_wforms_pair(parms,resamp_results)
  % plot pairs for each cond and each area
  fprintf('%s: plotting paired waveforms for %s...\n',...
    mfilename,parms.outstem);
  for r=1:parms.nrois
    roiname = parms.roinames{r};
    roistr = regexprep(roiname,'_',' ');
    for c=1:parms.nconds
      cond = parms.condnames{c};
      tmp_results = resamp_results(r,c);
      % set plotting parameters
      plot_parms = parms.plot_parms;
      plot_parms.outstem = sprintf('%s_%s_%s_wforms_pair',...
        parms.outstem,roiname,cond);
      if parms.title_flag
        plot_parms.title = sprintf('group %s %s %s responses',...
          parms.outstem,roistr,cond);
      end;
      if parms.legend_flag
        plot_parms.condnames = parms.outfix_list;
      end;
      plot_parms.ylim = parms.ylim_wform;
      plot_parms.colors = parms.pair_colors;
      % compile data
      wforms = zeros(parms.ntpoints,2);
      wforms(:,1) = tmp_results.wforms1_mean;
      wforms(:,2) = tmp_results.wforms2_mean;
      if parms.wforms_err_flag
        plot_parms.wforms_err = zeros(parms.ntpoints,2,2);
        plot_parms.wforms_err(:,1,:) = ...
          reshape(tmp_results.wforms1_bs_ci,[parms.ntpoints,1,2]);
        plot_parms.wforms_err(:,2,:) = ...
          reshape(tmp_results.wforms2_bs_ci,[parms.ntpoints,1,2]);
      else
        plot_parms.wforms_err = [];
      end;
      plot_parms.points = [];
      if parms.mark_peak_flag
        plot_parms.points(1).amplitude = tmp_results.peak1.amplitude;
        plot_parms.points(1).latency = tmp_results.peak1.latency;
        plot_parms.points(1).type = parms.resamp_parms.peak_pol;
        plot_parms.points(2).amplitude = tmp_results.peak2.amplitude;
        plot_parms.points(2).latency = tmp_results.peak2.latency;
        plot_parms.points(2).type = parms.resamp_parms.peak_pol;
      end;
      plot_parms.ranges = [];
      if parms.mark_onset_flag
        plot_parms.ranges(1).onset = tmp_results.onset1.latency;
        plot_parms.ranges(1).offset = ...
          tmp_results.onset1.latency + parms.onset_marker_width;
        plot_parms.ranges(2).onset = tmp_results.onset2.latency;
        plot_parms.ranges(2).offset = ...
          tmp_results.onset2.latency + parms.onset_marker_width;
        plot_parms.ranges_color = [];
      elseif parms.ranges_flag & ~isempty(tmp_results.ranges_sig)
        plot_parms.ranges(1).onset = ...
          plot_parms.time(tmp_results.ranges_sig.t_onset);
        plot_parms.ranges(1).offset = ...
          plot_parms.time(tmp_results.ranges_sig.t_offset);
        plot_parms.ranges(2).onset = [];
        plot_parms.ranges(2).offset = [];
        plot_parms.ranges_color = parms.ranges_color;
      end;
      % plot waveforms
      figure(parms.f); parms.f = parms.f + 1;
      plot_args = mmil_parms2args(plot_parms);
      ts_plot_wforms(wforms,plot_args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_wforms_conds(parms,resamp_results)
  % plot all conds for each area and each set
  fprintf('%s: plotting condition waveforms for %s...\n',...
    mfilename,parms.outstem);
  for s=1:2
    outfix = parms.outfix_list{s};
    for r=1:parms.nrois
      roiname = parms.roinames{r};
      roistr = regexprep(roiname,'_',' ');
      % set plotting parameters
      plot_parms = parms.plot_parms;
      plot_parms.outstem = ...
        sprintf('%s_%s_%s_wforms',parms.outstem,roiname,outfix);
      if parms.title_flag
        plot_parms.title = sprintf('group %s %s response',roistr,outfix);
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
        tmp_results = resamp_results(r,c);
        tmp_field = sprintf('wforms%d_mean',s);        
        wforms(:,c) = tmp_results.(tmp_field);
        if parms.wforms_err_flag
          tmp_field = sprintf('wforms%d_bs_ci',s);
          plot_parms.wforms_err(:,c,:) = reshape(tmp_results.(tmp_field),...
            [parms.ntpoints,1,2]);
        end;
        if parms.mark_peak_flag
          tmp_field = sprintf('peak%d',s);
          plot_parms.points(c).amplitude = tmp_results.(tmp_field).amplitude;
          plot_parms.points(c).latency = tmp_results.(tmp_field).latency;
          plot_parms.points(c).type = parms.resamp_parms.peak_pol;
        end;
        if parms.mark_onset_flag
          tmp_field = sprintf('onset%d',s);
          plot_parms.ranges(c).onset = tmp_results.(tmp_field).latency;
          plot_parms.ranges(c).offset = ...
            tmp_results.(tmp_field).latency + parms.onset_marker_width;
        end;
      end;
      % plot waveforms
      figure(parms.f); parms.f = parms.f + 1;
      plot_args = mmil_parms2args(plot_parms);
      ts_plot_wforms(wforms,plot_args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_responses(parms,resamp_results)
  % plot peak latency vs. condition value
  if parms.peak_flag
    parms = plot_latency_diff(parms,resamp_results,'peak');
    parms = plot_latency_pair(parms,resamp_results,'peak');
  end;
  % plot onset latency vs. condition value
  if parms.resamp_parms.onset_flag
    parms = plot_latency_diff(parms,resamp_results,'onset');
    parms = plot_latency_pair(parms,resamp_results,'onset');
  end;
  % plot peak amplitude vs. condition value
  if parms.peak_flag
    for normflag=0:1
      parms = plot_peak_amplitude_diff(parms,resamp_results,normflag);
      parms = plot_peak_amplitude_pair(parms,resamp_results,normflag);
    end;
  end;
  % plot area under the curve vs. condition value
  if parms.resamp_parms.auc_flag
    for normflag=0:1
      parms = plot_area_under_curve_diff(parms,resamp_results,normflag);
      parms = plot_area_under_curve_pair(parms,resamp_results,normflag);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_latency_diff(parms,resamp_results,type_str)
  fprintf('%s: plotting difference in %s latency for %s...\n',...
    mfilename,type_str,parms.outstem);
  rplot_parms = parms.rplot_parms;
  rplot_parms.outstem = sprintf('%s_%s_latency_diff',parms.outstem,type_str);
  if parms.title_flag
    rplot_parms.title = ...
      sprintf('difference of %s latency vs. %s',type_str,parms.cond_varname);
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = 'Latency Difference (msec)';
  end;
  if parms.legend_flag
    rplot_parms.roinames = regexprep(parms.roinames,'_',' ');
  end;
  rplot_parms.ylim = parms.ylim_latency_diff;
  rplot_parms.normflag = 0;
  rplot_parms.responses_err = zeros(parms.nrois,parms.nconds,2);
  responses = zeros(parms.nrois,parms.nconds);
  for r=1:parms.nrois
    for c=1:parms.nconds
      responses(r,c) = resamp_results(r,c).(type_str).latency;
      rplot_parms.responses_err(r,c,:) = resamp_results(r,c).(type_str).latency_ci;
    end;
  end;
  figure(parms.f); parms.f = parms.f + 1;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_latency_pair(parms,resamp_results,type_str)
  fprintf('%s: plotting %s latency for %s...\n',...
    mfilename,type_str,parms.outstem);    
  for r=1:parms.nrois
    roiname = parms.roinames{r};
    roistr = regexprep(roiname,'_',' ');
    rplot_parms = parms.rplot_parms;
    rplot_parms.outstem = sprintf('%s_%s_%s_latency',...
      parms.outstem,roiname,type_str);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('%s %s latency vs. %s',roistr,type_str,parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = 'Latency (msec)';
    end;
    if parms.legend_flag
      rplot_parms.roinames = parms.outfix_list;
    end;
    rplot_parms.ylim = parms.ylim_latency;
    rplot_parms.colors = parms.pair_colors;
    rplot_parms.normflag = 0;
    rplot_parms.responses_err = zeros(2,parms.nconds,2);
    responses = zeros(2,parms.nconds);
    for c=1:parms.nconds
      tmp_results = resamp_results(r,c);
      tmp_field = [type_str '1'];
      responses(1,c) = tmp_results.(tmp_field).latency;
      rplot_parms.responses_err(1,c,:) = tmp_results.(tmp_field).latency_ci;
      tmp_field = [type_str '2'];
      responses(2,c) = tmp_results.(tmp_field).latency;
      rplot_parms.responses_err(2,c,:) = tmp_results.(tmp_field).latency_ci;
    end;
    figure(parms.f); parms.f = parms.f + 1;
    rplot_args = mmil_parms2args(rplot_parms);
    ts_plot_responses(responses,rplot_args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_peak_amplitude_diff(parms,resamp_results,normflag)
  rplot_parms = parms.rplot_parms;
  if normflag
    fprintf('%s: plotting normalized peak amplitude difference for %s...\n',...
      mfilename,parms.outstem);
    rplot_parms.outstem = sprintf('%s_peak_amplitude_diff_norm',parms.outstem);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('normalized difference of peak amplitude vs. %s',...
          parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = ...
        sprintf('Normalized Difference of Peak Amplitude (%s)',parms.norm_units);
    end;
    rplot_parms.ylim = parms.ylim_peak_norm_diff;
    rplot_parms.normflag = 2;
  else
    fprintf('%s: plotting peak amplitude difference for %s...\n',...
      mfilename,parms.outstem);
    rplot_parms.outstem = sprintf('%s_peak_amplitude_diff',parms.outstem);
    if parms.title_flag
      rplot_parms.title = sprintf('difference of peak amplitude vs. %s',...
        parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Difference of Peak Amplitude (%s)',...
        parms.units);
    end;
    rplot_parms.ylim = parms.ylim_peak_amp_diff;
    rplot_parms.normflag = 0;
  end;
  if parms.legend_flag
    rplot_parms.roinames = regexprep(parms.roinames,'_',' ');
  end;
  rplot_parms.responses_err = zeros(parms.nrois,parms.nconds,2);
  responses = zeros(parms.nrois,parms.nconds);
  for r=1:parms.nrois
    for c=1:parms.nconds
      responses(r,c) = resamp_results(r,c).peak.amplitude;
      rplot_parms.responses_err(r,c,:) = resamp_results(r,c).peak.amplitude_ci;
    end;
  end;
  figure(parms.f); parms.f = parms.f + 1;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_peak_amplitude_pair(parms,resamp_results,normflag)
  if normflag
    fprintf('%s: plotting normalized peak amplitude for %s...\n',...
      mfilename,parms.outstem);
  else
    fprintf('%s: plotting peak amplitude for %s...\n',...
      mfilename,parms.outstem);
  end;
  for r=1:parms.nrois
    roiname = parms.roinames{r};
    roistr = regexprep(roiname,'_',' ');
    rplot_parms = parms.rplot_parms;
    if normflag
      rplot_parms.outstem = sprintf('%s_%s_peak_amplitude_norm',...
        parms.outstem,roiname);
      if parms.title_flag
        rplot_parms.title = ...
          sprintf('%s normalized peak amplitude vs. %s',...
            roistr,parms.cond_varname);
      end;
      if parms.ylabel_flag
        rplot_parms.ylabel = sprintf('Normalized Peak Amplitude (%s)',...
          parms.norm_units);
      end;
      rplot_parms.ylim = parms.ylim_peak_norm;
      rplot_parms.normflag = 2;
    else
      rplot_parms.outstem = sprintf('%s_%s_peak_amplitude',...
        parms.outstem,roiname);
      if parms.title_flag
        rplot_parms.title = sprintf('%s peak amplitude vs. %s',...
          roistr,parms.cond_varname);
      end;
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
    rplot_parms.colors = parms.pair_colors;
    rplot_parms.responses_err = zeros(2,parms.nconds,2);
    responses = zeros(2,parms.nconds);
    for c=1:parms.nconds
      responses(1,c) = resamp_results(r,c).peak1.amplitude;
      rplot_parms.responses_err(1,c,:) = resamp_results(r,c).peak1.amplitude_ci;
      responses(2,c) = resamp_results(r,c).peak2.amplitude;
      rplot_parms.responses_err(2,c,:) = resamp_results(r,c).peak2.amplitude_ci;
    end;
    figure(parms.f); parms.f = parms.f + 1;
    rplot_args = mmil_parms2args(rplot_parms);
    ts_plot_responses(responses,rplot_args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_area_under_curve_diff(parms,resamp_results,normflag)
  rplot_parms = parms.rplot_parms;
  if normflag
    fprintf('%s: plotting normalized area under curve difference for %s...\n',...
      mfilename,parms.outstem);
    rplot_parms.outstem = sprintf('%s_auc_diff_norm',parms.outstem);
    if parms.title_flag
      rplot_parms.title = ...
        sprintf('normalized difference of area under curve vs. %s',...
          parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = ...
        sprintf('Normalized Difference of Area Under Curve (%s)',...
          parms.norm_units);
    end;
    rplot_parms.ylim = parms.ylim_auc_norm_diff;
    rplot_parms.normflag = 2;
  else
    fprintf('%s: plotting area under curve difference for %s...\n',...
      mfilename,parms.outstem);
    rplot_parms.outstem = sprintf('%s_auc',parms.outstem);
    if parms.title_flag
      rplot_parms.title = sprintf('difference of area under curve vs. %s',...
        parms.cond_varname);
    end;
    if parms.ylabel_flag
      rplot_parms.ylabel = sprintf('Difference of Area Under Curve (%s)',...
        parms.units);
    end;
    rplot_parms.ylim = parms.ylim_auc_diff;
    rplot_parms.normflag = 0;
  end;
  if parms.legend_flag
    rplot_parms.roinames = regexprep(parms.roinames,'_',' ');
  end;
  rplot_parms.responses_err = zeros(parms.nrois,parms.nconds,2);
  responses = zeros(parms.nrois,parms.nconds);
  for r=1:parms.nrois
    for c=1:parms.nconds
      responses(r,c) = resamp_results(r,c).auc.mean;
      rplot_parms.responses_err(r,c,:) = resamp_results(r,c).auc.ci;
    end;
  end;
  figure(parms.f); parms.f = parms.f + 1;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_area_under_curve_pair(parms,resamp_results,normflag)
  if normflag
    fprintf('%s: plotting normalized area under curve for %s...\n',...
      mfilename,parms.outstem);
  else
    fprintf('%s: plotting area under curve for %s...\n',...
      mfilename,parms.outstem);
  end;
  for r=1:parms.nrois
    roiname = parms.roinames{r};
    roistr = regexprep(roiname,'_',' ');
    rplot_parms = parms.rplot_parms;
    if normflag
      rplot_parms.outstem = sprintf('%s_%s_auc_norm',parms.outstem,roiname);
      if parms.title_flag
        rplot_parms.title = ...
          sprintf('%s normalized area under curve vs. %s',...
            roistr,parms.cond_varname);
      end;
      if parms.ylabel_flag
        rplot_parms.ylabel = sprintf('Normalized Area Under Curve (%s)',...
          parms.norm_units);
      end;
      rplot_parms.ylim = parms.ylim_auc_norm;
      rplot_parms.normflag = 2;
    else
      rplot_parms.outstem = sprintf('%s_%s_auc',parms.outstem,roiname);
      if parms.title_flag
        rplot_parms.title = sprintf('%s area under curve vs. %s',...
          roistr,parms.cond_varname);
      end;
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
    rplot_parms.colors = parms.pair_colors;
    rplot_parms.responses_err = zeros(2,parms.nconds,2);
    responses = zeros(2,parms.nconds);
    for c=1:parms.nconds
      responses(1,c) = resamp_results(r,c).auc1.mean;
      rplot_parms.responses_err(1,c,:) = resamp_results(r,c).auc1.ci;
      responses(2,c) = resamp_results(r,c).auc2.mean;
      rplot_parms.responses_err(2,c,:) = resamp_results(r,c).auc2.ci;
    end;
    figure(parms.f); parms.f = parms.f + 1;
    rplot_args = mmil_parms2args(rplot_parms);
    ts_plot_responses(responses,rplot_args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_results(parms,resamp_results)
  outstem = [parms.outdir '/' parms.outstem];
  ts_summarize_resamp_results(resamp_results,'fstem_out',outstem);
return;

