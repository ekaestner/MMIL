function resamp_results = MEG_MMIL_RCSE_GroupAvg_DiffAreas_Resamp(rootdir,varargin)
%function resamp_results = MEG_MMIL_RCSE_GroupAvg_DiffAreas_Resamp(rootdir,[options])
%
% Purpose: Compare RCSE waveforms for different visual areas across subjects
%   using bootstrap resampling
%
% Required Input:
%   rootdir: root directory containing RCSE_GroupAvg files
%
% Optional Input ('key',arg pairs):
%   'outdir': output directory
%     {default = [pwd '/RCSE_GroupAvg_DiffAreas_Resamp']}
%   'outstem': output file stem
%     {default = 'RCSE_DiffAreas'}
%   'prefix': RCSE prefix
%     {default = 'RCSE_prior'}
%
% Output:
%   resamp_result: struct array with resampling results
%     including means, sem, confidence intervals, etc.
%
% Created:  03/08/12 by Don Hagler
% Last Mod: 12/09/13 by Don Hagler
%

% Purpose: Compare RCSE waveforms for different visual areas across subjects
%   using bootstrap resampling and permutation testing

%% todo: reorganize into subfunctions (like MEG_MMIL_RCSE_GroupAvg_Resamp)

%% todo: call ts_permute_wforms if ranges_flag = 1

if ~mmil_check_nargs(nargin,2), return; end;
clear resamp_results;
parms_filter = {...
  'rootdir',rootdir,[],...
  'outdir',[pwd '/RCSE_GroupAvg_DiffAreas_Resamp'],[],...
  'outstem','RCSE_DiffAreas',[],...
  'prefix','RCSE_prior',[],...
  'area_names',{'V1','V2','V3'},[],...
  'condnames',{'LC','MC','HC'},[],...
  'niters',2000,[],...
  'ylim_wform',[-20,10],[],...
  'ylim_wform_diff',[-5,5],[],...
  'ylim_latency_diff',[-100,100],[],...
  'ylim_peak_norm_diff',[-5,5],[],...
  'ylim_peak_amp_diff',[-10,10],[],...
  'ylim_auc_norm_diff',[-5,5],[],...
  'ylim_auc_diff',[-50,50],[],...
  'xlim_wform',[-100,300],[],...
  'ranges_color','y',[],...
  'plot_wforms_flag',true,[false,true],...
  'mark_peak_flag',true,[false,true],...
  'mark_onset_flag',true,[false,true],...
  'onset_marker_width',2,[1,1000],...
  'onset_kappa',3,[1,100],...
  'onset_baseline_range',[-100,0],[],...
  'area_colors',{'b','g','r'},[],...
  'comp_colors',{'r','k','g'},[],...
  'cond_colors',{'k','b','r'},[],...
  'fill_err_flag',true,[false true],...
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
  'ranges_flag',false,[false true],...
  'onset_flag',true,[false true],...
  'peak_flag',true,[false true],...
  'auc_flag',true,[false true],...
  'ylabel_flag',true,[false true],...
  'title_flag',true,[false true],...
  'legend_flag',true,[false true],...
  'markersize',6,[],...
  'axes_linewidth',1,[],...
  'csv_flag',true,[false true],...
  'alpha',0.05,[0.001,0.1],...
};
parms = mmil_args2parms(varargin,parms_filter);

if parms.normflag
  parms.units = parms.norm_units;
end;

if ~exist(parms.rootdir,'dir')
  error('rootdir %s not found',parms.rootdir);
end;

if ~parms.onset_flag
  parms.mark_onset_flag = false;
end;
if ~parms.peak_flag
  parms.mark_peak_flag = false;
end;

% parameters for resampling
resamp_parms = [];
resamp_parms.peak_flag = parms.peak_flag;
resamp_parms.onset_flag = parms.onset_flag;
resamp_parms.auc_flag = parms.auc_flag;
%resamp_parms.ranges_flag = parms.ranges_flag;
resamp_parms.niters = parms.niters;
resamp_parms.randseed_flag = 0;
resamp_parms.alpha = parms.alpha;
% for peaks
resamp_parms.peak_pol = parms.peak_pol;
resamp_parms.peak_range = parms.peak_range;
resamp_parms.peak_mindiff = parms.peak_mindiff;
resamp_parms.resfact = parms.resfact;
% for onset
resamp_parms.onset_baseline_range = parms.onset_baseline_range;
resamp_parms.onset_baseline_collapse_flag = 2;
resamp_parms.onset_minimum = 40;
resamp_parms.onset_polarity = -1;
resamp_parms.onset_baseline_flag = 2;
resamp_parms.onset_method = 'quartiles';
resamp_parms.onset_kappa = parms.onset_kappa;

% parameters for plotting waveforms
plot_parms = [];
if parms.ylabel_flag
  plot_parms.ylabel = sprintf('Response Amplitude (%s)',...
      parms.units);
end;
plot_parms.xlim = parms.xlim_wform;
plot_parms.legend_loc = parms.wforms_legend_loc;
plot_parms.outdir = parms.outdir;
plot_parms.baseline_flag = 0;
plot_parms.linewidth = parms.wforms_linewidth;
plot_parms.axes_linewidth = parms.axes_linewidth;
plot_parms.fill_err_flag = parms.fill_err_flag;
plot_parms.fill_alpha = 0.4;
plot_parms.errbar_interval = 5;
plot_parms.visible_flag = 0;
plot_parms.relative_err_flag = 0;
plot_parms.zero_line_flag = parms.wforms_zero_line_flag;
plot_parms.eps_flag = parms.eps_flag;
plot_parms.fig_size = parms.wforms_fig_size;
plot_parms.fontname = parms.fontname;
plot_parms.fontsize = parms.fontsize;
plot_parms.tif_dpi = parms.tif_dpi;

% for plotting responses
rplot_parms = [];
rplot_parms.outdir = parms.outdir;
rplot_parms.conditions = [0.15,0.71,0.95];
rplot_parms.xlabel = 'Luminance Contrast';
rplot_parms.xlim = [0,1];
rplot_parms.legend_loc = parms.resps_legend_loc;
rplot_parms.linewidth = parms.resps_linewidth;
rplot_parms.axes_linewidth = parms.axes_linewidth;
rplot_parms.visible_flag = 0;
rplot_parms.relative_err_flag = 0;
rplot_parms.zero_line_flag = parms.resps_zero_line_flag;
rplot_parms.eps_flag = parms.eps_flag;
rplot_parms.colors = parms.comp_colors;
rplot_parms.fig_size = parms.resps_fig_size;
rplot_parms.fontname = parms.fontname;
rplot_parms.fontsize = parms.fontsize;
rplot_parms.tif_dpi = parms.tif_dpi;
rplot_parms.markersize = parms.markersize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load results
fname = sprintf('%s/%s_GroupAvg_results.mat',parms.rootdir,parms.prefix);
results = load_results(fname);

nsubs = results.N;
ntpoints = length(results.time);
nareas = results.nareas;
nconds = results.nconditions;

plot_parms.time = results.time*1000;
resamp_parms.time = plot_parms.time;

resamp_args = mmil_parms2args(resamp_parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resample waveforms
aa = 1;
for a1=nareas:-1:1
  area1 = parms.area_names{a1};
  wforms1 = zeros(ntpoints,nconds,nsubs);
  for s=1:nsubs
    wforms1(:,:,s) = squeeze(results.subjdata(s).wform(:,a1,:));
  end;
  for a2=1:a1-1
    area2 = parms.area_names{a2};
    compstr = [area1 ' vs ' area2];
    wforms2 = zeros(ntpoints,nconds,nsubs);
    for s=1:nsubs
      wforms2(:,:,s) = squeeze(results.subjdata(s).wform(:,a2,:));
    end;
    wforms = cat(4,wforms1,wforms2);

    % normalize waveforms by max across time and conditions
    %  averaged across the pair of prefix+infix (to be unbiased)
    if parms.normflag
      max_vals = squeeze(max(squeeze(max(abs(wforms),[],1)),[],1));
      max_mat = repmat(reshape(max_vals,[1,1,nsubs,2]),[ntpoints,nconds,1,1]);
      wforms = wforms ./ max_mat;
    end;

    for c=1:nconds
      tmp_wforms = squeeze(wforms(:,c,:,:));
      tmp_results = ts_bootstrap_wforms(tmp_wforms,resamp_args{:});
      tmp_results.roi = compstr;
      tmp_results.cond = parms.condnames{c};
      resamp_results(aa,c) = tmp_results;
    end;
    aa = aa + 1;
  end;
end;
naa = size(resamp_results,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot waveforms

k = 1;
if parms.plot_wforms_flag
  aa = 1;
  for a1=nareas:-1:1
    area1 = parms.area_names{a1};
    for a2=1:a1-1
      area2 = parms.area_names{a2};
      infix = [area1 '_VS_' area2];
      compstr = [area1 ' vs ' area2];
      fprintf('%s: plotting difference waveforms for %s %s...\n',...
        mfilename,parms.outstem,compstr);
%% todo: remove unused variables
      wforms_mean = [];
      wforms_bs_ci = zeros(ntpoints,nconds,2);
      plot_sig_ranges = [];
      for c=1:nconds
        tmp_results = resamp_results(aa,c);
        if parms.ranges_flag & ~isempty(tmp_results.ranges_sig)
          fprintf('%s: significant range for %s, c=%d\n',...
            mfilename,compstr,c);
          plot_sig_ranges(c).onset = plot_parms.time(tmp_results.ranges_sig.t_onset);
          plot_sig_ranges(c).offset = plot_parms.time(tmp_results.ranges_sig.t_offset);
        else
          plot_sig_ranges(c).onset = [];
          plot_sig_ranges(c).offset = [];
        end;
        wforms_mean = cat(2,wforms_mean,tmp_results.wforms_mean);
        wforms_bs_ci(:,c,:) = reshape(tmp_results.wforms_bs_ci,[ntpoints,1,2]);
      end;

      figure(k); k=k+1;
      plot_parms.outstem = sprintf('%s_%s_diff_wforms',...
        parms.outstem,infix);
      if parms.title_flag
        plot_parms.title = sprintf('group %s response',compstr);
      end;
      if parms.wforms_err_flag
        plot_parms.wforms_err = wforms_bs_ci;
      else
        plot_parms.wforms_err = [];
      end;
      if parms.legend_flag
        plot_parms.condnames = parms.condnames;
      end;
      if parms.ylabel_flag
        plot_parms.ylabel = sprintf('Difference of Response Amplitude (%s)',...
            parms.units);
      end;
      plot_parms.ylim = parms.ylim_wform_diff;
      plot_parms.ranges = plot_sig_ranges; % show significant ranges
      plot_parms.ranges_color = [];
      plot_parms.points = [];
      plot_parms.colors = parms.cond_colors;
      plot_args = mmil_parms2args(plot_parms);
      ts_plot_wforms(wforms_mean,plot_args{:});

      for c=1:nconds
        figure(k); k=k+1;
        plot_parms.outstem = sprintf('%s_%s_cond%d_wforms',...
          parms.outstem,infix,c);
        if parms.title_flag
          plot_parms.title = sprintf('group %s response; contrast = %0.2f',...
            compstr,rplot_parms.conditions(c));
        end;
        tmp_wforms = ...
          cat(2,resamp_results(aa,c).wforms1_mean,resamp_results(aa,c).wforms2_mean);
        if parms.wforms_err_flag
          plot_parms.wforms_err = ...
            cat(2,resamp_results(aa,c).wforms1_bs_ci,resamp_results(aa,c).wforms2_bs_ci);
        end;
        if parms.legend_flag
          plot_parms.condnames = {area1,area2};
        end;
        if parms.ylabel_flag
          plot_parms.ylabel = sprintf('Response Amplitude (%s)',...
              parms.units);
        end;
        plot_parms.ylim = parms.ylim_wform;
        plot_ranges = [];
        if parms.mark_onset_flag
          plot_ranges(1).onset = resamp_results(aa,c).onset1.latency;
          plot_ranges(1).offset = resamp_results(aa,c).onset1.latency + parms.onset_marker_width;
          plot_ranges(2).onset = resamp_results(aa,c).onset2.latency;
          plot_ranges(2).offset = resamp_results(aa,c).onset2.latency + parms.onset_marker_width;
          plot_parms.ranges_color = [];
        elseif parms.ranges_flag & ~isempty(resamp_results(aa,c).ranges_sig)
          plot_ranges(1).onset = plot_parms.time(resamp_results(aa,c).ranges_sig.t_onset);
          plot_ranges(1).offset = plot_parms.time(resamp_results(aa,c).ranges_sig.t_offset);
          plot_ranges(2).onset = [];
          plot_ranges(2).offset = [];
          plot_parms.ranges_color = parms.ranges_color;
        end;
        plot_points = [];
        if parms.mark_peak_flag
          plot_points(1).amplitude = resamp_results(aa,c).peak1.amplitude;
          plot_points(1).latency = resamp_results(aa,c).peak1.latency;
          plot_points(1).type = resamp_parms.peak_pol;
          plot_points(2).amplitude = resamp_results(aa,c).peak2.amplitude;
          plot_points(2).latency = resamp_results(aa,c).peak2.latency;
          plot_points(2).type = resamp_parms.peak_pol;
        end;
        plot_parms.ranges = plot_ranges;
        plot_parms.points = plot_points;
        plot_parms.colors = parms.area_colors([a1,a2]);
        plot_args = mmil_parms2args(plot_parms);
        ts_plot_wforms(tmp_wforms,plot_args{:});
      end;
      aa = aa + 1;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot responses

if parms.legend_flag
  rplot_parms.roinames = cell(naa,1);
  aa = 1;
  for a1=nareas:-1:1
    for a2=1:a1-1
      rplot_parms.roinames{aa} = sprintf('%s vs %s',...
        parms.area_names{a1},parms.area_names{a2});
      aa = aa + 1;
    end;
  end;
end;

if parms.peak_flag
  % peak latency difference
  fprintf('%s: plotting difference in peak latency for %s...\n',...
    mfilename,parms.outstem);
  figure(k); k=k+1;
  rplot_parms.outstem = sprintf('%s_peak_latency_diff',parms.outstem);
  if parms.title_flag
    rplot_parms.title = 'difference of peak latency vs. contrast';
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = 'Difference of Peak Latency (msec)';
  end;
  rplot_parms.ylim = parms.ylim_latency_diff;
  rplot_parms.normflag = 0;

  responses = zeros(naa,nconds);
  rplot_parms.responses_err = zeros(naa,nconds,2);
  aa = 1;
  for a1=nareas:-1:1
    for a2=1:a1-1
      for c=1:nconds
        responses(aa,c) = resamp_results(aa,c).peak.latency;
        rplot_parms.responses_err(aa,c,:) = resamp_results(aa,c).peak.latency_ci;
      end;
      aa = aa + 1;
    end;
  end;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
end;

if resamp_parms.onset_flag
  % onset latency difference
  fprintf('%s: plotting onset latency difference for %s...\n',mfilename,parms.outstem);
  figure(k); k=k+1;
  rplot_parms.outstem = sprintf('%s_onset_latency_diff',parms.outstem);
  if parms.title_flag
    rplot_parms.title = 'difference of onset latency vs. contrast';
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = 'Difference of Onset Latency (msec)';
  end;
  rplot_parms.ylim = parms.ylim_latency_diff;
  rplot_parms.normflag = 0;
  responses = zeros(naa,nconds);
  rplot_parms.responses_err = zeros(nareas,nconds,2);
  aa = 1;
  for a1=nareas:-1:1
    for a2=1:a1-1
      for c=1:nconds
        responses(aa,c) = resamp_results(aa,c).onset.latency;
        rplot_parms.responses_err(aa,c,:) = resamp_results(aa,c).onset.latency_ci;
      end;
      aa = aa + 1;
    end;
  end;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
end;

if parms.peak_flag
  % normalized peak amplitude difference
  fprintf('%s: plotting normalized peak amplitude difference for %s...\n',...
    mfilename,parms.outstem);
  figure(k); k=k+1;
  rplot_parms.outstem = sprintf('%s_peak_amplitude_diff_norm',parms.outstem);
  if parms.title_flag
    rplot_parms.title = sprintf('normalized difference of peak amplitude vs. contrast');
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = sprintf('Normalized Difference of Peak Amplitude (%s)',...
      parms.norm_units);
  end;
  rplot_parms.ylim = parms.ylim_peak_norm_diff;
  rplot_parms.normflag = 2;
  responses = zeros(naa,nconds);
  rplot_parms.responses_err = zeros(nareas,nconds,2);
  aa = 1;
  for a1=nareas:-1:1
    for a2=1:a1-1
      for c=1:nconds
        responses(aa,c) = resamp_results(aa,c).peak.amplitude;
        rplot_parms.responses_err(aa,c,:) = resamp_results(aa,c).peak.amplitude_ci;
      end;
      aa = aa + 1;
    end;
  end;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});

  % raw peak amplitude difference
  figure(k); k=k+1;
  fprintf('%s: plotting peak amplitude difference for %s...\n',...
    mfilename,parms.outstem);
  rplot_parms.outstem = sprintf('%s_peak_amplitude_diff',parms.outstem);
  if parms.title_flag
    rplot_parms.title = sprintf('difference of peak amplitude vs. contrast');
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = sprintf('Difference of Peak Amplitude (%s)',...
      parms.units);
  end;
  rplot_parms.ylim = parms.ylim_peak_amp_diff;
  rplot_parms.normflag = 0;
  responses = zeros(naa,nconds);
  rplot_parms.responses_err = zeros(nareas,nconds,2);
  aa = 1;
  for a1=nareas:-1:1
    for a2=1:a1-1
      for c=1:nconds
        responses(aa,c) = resamp_results(aa,c).peak.amplitude;
        rplot_parms.responses_err(aa,c,:) = resamp_results(aa,c).peak.amplitude_ci;
      end;
      aa = aa + 1;
    end;
  end;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});

end;


if parms.auc_flag
  % normalized area under curve difference
  fprintf('%s: plotting normalized area under curve difference for %s...\n',...
    mfilename,parms.outstem);
  figure(k); k=k+1;
  rplot_parms.outstem = sprintf('%s_auc_diff_norm',parms.outstem);
  if parms.title_flag
    rplot_parms.title = sprintf('normalized difference of area under curve vs. contrast');
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = sprintf('Normalized Difference of Area Under Curve (%s)',...
      parms.norm_units);
  end;
  rplot_parms.ylim = parms.ylim_auc_norm_diff;
  rplot_parms.normflag = 2;
  responses = zeros(naa,nconds);
  rplot_parms.responses_err = zeros(nareas,nconds,2);
  aa = 1;
  for a1=nareas:-1:1
    for a2=1:a1-1
      for c=1:nconds
        responses(aa,c) = resamp_results(aa,c).auc.mean;
        rplot_parms.responses_err(aa,c,:) = resamp_results(aa,c).auc.ci;
      end;
      aa = aa + 1;
    end;
  end;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});

  % raw area under curve difference
  figure(k); k=k+1;
  fprintf('%s: plotting area under curve difference for %s...\n',...
    mfilename,parms.outstem);
  rplot_parms.outstem = sprintf('%s_auc_diff',parms.outstem);
  if parms.title_flag
    rplot_parms.title = sprintf('difference of area under curve vs. contrast');
  end;
  if parms.ylabel_flag
    rplot_parms.ylabel = sprintf('Difference of Area Under Curve (%s)',...
      parms.units);
  end;
  rplot_parms.ylim = parms.ylim_auc_diff;
  rplot_parms.normflag = 0;
  responses = zeros(naa,nconds);
  rplot_parms.responses_err = zeros(nareas,nconds,2);
  aa = 1;
  for a1=nareas:-1:1
    for a2=1:a1-1
      for c=1:nconds
        responses(aa,c) = resamp_results(aa,c).auc.mean;
        rplot_parms.responses_err(aa,c,:) = resamp_results(aa,c).auc.ci;
      end;
      aa = aa + 1;
    end;
  end;
  rplot_args = mmil_parms2args(rplot_parms);
  ts_plot_responses(responses,rplot_args{:});
end;

if parms.csv_flag
  fprintf('%s: summarize results in csv files for %s...\n',...
    mfilename,parms.outstem);
  summarize_results(parms,resamp_results);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_results(parms,resamp_results)
  outstem = [parms.outdir '/' parms.outstem];
  ts_summarize_resamp_results(resamp_results,'fstem_out',outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = load_results(fname)
  results = [];
  load(fname);
return;

