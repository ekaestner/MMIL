function resamp_results = MEG_MMIL_RCSE_GroupAvg_Err_Resamp(rootdir,varargin)
%function resamp_results = MEG_MMIL_RCSE_GroupAvg_Err_Resamp(rootdir,[options])
%
% Purpose: Analyze RCSE waveforms for a group of subjects
%   using bootstrap resampling
%
% Required Input:
%   rootdir: root directory containing RCSE_GroupAvg files
%
% Optional Input ('key',arg pairs):
%   'outdir': output directory
%     {default = [pwd '/RCSE_GroupAvg_Err_Resamp']}
%   'outstem': output file stem
%     {default = 'RCSE_Err'}
%   'prefix': RCSE prefix
%     {default = 'RCSE_prior'}
%
% Output:
%   resamp_result: struct array with resampling results
%     including means, sem, confidence intervals, etc.
%
% Created:  04/15/12 by Don Hagler
% Last Mod: 12/09/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: document optional parameters
%% todo: make more general
%%  (e.g. set_parms: conditions, xlabel; rplot_parms.title)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
resamp_results = [];

parms = check_input(rootdir,varargin);

fprintf('%s: loading RCSE GroupAvg results for %s...\n',mfilename,parms.outstem);
[results,parms] = load_results(parms);

parms = set_parms(parms);

fprintf('%s: resampling var, fit, and err waveforms for %s...\n',...
  mfilename,parms.outstem);
resamp_results = resamp_wforms(parms,results);

if parms.plot_wforms_flag
  fprintf('%s: plotting waveforms for %s...\n',mfilename,parms.outstem);
  parms = plot_wforms(parms,resamp_results);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(rootdir,options)
  parms_filter = {...
    'rootdir',rootdir,[],...
    'outdir',[pwd '/RCSE_GroupAvg_Err_Resamp'],[],...
    'outstem','RCSE_Err',[],...
    'prefix','RCSE_prior',[],...
    'condnames',{'var','fit','err'},[],...
    'niters',2000,[],...
    'ylim_wform',[0,1],[],...
    'xlim_wform',[-100,300],[],...
    'plot_wforms_flag',true,[false,true],...
    'cond_colors',{'k','b','r'},[],...
    'onset_kappa',3,[1,100],...
    'onset_baseline_range',[-100,0],[],...
    'wforms_err_flag',true,[false true],...
    'fill_err_flag',true,[false true],...
    'peak_pol',-1,[-1,1],...
    'peak_range',[0,200],[],...
    'peak_mindiff',0.25,[],...
    'resfact',1,[],...
    'norm_units','a.u.',[],...
    'wforms_fig_size',[],[],...
    'wforms_legend_loc','SouthEast',[],...
    'wforms_linewidth',1,[],...
    'wforms_zero_line_flag',false,[false,true],...
    'eps_flag',true,[false true],...
    'tif_dpi',300,[10,10000],...
    'fontname','Arial',[],...
    'fontsize',12,[],...
    'onset_flag',false,[false true],...
    'peak_flag',false,[false true],...
    'auc_flag',false,[false true],...
    'ylabel_flag',true,[false true],...
    'title_flag',true,[false true],...
    'legend_flag',true,[false true],...
    'markersize',6,[],...
    'axes_linewidth',1,[],...
    'chan_type','grad',{'grad','mag','eeg'},...
  };
  parms = mmil_args2parms(options,parms_filter);
  if ~exist(parms.rootdir,'dir')
    error('rootdir %s not found',parms.rootdir);
  end;
  parms.f = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_parms(parms)
  % parameters for resampling
  parms.resamp_parms = [];
  parms.resamp_parms.peak_flag = parms.peak_flag;
  parms.resamp_parms.onset_flag = parms.onset_flag;
  parms.resamp_parms.auc_flag = parms.auc_flag;
  parms.resamp_parms.niters = parms.niters;
  parms.resamp_parms.randseed_flag = 0;
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
  parms.plot_parms.fill_alpha = 0.4;
  parms.plot_parms.errbar_interval = 5;
  parms.plot_parms.visible_flag = 0;
  parms.plot_parms.relative_err_flag = 0;
  parms.plot_parms.zero_line_flag = parms.wforms_zero_line_flag;
  parms.plot_parms.eps_flag = parms.eps_flag;
  parms.plot_parms.fig_size = parms.wforms_fig_size;
  parms.plot_parms.fontname = parms.fontname;
  parms.plot_parms.fontsize = parms.fontsize;
  parms.plot_parms.tif_dpi = parms.tif_dpi;
  parms.plot_parms.time = parms.time;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results,parms] = load_results(parms);
  fname = sprintf('%s/%s_GroupAvg_results.mat',parms.rootdir,parms.prefix);
  tmp = load(fname);
  results = tmp.results;
  parms.nsubs = results.N;
  parms.ntpoints = length(results.time);
  parms.nareas = results.nareas;
  parms.nconds = results.nconditions;
  parms.time = results.time*1000;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resamp_results = resamp_wforms(parms,results)
  clear resamp_results;
  resamp_args = mmil_parms2args(parms.resamp_parms);
  for c=1:length(parms.condnames)
    cond = parms.condnames{c};
    % compile waveforms
    wforms = zeros(parms.ntpoints,parms.nsubs);
    for s=1:parms.nsubs
      fname = [cond '_' parms.chan_type];
      wforms(:,s) = results.subjdata(s).(fname);
    end;
    % resample waveforms
    tmp_results = ts_bootstrap_wforms(wforms,resamp_args{:});
    tmp_results.cond = cond;
    resamp_results(c) = tmp_results;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = plot_wforms(parms,resamp_results)
  % set plotting parameters
  plot_parms = parms.plot_parms;
  plot_parms.outstem = sprintf('%s_wforms',parms.outstem);
  if parms.title_flag
    plot_parms.title = sprintf('group average var, fit, and err');
  end;
  if parms.legend_flag
    plot_parms.condnames = parms.condnames;
  end;
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
  end;
  % plot waveforms
  figure(parms.f); parms.f = parms.f + 1;
  plot_args = mmil_parms2args(plot_parms);
  ts_plot_wforms(wforms,plot_args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

