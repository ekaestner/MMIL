function rc_plot_traces_datafit(prefix,conditions,chantype,scale_max)
%function rc_plot_traces_datafit(prefix,[conditions],[chantype],[scale_max])
%
% Purpose: multiplot MEG/EEG sensor waveforms for RCSE data, fit, and err
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Input:
%   conditions: vector of condition numbers to loop over
%     {default = 1}
%   chantype: e.g. 'grad1', 'grad2', 'mag', 'eeg'
%     {default = 'grad1'}
%   scale_max: y-axis max scale value
%     if 0, auto-scale
%     {default = 0}
%
% Early Mod: 12/22/06 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%% todo: varargin
%% todo: rootdir

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('conditions','var') || isempty(conditions), conditions = 1; end;
if ~exist('chantype','var') || isempty(chantype), chantype = 'grad1'; end;
if ~exist('scale_max','var') || isempty(scale_max), scale_max = 0; end;
imgsize = 1000;

matfile = sprintf('matfiles/%s_parms.mat',prefix);
load(matfile);
matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);
expmatfile = sprintf('matfiles/%s_fiterr.mat',prefix);
load(expmatfile);
% contains fit_data and err


for c=1:length(conditions)
  k=conditions(c);
  fprintf('%s: preparing to plot traces for condition %d...\n',mfilename,k);

  fprintf('%s: creating temporary structure combining data, exp, and err...\n',...
    mfilename);
  tmp_avg = avg_data;
  tmp_avg.averages=[];
  tmp_avg.averages=[avg_data.averages(k),...
                    fit_data.averages(k),...
                    err.averages(k)];

  figure;
  fprintf('%s: plotting %s traces...\n',mfilename,chantype);
  ts_multiplot_avg(tmp_avg,'conditions',[1:3],'chantype',chantype,...
    'scale_max',scale_max);
  plot_ax = gca;
  ax=axes('Units','Normal','Position',[0.075 0.075 0.89 0.89],'Visible','off');
  set(get(ax,'Title'),'Visible','on');
  title(sprintf('Actual, Expected, and Error %s plots cond %d',chantype,k));
  set(gcf,'position',[50 50 imgsize imgsize]);
  axes(plot_ax);
  print('-dtiff',sprintf('traces-%s-%s-cond%d.tif',...
    prefix,chantype,k));

end;

%close all;
