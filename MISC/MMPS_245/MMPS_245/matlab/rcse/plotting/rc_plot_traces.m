function rc_plot_traces(prefix,conditions,chantype,scale_max,linewidth)
%function rc_plot_traces(prefix,[conditions],[chantype],[scale_max],[linewidth])
%
% Purpose: multiplot MEG/EEG sensor waveforms for multiple conditions
%
% Required Input:
%   prefix: proc or RCSE prefix
%
% Optional Input:
%   conditions: vector of conditions numbers to overlay
%     {default = 1}
%   chantype: e.g. 'grad1', 'grad2', 'mag', 'eeg'
%     {default = 'grad1'}
%   scale_max: y-axis max scale value
%     if 0, auto-scale
%     {default = 0}
%   linewidth
%     {default = 1.5}
%
% Early Mod: 12/22/06 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%% todo: varargin
%% todo: rootdir
%% todo: merge with rc_plot_traces_datafit?

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('conditions','var') || isempty(conditions), conditions = 1; end;
if ~exist('chantype','var') || isempty(chantype), chantype = 'grad1'; end;
if ~exist('scale_max','var') || isempty(scale_max), scale_max = 0; end;
if ~exist('linewidth','var') || isempty(linewidth), linewidth = 1.5; end;
imgsize = 1000;

matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
load(matfile);

fprintf('%s: preparing to plot traces...\n',mfilename);

figure;
fprintf('%s: plotting %s traces...\n',mfilename,chantype);
ts_multiplot_avg(avg_data,'conditions',conditions,'chantype',chantype,...
  'scale_max',scale_max,'linewidth',linewidth);
plot_ax = gca;
print('-dtiff',sprintf('traces-%s-%s.tif',...
  prefix,chantype));

%close all;
