function ts_plot_responses(responses,varargin)
%function ts_plot_responses(responses,[options])
%
% Purpose: plot waveforms
%
% Required Input:
%   responses: matrix of values to be plotted (on y-axis)
%     e.g. amplitude, latency, auc, etc.
%     size should be [nareas,nconds]
%
% Optional Parameters:
%   'outdir': output directory for plot
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'responses'}
%   'responses_err': measurement error of responses (for error bars)
%     if size is [nareas,nconds], error is symmetric, relative to responses
%     if size is [nareas,nconds,2], error is an absolute confidence interval
%     {default = []}
%   'conditions': vector of condition values (x-axis)
%     e.g. contrast, spatial frequency, condition number, etc.
%     if not supplied, will be [1:nconds] with nconds=length(responses)
%     must have same length as responses
%     {default = []}
%   'xlim': x-axis limits for responses plot
%     if empty, use auto-scaling
%     {default = []}
%   'ylim': y-axis limits for responses plot
%     if empty, use auto-scaling
%     {default = []}
%   'relative_err_flag': [0|1] 0: absolute 1: relative
%     {default = 0}
%   'fill_err_flag': [0|1] whether to fill responses_err as shaded region
%     Otherwise use error bars
%     {default = 1}
%   'fill_alpha': transparency of shaded region for wforms_err
%     Only applies if fill_err_flag = 1
%     {default = 0.1}
%   'normflag': normalize responses
%     0: raw values
%     1: normalized to maximum condition for each area
%     2: normalized to max value for each area
%     {default = 0}
%   'colors': cell array of color codes for each area
%     {default: {'b','g','r','y','c','m','k'}}
%   'linewidth': plot line width 
%       may be vector with nareas elements
%     {default = 1}
%   'axes_linewidth': axes line width
%     {default = 1}
%   'markersize': plot marker size (points)
%     {default = 6}
%   'roinames': cell array of names of each area for legend
%     if not supplied, no legend displayed
%     {default = []}
%   'legend_loc': location of area name legend
%     {default = 'SouthEast'}
%   'label_flag': [0|1] display labels on plot
%     {default = 1}
%   'title': plot title
%     {default = []}
%   'xlabel': x-axis label
%     {default = []}
%   'ylabel': y-axis label
%     {default = []}
%   'logx_flag': [0|1] use log scale for x-axis
%     {default = 0}
%   'zero_line_flag': [0|1] draw black line at y = 0
%     {default = 0}
%   'fontname': font name for labels and axes
%     {default = 'Arial'}
%   'fontsize': font size for labels and axes
%     {default = 12}
%   'eps_flag': save plot in eps format (vector graphics)
%     {default = 0}
%   'tif_dpi': resolution of tif files (dots per inch)
%     {default = 300}
%   'fig_size': figure size in inches
%     if not specified, will use default
%     {default = []}
%   'visible_flag': [0|1] whether to display plot to screen
%     if 0, save only
%     {default = 1}
%
% Created:  02/09/12 by Don Hagler
% Last Mod: 03/15/14 by Don Hagler
%

%% todo: option to display colored * to mark significance (in column if multiple)

%% todo: allow multi-subject input; size(responses) = [nareas,nconds,nsubs]
%% todo: linear fit, power law, anonymized function?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
%% todo: check these, may not need all
parms_filter = {...
  'outdir',pwd,[],...
  'outstem','responses',[],...
  'responses_err',[],[],...
  'conditions',[],[],...
  'xlim',[],[],...
  'ylim',[],[],...
  'relative_err_flag',false,[false true],...
  'fill_err_flag',true,[false true],...
  'fill_alpha',0.1,[0 1],...
  'normflag',0,[0,1,2],...
  'colors',{'b','g','r','y','c','m','k'},[],...
  'linewidth',1,[],...
  'axes_linewidth',1,[],...
  'markersize',6,[],...
  'roinames',[],[],...
  'legend_loc','SouthEast',[],...
  'label_flag',true,[false true],...
  'title',[],[],...
  'xlabel',[],[],...
  'ylabel',[],[],...
  'logx_flag',false,[false true],...
  'zero_line_flag',false,[false,true],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
  'eps_flag',false,[false true],...
  'tif_dpi',300,[10,10000],...
  'fig_size',[],[],...
  'visible_flag',true,[false true],...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(responses,varargin,parms_filter);
plot_responses(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(responses,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  parms.nareas = size(responses,1);
  parms.nconds = size(responses,2);
  if isempty(parms.conditions)
    parms.conditions = [1:parms.nconds];
  end;
  [tmp,parms.ind_max_cond] = max(parms.conditions);
  if ~isempty(parms.responses_err)
    if size(parms.responses_err,1)~=parms.nareas
      error('number of areas in responses_err does not match responses');
    end;
    if size(parms.responses_err,2)~=parms.nconds
      error('number of conditions in responses_err does not match responses');
    end;
    if size(parms.responses_err,3)==1
      parms.responses_err = cat(3,responses-parms.responses_err,responses+parms.responses_err);
    elseif size(parms.responses_err,3)==2 & parms.relative_err_flag
      parms.responses_err = cat(3,responses+parms.responses_err(:,:,1),responses+parms.responses_err(:,:,2));
    elseif size(parms.responses_err,3)>2
      error('number of responses_err elements must be 1 or 2');
    end;
  end;
  parms.fontargs = {'FontSize',parms.fontsize,'FontName',parms.fontname};
  mmil_mkdir(parms.outdir);
  if mmil_isrelative(parms.outstem)
    parms.outstem = [parms.outdir '/' parms.outstem];
  end;
  if strcmp(parms.ylim,'none'), parms.ylim = []; end;
  if strcmp(parms.xlim,'none'), parms.xlim = []; end;
  if parms.nareas>1
    if length(parms.linewidth)==1
      parms.linewidth = parms.linewidth*ones(1,parms.nareas);
    elseif length(parms.linewidth)~=parms.nareas
      error('number of elements in linewidth (%d) does not match nareas (%d)',...
        length(parms.linewidth),parms.nareas);
    end;
    if length(parms.colors) < parms.nareas
      colors = parms.colors;
      k = 1;
      for c=1:parms.nareas
        if k>length(colors), k=1; end;
        parms.colors{c} = colors{k};
        k = k + 1;
      end;    
    end;
  end;
  if ~isempty(parms.roinames)
    if ~iscell(parms.roinames), parms.roinames = {parms.roinames}; end;
    if length(parms.roinames) ~= parms.nareas
      error('number of roinames (%d) does not match nareas (%d)',...
        length(parms.roinames),parms.nareas);
    end;
  end;
  parms.responses = responses;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_responses(parms)
  clf; hold on;
  % plot responses
  for a=1:parms.nareas
    resps = squeeze(parms.responses(a,:));
    if parms.normflag
      if parms.normflag==1
        normval = resps(parms.ind_max_cond);
      elseif parms.normflag==2
        [~,ind] = max(abs(resps));
        normval = resps(ind);
      end;
      resps = resps/(eps+normval);
    end;
    % plot amplitude vs. conditions
    if parms.linewidth(a)>0
      format_str = 'o-';
      tmp_args = {'Color',parms.colors{a},'LineWidth',parms.linewidth(a),...
        'MarkerSize',parms.markersize};
    else
      format_str = 'o';
      tmp_args = {'Color',parms.colors{a},'MarkerSize',parms.markersize};
    end;
    plot(parms.conditions,resps,format_str,tmp_args{:});
  end;
  % plot err
  if ~isempty(parms.responses_err)
    for a=1:parms.nareas
      resps = squeeze(parms.responses(a,:));
      lo = squeeze(parms.responses_err(a,:,1));
      hi = squeeze(parms.responses_err(a,:,2));
      if ~parms.fill_err_flag
        lo = resps - lo;
        hi = hi - resps;
      end;
      if parms.normflag
        if parms.normflag==1
          normval = resps(parms.ind_max_cond);
        elseif parms.normflag==2
          [~,ind] = max(abs(resps));
          normval = resps(ind);
        end;
        resps = resps/(eps+normval);
        lo = lo/(eps+normval);
        hi = hi/(eps+normval);
      end;
      % plot amplitude vs. conditions
      format_str = 'o';
      tmp_args = {'Color',parms.colors{a},'MarkerSize',parms.markersize};
      if ~parms.fill_err_flag
        errorbar(parms.conditions,resps,lo,hi,format_str,tmp_args{:});
      else
        tmp_conds = mmil_rowvec(parms.conditions);
        tmp_conds = [tmp_conds fliplr(tmp_conds)];
        tmp_int = [mmil_rowvec(hi) fliplr(mmil_rowvec(lo))];
        h = fill(tmp_conds,tmp_int,parms.colors{a});            
        set(h,'FaceAlpha',parms.fill_alpha,'EdgeAlpha',parms.fill_alpha);
      end;
    end;
  end;
  % draw zero line
  if parms.zero_line_flag
    line([parms.xlim(1),parms.xlim(2)],[0 0],'color',[0 0 0]);
  end;
  % set font for axes
  set(gca,parms.fontargs{:});
  % set linewidth for axes
  set(gca,'LineWidth',parms.axes_linewidth);
  % set limits for axes
  if ~isempty(parms.ylim), set(gca,'YLim',parms.ylim); end;
  if ~isempty(parms.xlim), set(gca,'XLim',parms.xlim); end;
  % add labels
  if parms.label_flag
    plot_response_labels(parms)
  end;
  save_plot(parms,parms.outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_response_labels(parms)
  xlabel(parms.xlabel,parms.fontargs{:});
  if ~isempty(parms.ylabel)
    ylabel(parms.ylabel,parms.fontargs{:});
  end;
  if ~isempty(parms.title)
    title(parms.title,parms.fontargs{:});
  end;
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
  if ~isempty(parms.roinames)
    legend(parms.roinames,parms.fontargs{:},'Location',parms.legend_loc);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_plot(parms,outstem)
  if ~parms.visible_flag, set(gcf,'Visible','off'); end;
  if ~isempty(parms.fig_size)
    if length(parms.fig_size)==1
      parms.fig_size = [parms.fig_size parms.fig_size];
    end;
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [parms.fig_size(1) parms.fig_size(2)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 parms.fig_size(1) parms.fig_size(2)]);
  end;
  print(gcf,'-dtiff',[outstem '.tif'],sprintf('-r %d',parms.tif_dpi));
  if parms.eps_flag, mmil_printeps(gcf,[outstem '.eps']); end;
  if ~parms.visible_flag, close(gcf); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
