function ts_plot_wforms(wforms,varargin)
%function ts_plot_wforms(wforms,[options])
%
% Purpose: plot waveforms
%
% Required Input:
%   wforms: 2D waveform matrix
%     size of wforms should be [ntpoints,nconds]
%
% Optional Parameters:
%   'outdir': output directory for plot
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'wforms'}
%   'wforms_err': waveform error
%     if size is [ntpoints,nconds], error is symmetric, relative to wforms
%     if size is [ntpoints,nconds,2], error is an absolute confidence interval
%     {default = []}
%   'ranges': struct array with info for drawing shaded time ranges
%     struct array must have nconds elements, and struct must have
%     fields 'onset' and 'offset' (each can be vector for multiple ranges)
%     onset and offset values should be in msec
%     {default = []}
%   'ranges_color': color string (e.g. 'k','r','g','b', etc.)
%     for drawing time ranges; if empty, will use 'colors'
%     {default = []}
%   'points': struct array with info for marking points (e.g. peaks)
%     struct array must have nconds elements, and struct must have
%      fields 'amplitude' and 'latency' (each can be vector)
%      latency values should be in msec
%     may also have field 'type' with values -1 or 1
%      will draw an X or O to mark points with negative or positive amplitudes
%     may also have field 'color' with color string
%     {default = []}
%   'points_color': color string (e.g. 'k','r','g','b', etc.)
%     for drawing points; if empty, will use 'colors'
%     {default = []}
%   'zero_line_flag': [0|1] draw black line at y = 0
%     {default = 0}
%   'xzero_line_flag': [0|1] draw black line at x = 0
%     {default = 0}
%   'time': time vector (msec)
%     if supplied, sfreq, t0, and t1 are ignored
%     length of time vector must match length of input wform
%     {default = []}
%   'sfreq': sampling frequency; used to calculate time vector
%     {default = 1000}
%   't0': start time of waveform (msec)
%     {default = -100}
%   't1': end time of waveform (msec)
%     {default = 300}
%   'xlim': x-axis limits for waveform plot
%     if empty, use entire time range
%     {default = []}
%   'ylim': y-axis limits for waveform plot
%     if empty, use auto-scaling
%     {default = []}
%   'relative_err_flag': [0|1] 0: absolute 1: relative
%     {default = 0}
%   'fill_err_flag': [0|1] whether to fill wforms_err as shaded region
%     Otherwise use error bars
%     {default = 1}
%   'fill_alpha': transparency of shaded region for wforms_err
%     Only applies if fill_err_flag = 1
%     {default = 0.1}
%   'errbar_interval': number of time samples between error bars
%     Only applies if fill_err_flag = 0
%     {default = 15}
%   'smooth_sigma': Gaussian blurring kernel applied to waveforms
%     {default = 0}
%   'offset': subtract this value from waveforms
%     {default = 0}
%   'baseline_flag': [0|1|2] subtract mean of baseline from waveforms
%     0: no baseline subtraction
%     1: subtract baseline of each condition
%     2: subtract baseline averaged across all conditions
%     {default = 0}
%   'baseline_t0': start time of baseline period
%     {default = -Inf}
%   'baseline_t1': start time of baseline period
%     {default = 0}
%   'colors': cell array of color codes for each condition
%     {default: {'b','g','r','y','c','m','k'}}
%   'linewidth': plot line width 
%     may be vector with nconds elements
%     {default = 1}
%   'axes_linewidth': axes line width
%     {default = 1}
%   'condnames': cell array of names of each condition for legend
%     if not supplied, no legend displayed
%     {default = []}
%   'legend_loc': location of area name legend
%     {default = 'SouthEast'}
%   'label_flag': [0|1] display labels on plot
%     {default = 1}
%   'title': plot title
%     {default = []}
%   'xlabel': x-axis label
%     {default = 'Time (msec)'}
%   'ylabel': y-axis label
%     {default = []}
%   'fontname': font name for labels and axes
%     {default = 'Arial'}
%   'fontsize': font size for labels and axes
%     {default = 12}
%   'tif_flag': save plot in tif format (bitmap)
%     {default = 1}
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
% Created:  01/24/12 by Don Hagler
% Last Mod: 01/24/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'outdir',pwd,[],...
  'outstem','wforms',[],...
  'wforms_err',[],[],...
  'ranges',[],[],...
  'ranges_color',[],[],...
  'points',[],[],...
  'points_color',[],[],...
  'zero_line_flag',false,[false,true],...
  'xzero_line_flag',false,[false,true],...
  'time',[],[],...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',300,[],...
  'xlim',[],[],...
  'ylim',[],[],...
  'relative_err_flag',false,[false true],...
  'fill_err_flag',true,[false true],...
  'fill_alpha',0.1,[0 1],...
  'errbar_interval',15,[],...
  'smooth_sigma',0,[0,1000],...
  'offset',0,[-Inf,Inf],...
  'baseline_flag',0,[0 1 2],...
  'baseline_t0',-Inf,[],...
  'baseline_t1',0,[],...
  'colors',{'b','g','r','y','c','m','k'},[],...
  'linewidth',1,[],...
  'axes_linewidth',1,[],...
  'condnames',[],[],...
  'legend_loc','SouthEast',[],...
  'label_flag',true,[false true],...
  'title',[],[],...
  'xlabel','Time (msec)',[],...
  'ylabel',[],[],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
  'tif_flag',true,[false true],...
  'eps_flag',false,[false true],...
  'tif_dpi',300,[10,10000],...
  'fig_size',[],[],...
  'visible_flag',true,[false true],...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(wforms,varargin,parms_filter);

if ~parms.tif_flag && ~parms.eps_flag && ~parms.visible_flag
  return;
end;

if parms.smooth_sigma
  wforms = smooth_waveforms(wforms,parms);
end;
if parms.baseline_flag || parms.offset
  [wforms,parms] = baseline_waveforms(wforms,parms);
end;
plot_waveforms(wforms,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(wforms,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  parms.ntpoints = size(wforms,1);
  parms.nconds = size(wforms,2);
  if ~isempty(parms.wforms_err)
    if size(parms.wforms_err,1)~=parms.ntpoints
      error('number of time points in wforms_err does not match wforms');
    end;
    if size(parms.wforms_err,2)~=parms.nconds
      error('number of conditions in wforms_err does not match wforms');
    end;
    if size(parms.wforms_err,3)==1
      parms.wforms_err = cat(3,wforms-parms.wforms_err,wforms+parms.wforms_err);
    elseif size(parms.wforms_err,3)==2 & parms.relative_err_flag
      parms.wforms_err = cat(3,wforms+parms.wforms_err(:,:,1),wforms+parms.wforms_err(:,:,2));
    elseif size(parms.wforms_err,3)>2
      error('number of wforms_err elements must be 1 or 2');
    end;
  end;
  parms.fontargs = {'FontSize',parms.fontsize,'FontName',parms.fontname};
  mmil_mkdir(parms.outdir);
  if mmil_isrelative(parms.outstem)
    parms.outstem = [parms.outdir '/' parms.outstem];
  end;
  if strcmp(parms.ylim,'none'), parms.ylim = []; end;
  if strcmp(parms.xlim,'none'), parms.xlim = []; end;
  if isempty(parms.time)
    sd = 1000/parms.sfreq;
    parms.time = [parms.t0:sd:parms.t1];
  end;
  if length(parms.time) ~= parms.ntpoints
    error('number of elements in time (%d) does not match wform (%d)',...
      length(parms.time),parms.ntpoints);
  end;
  parms.time = mmil_rowvec(parms.time);
  if parms.baseline_flag
    [tmp,parms.baseline_s0] = min(abs(parms.time - parms.baseline_t0));
    [tmp,parms.baseline_s1] = min(abs(parms.time - parms.baseline_t1));
  end;
  if isempty(parms.xlim)
    parms.xlim = [parms.time(1),parms.time(end)];
  end;
  if parms.nconds>1
    if length(parms.linewidth)==1
      parms.linewidth = parms.linewidth*ones(1,parms.nconds);
    elseif length(parms.linewidth)~=parms.nconds
      error('number of elements in linewidth (%d) does not match nconds (%d)',...
        length(parms.linewidth),parms.nconds);
    end;
    if length(parms.colors) < parms.nconds
      colors = parms.colors;
      k = 1;
      for c=1:parms.nconds
        if k>length(colors), k=1; end;
        parms.colors{c} = colors{k};
        k = k + 1;
      end;    
    end;
  end;
  if ~isempty(parms.condnames)
    if ~iscell(parms.condnames), parms.condnames = {parms.condnames}; end;
    if length(parms.condnames) ~= parms.nconds
      error('number of condnames (%d) does not match nconds (%d)',...
        length(parms.condnames),parms.nconds);
    end;
  end;
  % check that ranges has one element for each condition
  if ~isempty(parms.ranges)
    if ~isstruct(parms.ranges)
      error('ranges must be a struct or struct array');
    end;
    if length(parms.ranges) ~= parms.nconds
      error('ranges must have number of elements equal to nconds (%d)',...
        parms.nconds);    
    end;
    if ~isfield(parms.ranges,'onset') || ~isfield(parms.ranges,'offset')
      error('ranges must have fields ''onset'' and ''offset''');
    end;
    for c=1:parms.nconds
      onset = mmil_getfield(parms.ranges(c),'onset');
      offset = mmil_getfield(parms.ranges(c),'offset');
      if length(onset) ~= length(offset)
        error('ranges onset and offset vectors must have same length')
      end;      
      parms.ranges(c).nranges = length(onset);
      if ~isempty(parms.ranges_color)
        parms.ranges(c).color = parms.ranges_color;
      else
        parms.ranges(c).color = parms.colors{c};
      end;
    end;
  end;
  % check that points has one element for each condition
  if ~isempty(parms.points)
    if ~isstruct(parms.points)
      error('peaks must be a struct or struct array');
    end;
    if length(parms.points) ~= parms.nconds
      error('peaks must have number of elements equal to nconds (%d)',...
        parms.nconds);    
    end;
    if ~isfield(parms.points,'amplitude') || ~isfield(parms.points,'latency')
      error('points must have fields ''amplitude'' and ''latency''');
    end;
    for c=1:parms.nconds
      amplitude = mmil_getfield(parms.points(c),'amplitude');
      latency = mmil_getfield(parms.points(c),'latency');
      if length(amplitude) ~= length(latency)
        error('points amplitude and latency vectors must have same length')
      end;      
      if ~isfield(parms.points,'type') || isempty(parms.points(c).type)
        parms.points(c).type = sign(parms.points(c).amplitude);
      else
        if length(mmil_getfield(parms.points(c),'type'))~=length(amplitude)
          error('points type and amplitude vectors must have same length')
        end;
      end;
      parms.points(c).npoints = length(amplitude);
      if ~isempty(parms.points_color)
        parms.points(c).color = parms.points_color;
      elseif ~isfield(parms.points,'color') || isempty(parms.points(c).color)
        parms.points(c).color = parms.colors{c};
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wforms = smooth_waveforms(wforms,parms)
  for c=1:parms.nconds
    wforms(:,c) = mmil_smooth(squeeze(wforms(:,c)),parms.smooth_sigma);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wforms,parms] = baseline_waveforms(wforms,parms)
  if parms.baseline_flag == 2
    tmp_baseline = wforms(parms.baseline_s0:parms.baseline_s1,:);
    tmp_baseline = mean(tmp_baseline(:));
  end;
  for c=1:parms.nconds
    tmp_wform = wforms(:,c);
    if parms.baseline_flag == 1
      tmp_baseline = mean(tmp_wform(parms.baseline_s0:parms.baseline_s1));
    elseif parms.baseline_flag == 0;
      tmp_baseline = 0;
    end;
    if parms.offset
      tmp_baseline = tmp_baseline + parms.offset;
    end;
    % adjust wforms
    wforms(:,c) = tmp_wform - tmp_baseline;
    % adjust wforms_err
    if ~isempty(parms.wforms_err)
      parms.wforms_err(:,c,1) = parms.wforms_err(:,c,1) - tmp_baseline;
      parms.wforms_err(:,c,2) = parms.wforms_err(:,c,2) - tmp_baseline;
    end;
    % adjust points amplitudes
    if ~isempty(parms.points)
      for p=1:parms.points(c).npoints
        parms.points(c).amplitude(p) = ...
          parms.points(c).amplitude(p) - tmp_baseline;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_waveforms(wforms,parms)
  %figure;
  clf; hold on;
  % plot waveform
  for c=1:parms.nconds
    plot(parms.time,wforms(:,c),...
      'Color',parms.colors{c},'LineWidth',parms.linewidth(c));
  end;
  % plot shaded confidence interval or error bars
  if ~isempty(parms.wforms_err)
    for c=1:parms.nconds
      if parms.fill_err_flag
        % shaded confidence interval
        tmp_time = [parms.time fliplr(parms.time)];
        tmp_lo  = mmil_rowvec(parms.wforms_err(:,c,1));
        tmp_hi  = mmil_rowvec(parms.wforms_err(:,c,2));
        tmp_int = [tmp_hi fliplr(tmp_lo)];
        h = fill(tmp_time,tmp_int,parms.colors{c});            
        set(h,'FaceAlpha',parms.fill_alpha,'EdgeAlpha',parms.fill_alpha);
      else
        % plot error bars
        tmp_time = parms.time(1:parms.errbar_interval:end);
        tmp_vals = wforms(1:parms.errbar_interval:end,c);
        tmp_lo  = parms.wforms_err(1:parms.errbar_interval:end,c,1);
        tmp_hi  = parms.wforms_err(1:parms.errbar_interval:end,c,2);
        tmp_lo = tmp_vals - tmp_lo;
        tmp_hi = tmp_hi - tmp_vals;
        errorbar(tmp_time,tmp_vals,tmp_lo,tmp_hi,...
          '-','Color',parms.colors{c},'LineWidth',parms.linewidth(c));
      end;
    end;
  end;
  if isempty(parms.ylim)
    parms.ylim = get(gca,'ylim');
  end;
  % mark ranges
  if ~isempty(parms.ranges)
    for c=1:parms.nconds
      for r=1:parms.ranges(c).nranges
        t0 = parms.ranges(c).onset(r);
        t1 = parms.ranges(c).offset(r);
        col = parms.ranges(c).color;
        % shaded time range
        tmp_time = [t0:t1];
        tmp_lo = parms.ylim(1)*ones(size(tmp_time));
        tmp_hi = parms.ylim(2)*ones(size(tmp_time));
        tmp_int = [tmp_hi fliplr(tmp_lo)];
        tmp_time = [tmp_time fliplr(tmp_time)];
        h = fill(tmp_time,tmp_int,col);
        set(h,'FaceAlpha',parms.fill_alpha,'EdgeAlpha',parms.fill_alpha);
      end;
    end;
  end;
  % mark points (e.g. peaks)
  if ~isempty(parms.points)
    for c=1:parms.nconds
      for p=1:parms.points(c).npoints
        amplitude = parms.points(c).amplitude(p);
        latency = parms.points(c).latency(p);
        type = parms.points(c).type(p);
        col = parms.points(c).color;
        if type>0
          mark = 'O';
        else
          mark = 'X';
        end;
        plot(latency,amplitude,mark,'Color',col,'LineWidth',parms.linewidth(c));
      end;
    end;
  end;
  % draw zero line
  if parms.zero_line_flag
    line([parms.time(1),parms.time(end)],[0 0],'color',[0 0 0],...
      'LineWidth',parms.axes_linewidth);
  end;
  % draw xzero line
  if parms.xzero_line_flag
    line([0,0],[parms.ylim(1),parms.ylim(2)],'color',[0 0 0],...
      'LineWidth',parms.axes_linewidth);
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
    plot_waveform_labels(parms,parms.title,parms.ylabel);
  end;
  if parms.tif_flag || parms.eps_flag
    save_plot(parms,parms.outstem);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_waveform_labels(parms,toplabl,ylabl)
  xlabel(parms.xlabel,parms.fontargs{:});
  if ~isempty(ylabl)
    ylabel(ylabl,parms.fontargs{:});
  end;
  if ~isempty(toplabl)
    toplabl = regexprep(toplabl,'_',' ');
    title(toplabl,parms.fontargs{:});
  end;
  if ~isempty(parms.condnames)
    legend(parms.condnames,parms.fontargs{:},'Location',parms.legend_loc);
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
