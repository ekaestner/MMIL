function cfg = singleplotTFR(cfg, data)

% singleplotTFR plots the time-frequency representations of power of a
% single channel or the average over multiple channels.
%
% Use as:
%   singleplotTFR(cfg,data)
%
% The data can be a time-frequency representation of power that was
% computed using the FREQANALYSIS function. 
%
% The configuration can have the following parameters:
% cfg.xparam        = field to be plotted on x-axis, e.g. 'time' (default depends on data.dimord)
% cfg.yparam        = field to be plotted on y-axis, e.g. 'freq' (default depends on data.dimord)
% cfg.zparam        = field to be plotted on y-axis, e.g. 'powspcrtrm' (default depends on data.dimord)
% cfg.maskparameter = field in the data to be used for opacity masking of data 
%                     (not possible for mean over multiple channels)
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.zlim          = 'maxmin','absmax' or [zmin zmax] (default = 'maxmin')
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                     see CHANNELSELECTION for details
% cfg.fontsize      = font size of title (default = 8)
% cfg.colormap      = any sized colormap, see COLORMAP
% cfg.colorbar      = 'yes', 'no' (default = 'yes')
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas 
%                     can be selected by holding down the SHIFT key.
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = 'opengl')
% cfg.masknans      = 'yes' or 'no' (default = 'yes')
%
% See also:
%   singleplotER, multiplotER, multiplotTFR, topoplotER, topoplotTFR.

% Undocumented local options:
% cfg.channelname
% cfg.channelindex
%
% This function depends on FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype, documented

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: singleplotTFR.m,v $
% Revision 1.24  2008/01/29 19:43:33  sashae
% added option for trial selection; plot functions now also accept data with
% repetitions (either trials or subjects), the avg is computed and plotted
% removed some old code
%
% Revision 1.23  2007/06/19 15:57:51  ingnie
% fixed bug in cfg.maskparameter, thanks to Saskia
%
% Revision 1.22  2007/06/19 14:01:54  ingnie
% added cfg.maskparameter, changed some white spaces
%
% Revision 1.21  2007/06/14 12:23:48  ingnie
% added cfg.colormap option
%
% Revision 1.20  2007/04/26 09:58:46  ingnie
% default masknans to 'yes'
%
% Revision 1.19  2007/04/25 17:25:06  ingnie
% added cfg.masknans option
%
% Revision 1.18  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.17  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.16  2007/01/09 10:41:43  roboos
% Added the option cfg.renderer, default is opengl. Interactive plotting
% on linux/VNC sometimes does not work, using cfg.renderer='painters'
% seems to fix it.
%
% Revision 1.15  2006/10/24 12:09:54  ingnie
% added colorbar yes/no option
%
% Revision 1.14  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.13  2006/06/07 10:24:34  ingnie
% added error when no channels are selected
%
% Revision 1.12  2006/06/06 16:24:06  ingnie
% replaced cfg.channelindex and cfg.channelname with cfg.channel for consistency
% with other functions
%
% Revision 1.11  2006/05/30 14:18:02  ingnie
% fixed bug that appeared when plotting single channel, updated documentation
%
% Revision 1.10  2006/05/09 17:33:57  ingnie
% fixed bug that appeared when cfg.channelindex is more than one channel
%
% Revision 1.9  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.8  2006/04/20 09:57:53  roboos
% changed formatting of the code from DOS into UNIX, i.e. removed the <CR>
%
% Revision 1.7  2006/04/07 23:45:23  chrhes
% changed clf to cla at beginning of function to allow use in subplots
%
% Revision 1.6  2006/03/22 18:56:57  jansch
% removed hard-coded selection of the powerspectrum in the to be plotted data
%
% Revision 1.5  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
% 

cla

% For backward compatibility with old data structures:
data = checkdata(data);

% Set the defaults:
if ~isfield(cfg,'baseline'),        cfg.baseline = 'no';               end
if ~isfield(cfg,'baselinetype'),    cfg.baselinetype = 'absolute';     end
if ~isfield(cfg,'trials'),          cfg.trials = 'all';                end
if ~isfield(cfg,'xlim'),            cfg.xlim = 'maxmin';               end
if ~isfield(cfg,'ylim'),            cfg.ylim = 'maxmin';               end
if ~isfield(cfg,'zlim'),            cfg.zlim = 'maxmin';               end
if ~isfield(cfg,'fontsize'),        cfg.fontsize = 8;                  end
if ~isfield(cfg,'colorbar'),        cfg.colorbar = 'yes';              end
if ~isfield(cfg,'interactive'),     cfg.interactive = 'no';            end
if ~isfield(cfg,'renderer'),        cfg.renderer = 'opengl';           end
if ~isfield(cfg,'masknans'),        cfg.masknans = 'yes';              end
if ~isfield(cfg,'maskparameter'),   cfg.maskparameter = [];            end

% Set x/y/zparam defaults according to data.dimord value:
if strcmp(data.dimord, 'chan_freq_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
elseif strcmp(data.dimord, 'subj_chan_freq_time') || strcmp(data.dimord, 'rpt_chan_freq_time')
  if isfield(data, 'crsspctrm'),   data = rmfield(data, 'crsspctrm');  end % on the fly computation of coherence spectrum is not supported
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  tmpcfg.jackknife = 'no';
  data = freqdescriptives(tmpcfg, data);
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
end

%Pick the channel(s)
if ~isfield(cfg,'channel') % for backward compatibility
  if isfield(cfg,'channelindex') && isfield(cfg,'channelname') 
    warning('cfg.channelindex is old, please use cfg.channel instead');
    warning('cfg.channelname is old, please use cfg.channel instead');
    warning('cfg.channelname is used');
    cfg.channel = cfg.channelname;
    cfg = rmfield(cfg,'channelindex');
    cfg = rmfield(cfg,'channelname');
  elseif isfield(cfg,'channelindex') 
    warning('cfg.channelindex is old, please use cfg.channel instead')
    cfg.channel = cfg.channelindex;
    cfg = rmfield(cfg,'channelindex');
  elseif isfield(cfg,'channelname')
    warning('cfg.channelname is old, please use cfg.channel instead')
    cfg.channel = cfg.channelname;
    cfg = rmfield(cfg,'channelname');
  else
    % set the default
    cfg.channel = 'all';
  end
end

cfg.channel = channelselection(cfg.channel, data.label);
if isempty(cfg.channel)
  error('no channels selected');
else
  chansel = match_str(data.label, cfg.channel);
end

% cfg.maskparameter only possible for single channel
if length(chansel) > 1 && ~isempty(cfg.maskparameter)
  warning('no masking possible for average over multiple channels')
  masking = 0;
elseif length(chansel) == 1 && ~isempty(cfg.maskparameter)
  masking = 1;
end

% Apply baseline correction:
if ~strcmp(cfg.baseline, 'no')
  data = freqbaseline(cfg, data);
end

% Get physical x-axis range:
if strcmp(cfg.xlim,'maxmin')
  xmin = min(data.(cfg.xparam));
  xmax = max(data.(cfg.xparam));
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Find corresponding x-axis bins:
xidc = find(data.(cfg.xparam) >= xmin & data.(cfg.xparam) <= xmax);

% Align physical x-axis range to the array bins:
xmin = data.(cfg.xparam)(xidc(1));
xmax = data.(cfg.xparam)(xidc(end));

% Get physical y-axis range:
if strcmp(cfg.ylim,'maxmin')
  ymin = min(data.(cfg.yparam));
  ymax = max(data.(cfg.yparam));
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Find corresponding y-axis bins:
yidc = find(data.(cfg.yparam) >= ymin & data.(cfg.yparam) <= ymax);

% Align physical y-axis range to the array bins:
ymin = data.(cfg.yparam)(yidc(1));
ymax = data.(cfg.yparam)(yidc(end));

% Get TFR data averaged across selected channels, within the selected x/y-range:
dat = getsubfield(data, cfg.zparam);
TFR = squeeze(mean(dat(chansel,yidc,xidc), 1));
if ~isempty(cfg.maskparameter)
  mas = getsubfield(data, cfg.maskparameter);
  mdata = squeeze(mas(chansel,yidc,xidc));
end

% Get physical z-axis range (color axis):
if strcmp(cfg.zlim,'maxmin')
  zmin = min(TFR(:));
  zmax = max(TFR(:));
elseif strcmp(cfg.zlim,'absmax')
  zmin = -max(abs(TFR(:)));
  zmax = max(abs(TFR(:)));
else
  zmin = cfg.zlim(1);
  zmax = cfg.zlim(2);
end

% Draw plot:
hold on;
h = imagesc(data.time(xidc), data.freq(yidc), TFR, [zmin,zmax]);
% Mask Nan's and maskfield
if isequal(cfg.masknans,'yes') && isempty(cfg.maskparameter) || ~masking
  mask = ~isnan(TFR);
  mask = double(mask);
  set(h,'AlphaData',mask, 'AlphaDataMapping', 'scaled');
  alim([0 1]);
elseif isequal(cfg.masknans,'yes') && ~isempty(cfg.maskparameter) && masking
  mask = ~isnan(TFR);
  mask = mask .* mdata;
  mask = double(mask);
  set(h,'AlphaData',mask, 'AlphaDataMapping', 'scaled');
  alim([0 1]);
elseif isequal(cfg.masknans,'no') && ~isempty(cfg.maskparameter) && masking
  mask = mdata;
  mask = double(mask);
  set(h,'AlphaData',mask, 'AlphaDataMapping', 'scaled');
  alim([0 1]);
end
axis xy;

%% added by JSS on 20-Mar-2009
% add tick marks & labels if timesurfer standard foi are used
  mmilfoi = [2:12 14:2:24 25:5:55 70:10:200];
  plotfoi = data.(cfg.yparam)(yidc);
  if all(ismember(plotfoi,mmilfoi))
    xas = data.time(xidc);
    yas = data.freq(yidc);
    cnr = axis;
    yval          = [min(plotfoi) 12 25 55 200];
    [yval,jnk,yi] = intersect(yval,plotfoi);
%     ytick         = yas(sort(yi));
    ytick = sort((yi/max(yi))*max(yas));
    dy            = yi(2)-yi(1);
    xn = length(xas); dx = floor((xn-1) / 3); xh = xas(dx);  
    xi = [1:dx:xn];  
    xtick = xas(xi); xval = data.(cfg.xparam)(xidc(xi)); % (xt,xval) x-axis    
    for xx = 1:length(xval)
      xvals{xx} = sprintf('%g',round(10*xval(xx))/10);
%       ticklength = .02*(yas(dy)-yas(1));
      xticklength = (xas(end)-xas(1))/40;
      yticklength = (yas(end)-yas(1))/40;
      line([xtick(xx) xtick(xx)],[yas(1)-yticklength yas(1)],'color','k','linewidth',1);% yas(1)+.3*[-(yas(dy)-yas(1)) 0],'color','k','linewidth',.7)  
%       text(xtick(xx),yas(1)-(yas(dy)-yas(1)),xvals{xx},'fontsize',7);
      text(xtick(xx),yas(1)-2*yticklength,xvals{xx},'fontsize',10);
    end
    for yy = 1:length(yval)
      yvals{yy} = sprintf('%g',yval(yy));
      line([xas(1) xas(1)+xticklength],[ytick(yy) ytick(yy)],'color','k','linewidth',1)
%       text(xas(1)-(xas(dx)-xas(1)),ytick(yy),yvals{yy},'fontsize',7);
      text(xas(1)-2.5*xticklength,ytick(yy),yvals{yy},'fontsize',10);
    end
    axis off
  end
%%
% set colormap
if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('singleplotTFR(): Colormap must be a n x 3 matrix'); end
  colormap(cfg.colormap);
end;

if isequal(cfg.colorbar,'yes')
colorbar;
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  userData.hFigure = gcf;
  userData.hAxes = gca;
  for i=1:1 % no multiple selection regions
    userData.hSelection{i} = plot(mean([xmin xmax]),mean([ymin ymax]));
    set(userData.hSelection{i}, 'XData', mean([xmin xmax]));
    set(userData.hSelection{i}, 'YData', mean([ymin ymax]));
    set(userData.hSelection{i}, 'Color', [0 0 0]);
    set(userData.hSelection{i}, 'EraseMode', 'xor');
    set(userData.hSelection{i}, 'LineStyle', '--');
    set(userData.hSelection{i}, 'LineWidth', 1.5);
    set(userData.hSelection{i}, 'Visible', 'on');
    userData.range{i} = [];
  end
  userData.iSelection = 0;
  userData.plotType = 'singleplot';
  userData.selecting = 0;
  userData.selectionType = '';
  userData.selectAxes = 'xy';
  userData.lastClick = [];
  userData.cfg = cfg;
  userData.data = data;
  userData.chanX = [];
  userData.chanY = [];
  userData.chanLabels = [];
  tag = sprintf('%.5f', 10000 * rand(1));
  set(gcf, 'Renderer', cfg.renderer);
  set(gcf, 'Tag', tag);
  set(gcf, 'UserData', userData);
  set(gcf, 'WindowButtonMotionFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 0);']);
  set(gcf, 'WindowButtonDownFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 1);']);
  set(gcf, 'WindowButtonUpFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 2);']);
end

if length(chansel) == 1
  str = [char(cfg.channel) ' / ' num2str(chansel)];
else
  str = sprintf('mean(%0s)', join(',',cfg.channel));
end
if isfield(cfg,'title')
  title(sprintf('%s (%s)',cfg.title,str));
else
  title(str);
end

axis tight;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l = cellstrmatch(str,strlist)
l = [];
for k=1:length(strlist)
  if strcmp(char(str),char(strlist(k)))
    l = [l k];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = join(separator,cells)
if length(cells)==0
  t = '';
  return;
end
t = char(cells{1});

for i=2:length(cells)
  t = [t separator char(cells{i})];
end
