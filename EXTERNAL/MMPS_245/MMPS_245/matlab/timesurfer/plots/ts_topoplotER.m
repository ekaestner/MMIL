function ts_topoplotER(cfg, varargin)
%function ts_topoplotER(cfg, varargin)
%
% topoplotER plots the topographic distribution of event-related fields or 
% potentials or of the scalp distribution of oscillatory activity (power or
% coherence).
%
% Use as:
%   ts_topoplotER(cfg, data)
%
% Optional parameters & values:
%
% cfg.xparam        = field to be plotted on x-axis
% cfg.zparam        = field to be plotted as color
% cfg.layout        = specification of the layout, see below
% cfg.xlim          = 'maxmin' or [xmin xmax]         (default = 'maxmin')
% cfg.zlim          = 'maxmin' or [zmin zmax]         (default = 'maxmin')
% cfg.baseline      = 'yes','no' or [time1 time2]     (default = 'no')
% cfg.colormap      = any sized colormap
% cfg.style         = topoplot style                  (default = 'both')
%                     'straight' colormap only
%                     'contour' contour lines only
%                     'both' both colormap and contour lines
%                     'fill' constant color between lines
%                     'blank' just head and electrodes
% cdg.gridscale     = scaling grid size               (default = 67) 
% cfg.interplimits  = limits for interpolation        (default = 'head')
%                     'electrodes' to furthest electrode
%                     'head' to edge of head
% cfg.interpolation = interpolation method            (default = 'v4')
%                     'v4' MATLAB 4 griddata method
%                     'linear' triangle-based linear interpolation 
%                     'cubic' triangle-based cubic interpolation
%                     'nearest' nearest neighbor interpolation
% cfg.contournum    = number of contour lines         (default = 6)
% cfg.shading       = 'flat' or 'interp'              (default = 'flat')
% cfg.comment       = string of text                  (default = date)
%                     Add 'comment' to graph (according to COMNT in the layout)
% cfg.ecolor        = Marker color                    (default = [0 0 0] (black))
% cfg.emarker       = Marker symbol                   (default = 'o')
% cfg.emarkersize   = Marker size                     (default = 2)
% cfg.fontsize      = Font size of labels/comment     (default = 8 pt)
% cfg.hcolor        = Head lines color                (default = [0 0 0] (black))
% cfg.hlinewidth    = Head lines width                (default = 2)
% cfg.interactive   = Interactive plot 'yes', 'no'    (default = 'no')
%
% If you use the output of freqdescriptives as an input and if you want to visualise
% coherence with respect to a reference-channel, you should specify:
%
% cfg.zparam        =  'cohspctrm'
% cfg.cohrefchannel =  name of reference-channel
%
% The layout defines how the channels are arranged. You can specify the
% layout in a variety of ways:
%  - you can give the name of an ascii layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition, i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these and the data structure contains an
% electrode or gradiometer structure, that will be used for creating a
% layout. If you want to have more fine-grained control over the layout
% of the subplots, you should create your own layout file.
%
% See also:
%   topoplotTFR, singleplotER, singleplotTFR, multiplotER, multiplotTFR.
%
% The following additional cfg parameters are used when plotting 3-dimensional
% data (i.e. when topoplotTFR calls topoplotER):
% cfg.yparam          field to be plotted on y-axis
% cfg.ylim            'maxmin' or [ymin ymax]         (default = 'maxmin')
%
% Created:  ? from fieldtrip
% Last Mod: 02/09/12 by Don Hagler
%

cla

% Set some constants:
rmax = .5;

% Multiple data sets are not supported for topoplot:
if length(varargin)>1
  error('Multiple data sets are not supported for topoplotER/topoplotTFR.');
end

data = varargin{1};

% For backward compatibility with old data structures:
data = fixdimord(data);

% Set other config defaults:
if ~isfield(cfg, 'xlim'),          cfg.xlim = 'maxmin';       end
if ~isfield(cfg, 'ylim'),          cfg.ylim = 'maxmin';       end
if ~isfield(cfg, 'zlim'),          cfg.zlim = 'maxmin';       end
if ~isfield(cfg, 'baseline'),      cfg.baseline = 'no';       end;
if ~isfield(cfg, 'style'),         cfg.style = 'both';        end
if ~isfield(cfg, 'gridscale'),     cfg.gridscale = 67;        end
if ~isfield(cfg, 'interplimits'),  cfg.interplimits = 'head'; end
if ~isfield(cfg, 'interpolation'), cfg.interpolation = 'v4';  end
if ~isfield(cfg, 'contournum'),    cfg.contournum = 6;        end
if ~isfield(cfg, 'shading'),       cfg.shading = 'flat';      end
if ~isfield(cfg, 'comment'),       cfg.comment = date;        end
if ~isfield(cfg, 'ecolor'),        cfg.ecolor = [0 0 0];      end
if ~isfield(cfg, 'emarker'),       cfg.emarker = 'o';         end
if ~isfield(cfg, 'emarkersize'),   cfg.emarkersize = 2;       end
if ~isfield(cfg, 'fontsize'),      cfg.fontsize = 8;          end
if ~isfield(cfg, 'hcolor'),        cfg.hcolor = [0 0 0];      end
if ~isfield(cfg, 'hlinewidth'),    cfg.hlinewidth = 2;        end
if ~isfield(cfg, 'interactive'),   cfg.interactive = 'no';    end

% Set x/y/zparam defaults according to data.dimord value:
if strcmp(data.dimord, 'chan_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';             end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='avg';          end
elseif strcmp(data.dimord, 'chan_freq')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='freq';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';             end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';    end
elseif strcmp(data.dimord, 'chan_freq_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='freq';         end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';    end
elseif strcmp(data.dimord, 'chan_comp')
  % Add a pseudo-axis with the component numbers:
  data.comp = 1:size(data.topo,2);
  % Rename the field with topographic label information:
  data.label = data.topolabel;
  if ~isfield(cfg, 'xparam'),      cfg.xparam='comp';         end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';             end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='topo';         end
end

% Old style coherence plotting with cohtargetchannel is no longer supported:
if isfield(cfg,'cohtargetchannel'), 
  error('cfg.cohtargetchannel is obsolete, check the documentation for help about coherence plotting.');
end

% Create time-series of small topoplots:
if ~ischar(cfg.xlim) & length(cfg.xlim)>2
  % Switch off interactive mode:
  cfg.interactive = 'no';
  xlims = cfg.xlim;
  % Iteratively call topoplotER with different xlim values:
  for i=1:length(xlims)-1
    subplot(ceil(sqrt(length(xlims)-1)), ceil(sqrt(length(xlims)-1)), i);
    cfg.xlim = xlims(i:i+1);
    ts_topoplotER(cfg, data);
  end
  return
end

% Check for unconverted coherence spectrum data:
if (strcmp(cfg.zparam,'cohspctrm')) && isfield(data, 'labelcmb')
  % A reference channel is required:
  if ~isfield(cfg,'cohrefchannel'),
    error('no reference channel specified');
  end
  % Convert 2-dimensional channel matrix to a single dimension:
  sel1           = strmatch(cfg.cohrefchannel, data.labelcmb(:,2));
  sel2           = strmatch(cfg.cohrefchannel, data.labelcmb(:,1));
  fprintf('selected %d channels for coherence\n', length(sel1)+length(sel2));
  data.cohspctrm = abs(data.cohspctrm([sel1;sel2],:,:));
  data.label     = [data.labelcmb(sel1,1);data.labelcmb(sel2,2)];
  data           = rmfield(data, 'labelcmb');
end

% Apply baseline correction:
if ~strcmp(cfg.baseline, 'no')
  if strcmp(cfg.xparam, 'time')
    data = timelockbaseline(cfg, data);
  elseif strcmp(cfg.xparam, 'freq')
    data = freqbaseline(cfg, data);
  end
end

% Get physical min/max range of x:
if strcmp(cfg.xlim,'maxmin')
  xmin = min(data.(cfg.xparam)(:));
  xmax = max(data.(cfg.xparam)(:));
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Replace value with the index of the nearest bin
xmin = nearest(data.(cfg.xparam)(:), xmin);
xmax = nearest(data.(cfg.xparam)(:), xmax);

if ~isempty(cfg.yparam)
  if strcmp(cfg.ylim,'maxmin')
    ymin = min(data.(cfg.yparam)(:));
    ymax = max(data.(cfg.yparam)(:));
  else
    ymin = cfg.ylim(1);
    ymax = cfg.ylim(2);
  end
  
  % Replace value with the index of the nearest bin:
  ymin = nearest(data.(cfg.yparam)(:), ymin);
  ymax = nearest(data.(cfg.yparam)(:), ymax);
end

dat = data.(cfg.zparam);
if ~isempty(cfg.yparam)
  dat = dat(:,ymin:ymax,xmin:xmax);
  dat = ts_nan_mean(ts_nan_mean(dat, 2), 3);
else
  dat = dat(:,xmin:xmax);
  dat = ts_nan_mean(dat, 2);
end
dat = dat(:);

% Read or create the layout that will be used for plotting
if ~isfield(cfg, 'layout')
  if isfield(cfg, 'layoutname')
    cfg.layout = cfg.layoutname;    % backward compatible cfg fieldname
  elseif isfield(data,'grad')
    cfg.layout = data.grad;         % create layout from gradiometer definition
  elseif isfield(data, 'elec')
    cfg.layout = data.elec;         % create layout from electrode definition
  else
    cfg.layout = 'CTF151s.lay';     % revert to the FCDC default
  end
end
lay = createlayout(cfg.layout);

% Select the channels in the data that match with the layout:
[seldat, sellay] = match_str(data.label, lay.label);
datavector = dat(seldat);
chanX = lay.prj(sellay,1);
chanY = lay.prj(sellay,2);
chanLabels = lay.label(sellay);

% Scale the channel locations between -0.45 and +0.45, which fits with the head:
chanX = 0.9*((chanX-min(chanX))/(max(chanX)-min(chanX))-0.5);
chanY = 0.9*((chanY-min(chanY))/(max(chanY)-min(chanY))-0.5);

hold on;

% Apply custom colormap:
if isfield(cfg, 'colormap')
	if size(cfg.colormap,2)~=3
		error('topoplot(): Colormap must be a n x 3 matrix');
	end
	colormap(cfg.colormap);
end

% Draw the head:
drawHead(cfg, rmax);

% Draw topoplot:
if ~strcmp(cfg.style, 'blank')
  drawTopoplot(cfg, datavector, rmax, chanX, chanY);
  % topoplot(cfg,chanX,chanY,datavector,chanLabels)
end

% Draw electrodes:
plot(chanX, chanY, cfg.emarker, ...
  'Color', cfg.ecolor, 'markersize', cfg.emarkersize);

if 0
% Write comment:
i = strmatch('COMNT', lay.label);
if length(i)==1
  comment = cfg.comment;
  if ~isempty(cfg.xparam)
    comment = sprintf('%0s\nxlim (%0s)=[%.3g %.3g]', comment, cfg.xparam, data.(cfg.xparam)(xmin), data.(cfg.xparam)(xmax));
  end
  if ~isempty(cfg.yparam)
    comment = sprintf('%0s\nylim (%0s)=[%.3g %.3g]', comment, cfg.yparam, data.(cfg.yparam)(ymin), data.(cfg.yparam)(ymax));
  end
  if ~isempty(cfg.zparam)
    comment = sprintf('%0s\nzlim (%0s)=[%.3g %.3g]', comment, cfg.zparam, cfg.zlim(1), cfg.zlim(2));
  end
  text(lay.prj(i,1), lay.prj(i,2), comment, 'Fontsize', cfg.fontsize);
end
end;

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  userData.hFigure = gcf;
  userData.hAxes = gca;
  for i=1:10
    userData.hSelection{i} = plot(0,0);
    set(userData.hSelection{i}, ...
      'XData', [0], ...
      'YData', [0], ...
      'Color', [0 0 0], ...
      'EraseMode', 'xor', ...
      'LineStyle', '--', ...
      'LineWidth', 1.5, ...
      'Visible', 'on');
    userData.range{i} = [];
  end
  userData.iSelection = 0;
  userData.plotType = 'topoplot';
  userData.selecting = 0;
  userData.selectionType = '';
  userData.selectAxes = 'z';
  userData.lastClick = [];
  userData.cfg = cfg;
  userData.data = data;
  userData.chanX = chanX;
  userData.chanY = chanY;
  userData.chanLabels = chanLabels;
  tag = sprintf('%.5f', 10000 * rand(1));
  set(gcf, 'Renderer', 'opengl');
  set(gcf, 'Tag', tag);
  set(gcf, 'UserData', userData);
  set(gcf, 'WindowButtonMotionFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 0);']);
  set(gcf, 'WindowButtonDownFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 1);']);
  set(gcf, 'WindowButtonUpFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 2);']);
end

axis off;
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
function drawHead(cfg, rmax)

l = 0:2*pi/100:2*pi;
basex = .18*rmax;
tip = rmax*1.15;
base = rmax-.004;
EarX = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
EarY = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];

% Plot head, ears, and nose:
plot(cos(l).*rmax, sin(l).*rmax, ...
  'Color', cfg.hcolor, 'Linestyle', '-', 'LineWidth', cfg.hlinewidth);
plot([.18*rmax;0;-.18*rmax], [base;tip;base], ...
  'Color', cfg.hcolor, 'LineWidth', cfg.hlinewidth);
plot(EarX, EarY, 'color', cfg.hcolor, 'LineWidth', cfg.hlinewidth);
plot(-EarX, EarY, 'color', cfg.hcolor, 'LineWidth', cfg.hlinewidth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawTopoplot(cfg, datavector, rmax, chanX, chanY)

% find limits for interpolation:
if strcmp(cfg.interplimits, 'head')
  xmin = min(-.5, min(chanX)); xmax = max(0.5, max(chanX));
  ymin = min(-.5, min(chanY)); ymax = max(0.5, max(chanY));
else
  xmin = max(-.5, min(chanX)); xmax = min(0.5, max(chanX));
  ymin = max(-.5, min(chanY)); ymax = min(0.5, max(chanY));
end

xi = linspace(xmin, xmax, cfg.gridscale);   % x-axis description (row vector)
yi = linspace(ymin, ymax, cfg.gridscale);   % y-axis description (row vector)

[Xi,Yi,Zi] = griddata(chanX, chanY, datavector, yi', xi, cfg.interpolation); % Interpolate data

% Take data within head:
mask = (sqrt(Xi.^2 + Yi.^2) <= rmax);
Zi(find(mask == 0)) = NaN;

% Calculate colormap limits:
m = size(colormap, 1);
if isstr(cfg.zlim)
  if strcmp(cfg.zlim, 'absmax')
    amin = -max(max(abs(Zi)));
    amax = max(max(abs(Zi)));
  elseif strcmp(cfg.zlim, 'maxmin')
    amin = min(min(Zi));
    amax = max(max(Zi));
  end
else
  amin = cfg.zlim(1);
  amax = cfg.zlim(2);
end
delta = xi(2)-xi(1); % length of grid entry

% Draw topoplot on head according to cfg.style:
if strcmp(cfg.style, 'contour')
  contour(Xi, Yi, Zi, cfg.contournum, 'k');
elseif strcmp(cfg.style, 'both')
  surface(Xi-delta/2, Yi-delta/2, zeros(size(Zi)), Zi, 'EdgeColor', 'none',...
    'FaceColor', cfg.shading);
  contour(Xi, Yi, Zi, cfg.contournum, 'k');
elseif strcmp(cfg.style, 'straight')
  surface(Xi-delta/2, Yi-delta/2, zeros(size(Zi)), Zi, 'EdgeColor', 'none',...
    'FaceColor', cfg.shading);
elseif strcmp(cfg.style, 'fill')
  contourf(Xi, Yi, Zi, cfg.contournum,'k');
else
  error(sprintf('Invalid style (cfg.style=''%0s'')', cfg.style));
end

% Set color axis:
caxis([amin amax]);
