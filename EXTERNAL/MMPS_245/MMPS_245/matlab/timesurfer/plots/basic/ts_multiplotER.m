function ts_multiplotER(cfg, varargin)

% multiplotER makes plots of event-related fields or potentials versus time
% or of oscillatory activity (power or coherence) versus frequency. Multiple
% datasets can be overlayed.  The plots are arranged according to their
% location specified in the LAYOUTFILE.
%
% Use as:
%   multiplotER(cfg, data)
%   multiplotER(cfg, data, data2, ..., dataN)
%
% The data can be an ERP/ERF produced by TIMELOCKANALYSIS, a
% powerspectrum produced by FREQANALYSIS or a coherencespectrum produced
% by FREQDESCRIPTIVES. If you specify multiple datasets they must contain
% the same channels, etc.
%
% The cfg structure contains all configurations details for the plotting
%
% cfg.layout        = specification of the layout, see below
% cfg.xparam        = field to be plotted on x-axis      (default = 'time')
% cfg.yparam        = field to be plotted on y-axis      (default = 'avg')
% cfg.xlim          = 'maxmin' or [xmin xmax]            (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax]            (default = 'maxmin')
% cfg.baseline      = 'yes','no' or [time1 time2]        (default = 'no')
% cfg.comment       = string of text                     (default = date + colors)
%                     Add 'comment' to graph (according to COMNT in the layout)
% cfg.axes          = 'yes', 'no'                        (default = 'yes')
%                     Draw x- and y-axes for each graph
% cfg.box           = 'yes', 'no'                        (default = 'no')
%                     Draw a box around each graph
% cfg.showlabels    = 'yes', 'no'                        (default = 'no')
% cfg.fontsize      = size of font                       (default = 8)
% cfg.linewidth     = trace line width                   (default = 1.5)
% cfg.interactive   = 'yes', 'no'                        (default = 'no')
%
% If you use the output of freqdescriptives as an input and if you want to visualise
% coherence with respect to a reference-channel, you should specify:
%
% cfg.yparam        =  'cohspctrm'
% cfg.cohrefchannel =  name of reference-channel
%
% The layout defines how the channels are arranged and what the size of each
% subplot is. You can specify the layout in a variety of ways:
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
%   multiplotTFR, singleplotER, singleplotTFR, topoplotER, topoplotTFR.

% Copyright (C) 2003, Ole Jensen
%
% $Log: multiplotER.m,v $
% Revision 1.15  2006/03/02 13:54:47  jansch
% fixed multiple small bugs
%
% Revision 1.14  2006/02/28 12:43:15  roboos
% made plotting of coherence consistent between all xxxplotER functions
% made baselining consistent, use cfg.xparam to decide between freqbaseline and timelockbaseline
%
% Revision 1.13  2006/02/27 15:03:03  denpas
% many changes, most important is added interactive functionality
% made data selection consistent between different plot functions
% changed dimord for consistency
%
% Revision 1.11  2005/08/18 12:15:41  jansch
% added support to plot subfields
%
% Revision 1.10  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.9  2005/04/29 12:46:44  roboos
% small change in help
%
% Revision 1.8  2005/04/29 12:40:43  roboos
% cleaned up and updated the help
%
% Revision 1.7  2005/04/12 10:06:57  olejen
% clf added
%
% Revision 1.6  2005/04/06 07:40:56  jansch
% included option cfg.cohrefchannel. updated help.
%
% Revision 1.5  2005/02/07 17:12:00  roboos
% changed handling of layout files (using new function createlayout), now also supports automatic layout creation based on gradiometer/electrode definition in data, updated help, cleaned up indentation
%
% Revision 1.4  2005/01/27 09:31:49  roboos
% applied autoindentation on code, removed many empty lines and spaces,
% replaced layoutfile reading with read_lay, applied doudavs code to all input arguments
% implemented automatic detection of arguments to plot,
% updated help, changed input from p1, p2, p3... to varargin
%
% Revision 1.3  2004/09/24 15:54:54  roboos
% included the suggested improvements by Doug Davidson: added option cfg.cohtargetchannel
% and updated the help
%
% Revision 1.2  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details

if nargin<2
  help(mfilename);
  return;
end;

clf

%% todo: initially make plots small, then resize
%  (trick to avoid plots getting deleted)

% For backward compatibility with old data structures:
for i=1:length(varargin)
  varargin{i} = fixdimord(varargin{i});
end

% set the defaults:
if ~isfield(cfg,'baseline'),    cfg.baseline    = 'no';                        end
if ~isfield(cfg,'xlim'),        cfg.xlim        = 'maxmin';                    end
if ~isfield(cfg,'ylim'),        cfg.ylim        = 'maxmin';                    end
if ~isfield(cfg,'comment'),     cfg.comment     = strcat([date '\n']);         end
if ~isfield(cfg,'axes'),        cfg.axes        = 'yes';                       end
if ~isfield(cfg,'showlabels'),  cfg.showlabels  = 'no';                        end
if ~isfield(cfg,'box'),         cfg.box         = 'no';                        end
if ~isfield(cfg,'fontsize'),    cfg.fontsize    = 8;                           end
if ~isfield(cfg,'linewidth'),   cfg.linewidth   = 1.5;                         end
if ~isfield(cfg,'graphcolor')   cfg.graphcolor  = ['brgkywrgbkywrgbkywrgbkyw'];end
if ~isfield(cfg,'interactive'), cfg.interactive = 'no';                        end

GRAPHCOLOR = ['k' cfg.graphcolor ];

% Set x/y/zparam defaults according to varargin{1}.dimord value:
if strcmp(varargin{1}.dimord, 'chan_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';                      end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='avg';                   end
elseif strcmp(varargin{1}.dimord, 'chan_freq')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='freq';                  end
  if ~isfield(cfg, 'yparam'),      cfg.yparam='';                      end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
end

% Old style coherence plotting with cohtargetchannel is no longer supported:
if isfield(cfg,'cohtargetchannel'), 
  error('cfg.cohtargetchannel is obsolete, check the documentation for help about coherence plotting.');
end

for k=1:length(varargin)
  % Check for unconverted coherence spectrum data:
  if (strcmp(cfg.zparam,'cohspctrm')) & (isfield(varargin{k}, 'labelcmb'))
    % A reference channel is required:
    if ~isfield(cfg,'cohrefchannel'),
      error('no reference channel specified');
    end
    % Convert 2-dimensional channel matrix to a single dimension:
    sel1                  = strmatch(cfg.cohrefchannel, varargin{k}.labelcmb(:,2));
    sel2                  = strmatch(cfg.cohrefchannel, varargin{k}.labelcmb(:,1));
    fprintf('selected %d channels for coherence\n', length(sel1)+length(sel2));
    varargin{k}.cohspctrm = varargin{k}.cohspctrm([sel1;sel2],:,:);
    varargin{k}.labelcmb  = varargin{k}.labelcmb([sel1;sel2],:);
    varargin{k}.label     = [varargin{k}.labelcmb(sel1,1);varargin{k}.labelcmb(sel2,2)];
    varargin{k}           = rmfield(varargin{k}, 'labelcmb');
  end

  % Apply baseline correction:
  if ~strcmp(cfg.baseline, 'no')
    if strcmp(cfg.yparam, 'time')
      varargin{k} = timelockbaseline(cfg, varargin{k});
    elseif strcmp(cfg.yparam, 'freq')
      varargin{k} = freqbaseline(cfg, varargin{k});
    end
  end
end

% Get physical x-axis range:
if strcmp(cfg.xlim,'maxmin')
  % Find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:length(varargin)
    xmin = min([xmin varargin{i}.(cfg.xparam)]);
    xmax = max([xmax varargin{i}.(cfg.xparam)]);
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Find corresponding x-axis bins:
xidc = find(varargin{1}.(cfg.xparam) >= xmin & varargin{1}.(cfg.xparam) <= xmax);

% Align physical x-axis range to the array bins:
xmin = varargin{1}.(cfg.xparam)(xidc(1));
xmax = varargin{1}.(cfg.xparam)(xidc(end));


% Get physical y-axis range (ylim / zparam):
if strcmp(cfg.ylim,'maxmin')
  % Find maxmin throughout all varargins:
  ymin = [];
  ymax = [];
  for i=1:length(varargin)
    % dh - fixed bug here 06/07/06
    ymin = min([ymin;varargin{i}.(cfg.zparam)(:)]);
    ymax = max([ymax;varargin{i}.(cfg.zparam)(:)]);
  end
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% Read or create the layout that will be used for plotting:
if ~isfield(cfg, 'layout')
  if isfield(cfg, 'layoutname')
    cfg.layout = cfg.layoutname;       % backward compatible cfg fieldname
  elseif isfield(varargin{1}, 'grad')
    cfg.layout = varargin{1}.grad;            % create layout from gradiometer definition
  elseif isfield(varargin{1}, 'elec')
    cfg.layout = varargin{1}.elec;            % create layout from electrode definition
  else
    cfg.layout = 'CTF151s.lay';        % revert to the FCDC default
  end
end

lay = createlayout(cfg.layout);

% convert the layout to Ole's style of variable names
X      = lay.prj(:,1);
Y      = lay.prj(:,2);
Width  = lay.box(:,1);
Height = lay.box(:,2);
Lbl    = lay.label;

% Create empty channel coordinates and labels arrays:
chanX(1:length(Lbl)) = NaN;
chanY(1:length(Lbl)) = NaN;
chanLabels = cell(1,length(Lbl));

hold on;
colorLabels = [];

% Plot each data set:
for k=1:length(varargin)
  P          = getsubfield(varargin{k}, cfg.zparam);
  Labels     = getfield(varargin{k}, 'label');

  if length(varargin) > 1
%    colorLabels = [colorLabels inputname(k+1) '=' GRAPHCOLOR(k+1) '\n'];
    colorLabels = [colorLabels sprintf('cond-%d',k) '=' GRAPHCOLOR(k+1) '\n'];
  end

  style = GRAPHCOLOR(k+1);
  
  for m=1:length(Lbl)
    l = cellstrmatch(Lbl(m),Labels);
    if ~isempty(l)
      % Plot ER:  
      plotWnd(varargin{k}.(cfg.xparam),P(l,:),xidc,[xmin xmax],[ymin ymax], ...
        X(m), ...
        Y(m), ...
        Width(m), ...
        Height(m), ...
        Lbl(m), ...
        cfg,style);

      % Keep ER plot coordinates (at centre of ER plot), and channel labels (will be stored in the figure's UserData struct):
      chanX(m) = X(m) + 0.5 * Width(m);
      chanY(m) = Y(m) + 0.5 * Height(m);
      chanLabels{m} = Lbl{m};
    end
  end
end

% Add the colors of the different datasets to the comment:
cfg.comment = [cfg.comment colorLabels];

% Write comment text:
l = cellstrmatch('COMNT',Lbl);
if ~isempty(l)
  text(X(l),Y(l),sprintf(cfg.comment),'Fontsize',cfg.fontsize);
end

% Plot scales:
l = cellstrmatch('SCALE',Lbl);
if ~isempty(l)
  plotScales([xmin xmax],[ymin ymax],X(l),Y(l),Width(1),Height(1),cfg)
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  userData.hFigure = gcf;
  userData.hAxes = gca;
  for i=1:10
    userData.hSelection{i} = plot(mean(X), mean(Y));
    set(userData.hSelection{i}, 'XData', [mean(X)]);
    set(userData.hSelection{i}, 'YData', [mean(Y)]);
    set(userData.hSelection{i}, 'Color', [0 0 0]);
    set(userData.hSelection{i}, 'EraseMode', 'xor');
%    set(userData.hSelection{i}, 'LineStyle', '--');
    set(userData.hSelection{i}, 'LineStyle', '-');
    set(userData.hSelection{i}, 'LineWidth', cfg.linewidth)
    set(userData.hSelection{i}, 'Visible', 'on');
    userData.range{i} = [];
  end
  userData.iSelection = 0;
  userData.plotType = 'multiplot';
  userData.selecting = 0;
  userData.selectionType = '';
  userData.selectAxes = 'z';
  userData.lastClick = [];
  userData.cfg = cfg;
  userData.data = varargin;
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

axis tight
axis off
orient landscape
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScales(xlim,ylim,xpos,ypos,width,height,cfg)
x1 =  xpos;
x2 =  xpos+width;
y1 =  ypos;
y2 =  ypos+width;
plot([xpos xpos+width xpos+width xpos xpos],[ypos ypos ypos+height ypos+height ypos],'b');
if xlim(1) <=  0 && xlim(2) >= 0
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  plot(xs,ys,'b');
end

if ylim(1) <= 0 && ylim(2) >= 0
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  plot(xs,ys,'b');
end

text( x1,y1,num2str(xlim(1),3),'rotation',90,'HorizontalAlignment','Right','VerticalAlignment','middle','Fontsize',cfg.fontsize);
text( x2,y1,num2str(xlim(2),3),'rotation',90,'HorizontalAlignment','Right','VerticalAlignment','middle','Fontsize',cfg.fontsize);
text( x2,y1,num2str(ylim(1),3),'HorizontalAlignment','Left','VerticalAlignment','bottom','Fontsize',cfg.fontsize);
text( x2,y2,num2str(ylim(2),3),'HorizontalAlignment','Left','VerticalAlignment','bottom','Fontsize',cfg.fontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotWnd(x,y,xidc,xlim,ylim,xpos,ypos,width,height,label,cfg,style)
set(gca,'FontSize',cfg.fontsize);

x = x(xidc);
y = y(xidc);

% Clip out of bounds y values:
y(find(y > 1.5*ylim(2))) = 1.5*ylim(2);
y(find(y < 1.5*ylim(1))) = 1.5*ylim(1);

xs = xpos+width*(x-xlim(1))/(xlim(2)-xlim(1));
ys = ypos+height*(y-ylim(1))/(ylim(2)-ylim(1));
plot(xs,ys,style,'LineWidth',cfg.linewidth)

if strcmp(cfg.showlabels,'yes')
  text(xpos,ypos+1.0*height,label,'Fontsize',cfg.fontsize)
end

% Draw axes:
if strcmp(cfg.axes,'yes')
  % Draw y axis if xlim crosses 0:
  if (xlim(1) <= 0) & (xlim(2) >= 0)
    xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
    ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
    plot(xs,ys,'k');
  end
  
  % Draw x axis if ylim crosses 0:
  if (ylim(1) <= 0) & (ylim(2) >= 0)
    xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
    ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
    plot(xs,ys,'k');
  end
end

% Draw box around plot:
if strcmp(cfg.box,'yes')
  plot([xpos xpos+width xpos+width xpos xpos],[ypos ypos ypos+height ypos+height ypos],'b');
end

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
function z = pwrspctm2cohspctrm(cfg, freq)
% This creates a copy of freq with new entries in the freq.pwrspctrm. The
% TFR corresponding to cfg.cohtargetchannel is set to ones, but all other
% channels are replaced with the coherence TFRs from freq.cohspctrm.
%
% [z] = pwrspctrm2cohspctrm(cfg, freq);
%
% freq should be organised in a structure as obtained from
% the freqdescriptives function.
%
% doudav; 24.09.04
z = freq;
[a a1] = match_str( cfg.cohtargetchannel, z.label);
[a b1] = match_str( cfg.cohtargetchannel, z.labelcmb(:,1));
[a b2] = match_str( cfg.cohtargetchannel, z.labelcmb(:,2));
[a c1] = match_str( z.labelcmb(b1,2), z.label);
[a c2] = match_str( z.labelcmb(b2,1), z.label);
targets  = [c1; c2];
targets2 = [b1; b2];
z.powspctrm(a1,:,:) = ones(size(z.powspctrm(a1,:,:)));
for i=1:size(targets,1),
  z.powspctrm(targets(i),:,:) = squeeze(z.cohspctrm(targets2(i),:,:));
end
clear a a1 b1 b2 c1 c2;
