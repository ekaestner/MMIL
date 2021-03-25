function [handle] = topoplot(varargin)

% TOPOPLOT plots a topographic map of an EEG or MEG field as a 2-D
% circular view (looking down at the top of the head) using interpolation
% on a fine cartesian grid. 
%
% This function is called by topoplotER or topoplotTFR
%
% You can also call this function directly as follows:
%         topoplot(cfg, datavector)
%         topoplot(cfg, X, Y, datavector)
%         topoplot(cfg, X, Y, datavector, Labels)
%         topoplot(datavector,'Key1','Value1','Key2','Value2',...)
%         
% Inputs can be either:
%     datavector  = vector of values to be plotted as color
%     cfg         = configuration structure containing the (optional) parameters
%     X           = x-coordinates for channels in datavector
%     Y           = y-coordinates for channels in datavector
%     Labels      = labels for channels in datavector
% or the inputs can be key-value pairs containing the (optional) parameters. 
% Every cfg field can be specified using the fieldname as a key in the
% key-value pairs.
%
% if X, Y and Labels are given, cfg.layout is NOT used. If X, Y, and Labels
% are not given, cfg.layout must be given and it is assumed that the
% channels in the datavector exactly mach the channels in the layout. 
%
% The layout defines how the channels will be arranged in the 2-D plane. 
% You can specify the layout in a variety of ways:
%  - you can give the name of an ascii layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition,  i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these, and if the data structure contains an
% electrode or gradiometer structure, that will be used for creating a
% layout.
% 
% Optional Parameters and Values
% 
% cfg.colormap        = any sized colormap, see COLORMAP
% cfg.colorbar        = 'yes'
%                       'no' (default)
%                       'North'              inside plot box near top
%                       'South'              inside bottom
%                       'East'               inside right
%                       'West'               inside left
%                       'NorthOutside'       outside plot box near top
%                       'SouthOutside'       outside bottom
%                       'EastOutside'        outside right
%                       'WestOutside'        outside left
% cfg.interplimits    = limits for interpolation (default = 'head')
%                       'electrodes' to furthest electrode
%                       'head' to edge of head
% cfg.gridscale       = scaling grid size (default = 67)
%                       determines resolution of figure   
% cfg.maplimits       = 'absmax' +/- the absolute-max (default = 'absmax')
%                       'maxmin' scale to data range
%                       [clim1, clim2] user-defined lo/hi
% cfg.style           = topoplot style (default = 'both')
%                       'straight' colormap only
%                       'contour' contour lines only
%                       'both' (default) both colormap and contour lines
%                       'fill' constant color between lines
%                       'blank' just head and electrodes
% cfg.contournum      = number of contour lines (default = 6), see CONTOUR
% cfg.shading         = 'flat' 'interp' (default = 'flat')
% cfg.interpolation   = 'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
% cfg.headcolor       = Color of head cartoon (default = [0,0,0])
% cfg.hlinewidth      = number, Linewidth of the drawn head, nose and ears (default = 2) 		
% cfg.contcolor       = Contourline color (default = [0 0 0])	
% cfg.electrodes      = 'on','off','labels','numbers','highlights' or 'dotnum' (default = 'on')
% cfg.emarker         = Marker symbol (default = 'o')		
% cfg.ecolor          = Marker color (default = [0 0 0] (black))
% cfg.emarkersize     = Marker size (default = 2)	
% cfg.efontsize       = Font size of electrode labels/numbers (default = 8 pt)
%                       when cfg.electrodes = 'numbers' or 'labels'
% cfg.comment         =  string of text 
% cfg.commentpos      = position of comment (default = 'leftbottom')
%                       'lefttop' 'leftbottom' 'middletop' 'middlebottom' 'righttop' 'rightbottom'
%                       or [x y] coordinates 
%                       or 'title' to place comment as title 
% cfg.fontsize        = Font size of comment (default = 8 pt)
% cfg.highlight       = 'off' or the channel numbers you want to highlight (default = 'off').
%                       These numbers should correspond with the channels in the data, not in
%                       the layout file.
% cfg.hlmarker        = Highlight marker symbol (default = 'o')  
% cfg.hlcolor         = Highlight marker color (default = [0 0 0] (black)) 
% cfg.hlmarkersize    = Highlight marker size (default = 6) 
% cfg.hllinewidth     = Highlight marker linewidth (default = 3) 
% cfg.outline         = 'scalp' or 'ECog' (default = 'scalp')
% 
% Note: topoplot() only works when map limits are >= the max and min 
%                             interpolated data values.

% Undocumented local options:
% cfg.efsize
% cfg.electcolor
% cfg.electrod
% cfg.emsize
% cfg.grid
% cfg.hcolor
% cfg.headlimits
% cfg.interpolate
% cfg.maxchans
% cfg.showlabels
% cfg.zlim
% cfg.mask for opacity masking, e.g. with statistical significance

% Copyright (C) 1996, Andy Spydell, Colin Humphries & Arnaud Delorme, CNL / Salk Institute
% Copyright (C) 2004-2006, F.C. Donders Centre
% New implementation by Geerten Kramer, based on versions of Ole Jensen and Jan-Mathijs Schoffelen
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: topoplot.m,v $
% Revision 1.33  2008/05/29 13:45:40  roboos
% added hack for Miriam, should be finished by ingnie
%
% Revision 1.32  2008/04/24 10:29:38  roboos
% Added cfg.outline, which can be used to toggle between scalp-mode (default) and ecog mode. This affects the masking and the outline of the head.
%
% Revision 1.31  2007/10/29 16:08:22  marvger
% added the possibility to define the colorbar location (with documentation)
%
% Revision 1.30  2007/10/29 16:01:00  marvger
% suppressed colorbar when all data is equal (this led to an error)
%
% Revision 1.29  2007/10/24 06:13:08  roboos
% improved the documentation for highlights, thanks to Nicholas
%
% Revision 1.28  2007/08/08 07:02:25  roboos
% extended the detection of eeglab-style inputs, see mail from and to Arno on 8 Aug 2007
%
% Revision 1.27  2007/07/24 16:23:49  ingnie
% fixed cfg.contcolor, did not work previously now only works when is string
%
% Revision 1.26  2007/03/21 12:44:29  roboos
% added detection of EEGLAB-style input arguments and give a long error message that explains the path setting
%
% Revision 1.25  2007/03/14 08:43:12  roboos
% replaced call to createlayout to prepare_layout, made some small changes to the lay structure

% Try to detect EEGLAB-style input and give an informative error
% message. The EEGLAB documentation describes the usage as 
%        >>  topoplot(datavector, EEG.chanlocs);   % plot a map using an EEG chanlocs structure
%        >>  topoplot(datavector, 'my_chan.locs'); % read a channel locations file and plot a map
%        >>  topoplot('example');                  % give an example of an electrode location file
%        >>  [h grid_or_val plotrad_or_grid, xmesh, ymesh]= ...
%                           topoplot(datavector, chan_locs, 'Input1','Value1', ...);

if nargin==2 && isvector(varargin{1}) && isstruct(varargin{2}) && isfield(varargin{2}, 'labels')
  eeglab = 1;
elseif nargin==2 && isvector(varargin{1}) && ischar(varargin{2})
  eeglab = 1;
elseif nargin==1 && isequal(varargin{1}, 'example')
  eeglab = 1;
elseif nargin>2 && isvector(varargin{1}) && mod(nargin,2)==0 && isstruct(varargin{2}) && isfield(varargin{2}, 'labels')
  eeglab = 1;
else
  eeglab = 0;
end

if eeglab
  % the input resembles the input of the EEGLAB topoplot function
  error('Unrecognized input, please look at "help topoplot", "which topoplot" and "path". The input looks as if you expect the EEGLAB version of the topoplot function. Your path settings may be incorrect and the FieldTrip and EEGLAB version of the "topoplot" function may be confused.')
end

% deal with the different types of input syntax that this function can get
if mod(nargin,2) && isnumeric(varargin{1}) && isstr(varargin{2})
  % topoplot(data, key, val, ...)
  cfg  = keyval2cfg(varargin(2:end));
  data = varargin{1};
  OldStyleSyntax=0;
elseif nargin==2
  % topoplot(cfg,data)
  OldStyleSyntax=0;
  cfg  = varargin{1};
  data = varargin{2};
  err  = 0;
  if ~isempty(cfg),
    err  = err + ~isstruct(cfg);
  end
  err  = err + ~isnumeric(data);
  if err
    errmsg=['\n'];
    errmsg=[errmsg,'When two input arguments are supplied, the following syntax should be used:\n'];
    errmsg=[errmsg,'topoplot(cfg,datavector);\n'];
    error(sprintf(errmsg));
  end;
elseif nargin==4
  % topoplot(cfg,X,Y,data)
  OldStyleSyntax=1;
  cfg  = varargin{1};
  X    = varargin{2};
  Y    = varargin{3};
  data = varargin{4};
  err  = 0;
  if ~isempty(cfg),
    err  = err + ~isstruct(cfg);
  end
  err  = err + ~isnumeric(data);
  err  = err + ~isnumeric(X);
  err  = err + ~isnumeric(Y);
  if err
    errmsg=['\n'];
    errmsg=[errmsg,'When four input arguments are supplied, the following syntax should be used:\n'];
    errmsg=[errmsg,'topoplot(cfg,X,Y,datavector);\n'];
    error(sprintf(errmsg));
  end;
elseif nargin==5
  % topoplot(cfg,X,Y,data,labels)
  OldStyleSyntax=1;
  cfg    = varargin{1};
  X      = varargin{2};
  Y      = varargin{3};
  data   = varargin{4};
  labels = varargin{5};
  err    = 0;
  if ~isempty(cfg), 
    err  = err + ~isstruct(cfg);
  end
  err    = err + ~isnumeric(data);
  err    = err + ~isnumeric(X);
  err    = err + ~isnumeric(Y);
  err    = err + ~iscell(labels);
  if err
    errmsg=['\n'];
    errmsg=[errmsg,'When five input arguments are supplied, the following syntax should be used:\n'];
    errmsg=[errmsg,'topoplot(cfg,X,Y,datavector,Labels);\n'];
    error(sprintf(errmsg));
  end;
else
  error('unrecognized input, please look at the help of this function')
end

% set the defaults
if ~isfield(cfg, 'maxchans')      cfg.maxchans = 256;       end;
if ~isfield(cfg, 'maplimits')     cfg.maplimits = 'absmax'; end; % absmax, maxmin, [values]
if ~isfield(cfg, 'interplimits')  cfg.interplimits ='head'; end; % head, electrodes
if ~isfield(cfg, 'grid_scale')    cfg.grid_scale = 67;      end; % 67 in original
if ~isfield(cfg, 'contournum')    cfg.contournum = 6;       end;
if ~isfield(cfg, 'colorbar')      cfg.colorbar = 'no';      end;
if ~isfield(cfg, 'style')         cfg.style = 'both';       end; % both,straight,fill,contour,blank
if ~isfield(cfg, 'hcolor')        cfg.hcolor = [0 0 0];     end;
if ~isfield(cfg, 'contcolor')     cfg.contcolor = 'k';      end;
if ~isfield(cfg, 'hlinewidth')    cfg.hlinewidth = 2;       end;
if ~isfield(cfg, 'shading')       cfg.shading = 'flat';     end; % flat or interp
if ~isfield(cfg, 'interpolation') cfg.interpolation = 'v4'; end;
if ~isfield(cfg, 'fontsize'),     cfg.fontsize = 8;         end;
if ~isfield(cfg, 'commentpos'),   cfg.commentpos = 'leftbottom';    end;
if ~isfield(cfg, 'mask'),         cfg.mask = [];            end;
if ~isfield(cfg, 'outline'),      cfg.outline = 'scalp';    end; % scalp or ecog

if ~isfield(cfg,'layout')
  if ~OldStyleSyntax
    error('Specify at least the field or key "layout".');
  end;
end;
if ~isstr(cfg.contcolor)     cfg.contcolor = 'k'; warning('cfg.contcolor must be string, put to ''k''');   end;

if ~isfield(cfg,'electrodes')   cfg.electrodes = 'on'; end; % on,off,label,numbers or highlights
if ~isfield(cfg,'showlabels') % for compatibility with OLDSTYLE
  cfg.showlabels = '';
else
  cfg.electrodes = '';
end; 

if ~isfield(cfg,'emarker')      cfg.emarker = 'o';     end;
if ~isfield(cfg,'ecolor')       cfg.ecolor = [0 0 0];  end;
if ~isfield(cfg,'emarkersize')  cfg.emarkersize = 2;   end;
if ~isfield(cfg,'efsize')       cfg.efsize = get(0,'DefaultAxesFontSize');end;

if ~isfield(cfg,'highlight')    cfg.highlight = 'off'; end; % 'off' or the electrodenumbers.
if ~isfield(cfg,'hlmarker')     cfg.hlmarker = 'o';    end;
if ~isfield(cfg,'hlcolor')      cfg.hlcolor = [0 0 0]; end;
if ~isfield(cfg,'hlmarkersize') cfg.hlmarkersize = 6;  end;
if ~isfield(cfg,'hllinewidth')	cfg.hllinewidth = 3;   end;

if isfield(cfg,'colormap')
  if size(cfg.colormap,2)~=3, error('topoplot(): Colormap must be a n x 3 matrix'); end
  colormap(cfg.colormap);
end;

if isfield(cfg,'headlimits')
  cfg.interplimits = cfg.headlimits;
  cfg              = rmfield(cfg,'headlimits');
  if ~isstr(cfg.interplimits), error('topoplot(): interplimits value must be a string'); end
  cfg.interplimits = lower(cfg.interplimits);
  if ~strcmp(cfg.interplimits,'electrodes') & ~strcmp(cfg.interplimits,'head'),
    error('topoplot(): Incorrect value for interplimits');
  end
end;
	
if isfield(cfg,'gridscale')
  cfg.grid_scale = cfg.gridscale;
  cfg            = rmfield(cfg,'gridscale');
end;

if isfield(cfg,'interpolate')
  cfg.interpolation = lower(cfg.interpolate);
  cfg               = rmfield(cfg,'interpolate');
end;

if isfield(cfg,'numcontour')
  cfg.contournum = cfg.numcontour;
  cfg            = rmfield(cfg,'numcontour');
end;

if isfield(cfg,'electrod')
	cfg.electrodes = lower(cfg.electrod);
	cfg            = rmfield(cfg,'electrod');
end;

if isfield(cfg,'headcolor')
  cfg.hcolor = cfg.headcolor;
  cfg        = rmfield(cfg,'headcolor');
end;

if isfield(cfg,'electcolor')
  cfg.ecolor = cfg.electcolor;
  cfg        = rmfield(cfg,'electcolor');
end;

if isfield(cfg,'emsize')
  cfg.emarkersize = cfg.emsize;
  cfg             = rmfield(cfg,'emsize');
end;

if isfield(cfg,'efontsize') 
  cfg.efsize = cfg.efontsize;
  cfg        = rmfield(cfg,'efontsize');
end;

if isfield(cfg,'shading')
  cfg.shading = lower(cfg.shading);
  if ~any(strcmp(cfg.shading,{'flat','interp'})), error('Invalid Shading Parameter'); end
end

if isfield(cfg,'zlim') 
  cfg.maplimits = cfg.zlim;
  cfg           = rmfield(cfg,'zlim');
end;

[r,c] = size(data);
if r>1 && c>1 ,
  error('topoplot(): data should be a single vector\n');
end

% create layout from cfg.layout to find matching X and Y coordinates and labels to the datavector
if ~OldStyleSyntax
  lay = prepare_layout(cfg);
  X = lay.pos(:,1);
  Y = lay.pos(:,2);
  labels = lay.label;
end
if exist('labels', 'var'),
  ind_SCALE = strmatch('SCALE', labels);
  ind_COMNT = strmatch('COMNT', labels);
  if length(ind_SCALE) == 1;
    X_SCALE           = X(ind_SCALE); 
    Y_SCALE           = Y(ind_SCALE);
    X(ind_SCALE)      = [];
    Y(ind_SCALE)      = [];
    labels(ind_SCALE) = [];
  end
  if length(ind_COMNT) == 1;
    ind_COMNT         = strmatch('COMNT', labels); %again; may be changed
    X_COMNT           = X(ind_COMNT);
    Y_COMNT           = Y(ind_COMNT);
    X(ind_COMNT)      = [];
    Y(ind_COMNT)      = [];
    labels(ind_COMNT) = [];
  end
end
if length(data)~=length(X)
  error('topoplot(): data vector must be the same size as layout-file')
end

if strcmpi(cfg.outline, 'scalp')
  % Scale the data to a circle with x-axis and y-axis: -0.45 to 0.45
  y = 0.9*((X-min(X))/(max(X)-min(X))-0.5); %ATTENTION x becomes y and y becomes x, in griddata is also reversed which makes x x and y y again
  x = 0.9*((Y-min(Y))/(max(Y)-min(Y))-0.5);
elseif strcmpi(cfg.outline, 'ecog') || strcmp(cfg.outline, 'miriam')
  y = X; %ATTENTION x becomes y and y becomes x, in griddata is also reversed which makes x x and y y again
  x = Y;
end

% Set coordinates for comment
if strcmp(cfg.commentpos,'lefttop') 
  x_COMNT = -0.7; 
  y_COMNT =  0.6;
  HorAlign = 'left';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'leftbottom') 
  x_COMNT = -0.6; 
  y_COMNT = -0.6;
  HorAlign = 'left';
  VerAlign = 'bottom';
elseif strcmp(cfg.commentpos,'middletop') 
  x_COMNT =  0; 
  y_COMNT =  0.75;
  HorAlign = 'center';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'middlebottom') 
  x_COMNT =  0; 
  y_COMNT = -0.7;
  HorAlign = 'center';
  VerAlign = 'bottom';
elseif strcmp(cfg.commentpos,'righttop') 
  x_COMNT =  0.65; 
  y_COMNT =  0.6;
  HorAlign = 'right';
  VerAlign = 'top';
elseif strcmp(cfg.commentpos,'rightbottom') 
  x_COMNT =  0.6; 
  y_COMNT = -0.6;
  HorAlign = 'right';
  VerAlign = 'bottom';
elseif isnumeric(cfg.commentpos)
  x_COMNT = cfg.commentpos(1);
  y_COMNT = cfg.commentpos(2);
  HorAlign = 'left';
  VerAlign = 'middle';
  x_COMNT = 0.9*((x_COMNT-min(X))/(max(X)-min(X))-0.5);
  y_COMNT = 0.9*((y_COMNT-min(Y))/(max(Y)-min(Y))-0.5);
end

rmax = .5;

ha = gca;
cla
hold on

if ~strcmp(cfg.style,'blank')
  % find limits for interpolation:
  if strcmp(cfg.interplimits,'head')
    xmin = min(-.5,min(x)); xmax = max(0.5,max(x));
    ymin = min(-.5,min(y)); ymax = max(0.5,max(y));
  else
    xmin = max(-.5,min(x)); xmax = min(0.5,max(x));
    ymin = max(-.5,min(y)); ymax = min(0.5,max(y));
  end

  xi         = linspace(xmin,xmax,cfg.grid_scale);   % x-axis description (row vector)
  yi         = linspace(ymin,ymax,cfg.grid_scale);   % y-axis description (row vector)
  [Xi,Yi,Zi] = griddata(y, x, data, yi', xi, cfg.interpolation); % Interpolate data
  % [Xi,Yi,Zi] = griddata(y,x,data,yi',xi,'invdist'); % Interpolate data
%keyboard  
  if strcmpi(cfg.outline, 'scalp') || strcmp(cfg.outline, 'miriam')
    % Take data within head
    mask   = (sqrt(Xi.^2+Yi.^2) <= rmax);
    ii     = find(mask == 0);
    Zi(ii) = NaN;
%		fprintf('[min max] of interpolated data: [%g %g]\n',min(Zi),max(Zi));
  elseif strcmpi(cfg.outline, 'ecog')
    % FIXME masking of ECoG data with the outline of the electrode grid is not yet implemented
    warning('masking of ECoG data with the outline of the electrode grid is not yet implemented');
  end

  % calculate colormap limits
  m = size(colormap,1);
  if isstr(cfg.maplimits)
    if strcmp(cfg.maplimits,'absmax')
      amin = -max(max(abs(Zi)));
      amax = max(max(abs(Zi)));
    elseif strcmp(cfg.maplimits,'maxmin')
      amin = min(min(Zi));
      amax = max(max(Zi));
    end
  else
    amin = cfg.maplimits(1);
    amax = cfg.maplimits(2);
  end
  delta = xi(2)-xi(1); % length of grid entry
%keyboard
  if ~isempty(cfg.mask),
    % create mask, if requested
    [maskX,maskY,maskZ] = griddata(y, x, double(cfg.mask), yi', xi, cfg.interpolation);
    % mask should be scaled between 0 and 1, clip the values that ly outside that range
    maskZ(isnan(maskZ)) = 0;
    maskZ(isinf(maskZ)) = 0;
    maskZ(maskZ<0) = 0;
    maskZ(maskZ>1) = 1;
  end 
%keyboard
  % Draw topoplot on head
  if strcmp(cfg.style,'contour')
    contour(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
  elseif strcmp(cfg.style,'both')
    contour(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
    h = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',cfg.shading);
    if exist('maskZ','var'),
      set(h, 'AlphaData', maskZ);
      alim([0 1]);
      set(h, 'FaceAlpha', 'interp');
    end
  elseif strcmp(cfg.style,'straight')
    h = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',cfg.shading);
    if exist('maskZ','var'),
      set(h, 'AlphaData', maskZ);
      alim([0 1]);
      set(h, 'FaceAlpha', 'interp');
    end
  elseif strcmp(cfg.style,'fill')
    contourf(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
  else
    error('Invalid style')
  end
  caxis([amin amax]) % set coloraxis
end

% Plot electrodes:
if strcmp(cfg.electrodes,'on')||strcmp(cfg.showlabels,'markers')
  if ischar(cfg.highlight) 
    hp2 = plot(y,x,cfg.emarker,'Color',cfg.ecolor,'markersize',cfg.emarkersize);
  elseif isnumeric(cfg.highlight) 
    normal = setdiff(1:length(X), cfg.highlight);
    hp2    = plot(y(normal),        x(normal),        cfg.emarker,  'Color', cfg.ecolor,  'markersize', cfg.emarkersize);
    hp2    = plot(y(cfg.highlight), x(cfg.highlight), cfg.hlmarker, 'Color', cfg.hlcolor, 'markersize', cfg.hlmarkersize, ...
                                                                                          'linewidth',  cfg.hllinewidth);
  elseif iscell(cfg.highlight)
    hp2 = plot(y,x,cfg.emarker,'Color',cfg.ecolor,'markersize',cfg.emarkersize);
    for iCell = 1:length(cfg.highlight)
    hp2    = plot(y(cfg.highlight{iCell}), x(cfg.highlight{iCell}), cfg.hlmarker{iCell}, 'Color', cfg.hlcolor{iCell},...
                                             'markersize', cfg.hlmarkersize{iCell},'linewidth',  cfg.hllinewidth{iCell});
    end    
  else
    error('Unknown highlight type');
  end;
elseif any(strcmp(cfg.electrodes,{'highlights','highlight'}))
  if isnumeric(cfg.highlight) 
    hp2 = plot(y(cfg.highlight), x(cfg.highlight), cfg.hlmarker, 'Color', cfg.hlcolor, 'markersize', cfg.hlmarkersize, ...
                                                                                       'linewidth',cfg.hllinewidth);
  else
    error('Unknown highlight type');
  end;
elseif strcmp(cfg.electrodes,'labels') || strcmp(cfg.showlabels,'yes') 
  for i = 1:size(labels,1) 
    text(y(i), x(i), labels(i,:), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                                  'Color', cfg.ecolor, 'FontSize', cfg.efsize);
  end
elseif strcmp(cfg.electrodes,'numbers') || strcmp(cfg.showlabels,'numbers') 
  for i = 1:size(labels,1)
    text(y(i), x(i), int2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                                 'Color', cfg.ecolor, 'FontSize',cfg.efsize);
  end
elseif strcmp(cfg.electrodes,'dotnum') 
  for i = 1:size(labels,1)
    text(y(i), x(i), int2str(i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                                 'Color', cfg.ecolor, 'FontSize', cfg.efsize);
  end
  if ischar(cfg.highlight) 
    hp2 = plot(y, x, cfg.emarker, 'Color', cfg.ecolor, 'markersize', cfg.emarkersize);
  elseif isnumeric(cfg.highlight) 
    normal = setdiff(1:length(X), cfg.highlight);
    hp2    = plot(y(normal)       , x(normal),        cfg.emarker,  'Color', cfg.ecolor,  'markersize', cfg.emarkersize);
    hp2    = plot(y(cfg.highlight), x(cfg.highlight), cfg.hlmarker, 'Color', cfg.hlcolor, 'markersize', cfg.hlmarkersize, ...
                                                                    'linewidth', cfg.hllinewidth);
  else
    error('Unknown highlight type');
  end;
end

if strcmp(cfg.outline, 'scalp') || strcmp(cfg.outline, 'miriam')
  % Define the outline of the head, ears and nose:
  l     = 0:2*pi/100:2*pi;
  tip   = rmax*1.15; base = rmax-.004;
  EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
  EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
  % Plot head, ears, and nose:
  plot(cos(l).*rmax, sin(l).*rmax, 'color', cfg.hcolor, 'Linestyle', '-', 'LineWidth', cfg.hlinewidth);
  plot([0.18*rmax;0;-0.18*rmax], [base;tip;base], 'Color', cfg.hcolor, 'LineWidth', cfg.hlinewidth);
  plot( EarX, EarY, 'color', cfg.hcolor, 'LineWidth', cfg.hlinewidth)
  plot(-EarX, EarY, 'color', cfg.hcolor, 'LineWidth', cfg.hlinewidth)
elseif strcmpi(cfg.outline, 'ecog')
  % FIXME drawing an outline of the ECoG electrode grid has not been
  % implemented yet
  warning('drawing an outline of the ECoG electrode grid has not been implemented yet')
end

% Write comment:
if isfield(cfg, 'comment') 
  if strcmp(cfg.commentpos, 'title')
    title(cfg.comment, 'Fontsize', cfg.fontsize);
  else
    text(x_COMNT, y_COMNT, cfg.comment, 'Fontsize', cfg.fontsize, 'HorizontalAlignment', HorAlign, 'VerticalAlignment', VerAlign);
  end
end

% plot colorbar:
if isfield(cfg, 'colorbar') && ~all(data == data(1))
    if strcmp(cfg.colorbar, 'yes')
        colorbar;
    elseif ~strcmp(cfg.colorbar, 'no')
        colorbar(cfg.colorbar);
    end
end

hold off
axis off
if strcmpi(cfg.outline, 'scalp') || strcmp(cfg.outline, 'miriam')
  xlim([-.6 .6]);
  ylim([-.6 .6]);
elseif strcmpi(cfg.outline, 'ecog')
  % do not apply any scaling
end
