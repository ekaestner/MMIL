function topo(data,varargin)
% zlim: maxabs, maxmin, or [zmin zmax]
cfg     = mmil_args2parms( varargin, ...
                         { 'emarker'      ,'.',[],...
                           'ecolor'       ,[0 0 0],[],...
                           'emarkersize'  ,2,[],...
                           'electrodes'   ,[],[],...
                           'highlight'    ,[],[],...
                           'hlmarker'     ,'o',[],...
                           'hlcolor'      ,[0 0 0],[],...
                           'hlmarkersize' ,6,[],...
                           'hllinewidth'  ,3,[],...
                           'hcolor'       ,[0 0 0],[],...
                           'hlinewidth'   ,2,[],...
                           'fontsize'     ,6,[],...
                           'efsize'       ,get(0,'DefaultAxesFontSize'),[],...
                           'method'       ,'stereographic',[],...
                           'interplimits' ,'head',[],...
                           'interp'       ,'v4',[],...   
                           'gridscale'    ,100,[],...
                           'fig'          ,[],[],...
                           'nrows'        ,1,[],...
                           'ncols'        ,1,[],...
                           'xlim'         ,[],[],...
                           'toilim'       ,[],[],...
                           'chantype'     ,[],[],...
                           'layout'       ,[],[],...
                           'shading'      ,'flat',{'flat','interp'},...
                           'numcontours'  ,6,[],...
                           'showlabels'   ,0,[],...
                           'badchans'     ,[],[],...
                           'badlabels'    ,[],[],...
                           'iEEG_flag'    ,0,{0,1},...
                           'title',[],[],...
                           'avgovertime'  ,'yes',[],...
                           'zlim',[],[],...
                         }, ...
                         false );
cfg.nframes = cfg.nrows * cfg.ncols;
if ~isempty(cfg.xlim)   , cfg.toilim   = cfg.xlim; end
if isempty(cfg.chantype), cfg.chantype = data.sensor_info(1).typestring; end
if ischar(cfg.badchans) , cfg.badchans = strmatch(cfg.badchans,{data.sensor_info.label}); end
if strcmp(cfg.interplimits,'grid'), cfg.iEEG_flag = 1; end

data = ts_data_selection(data,'toilim',cfg.toilim,'chantype',cfg.chantype,'verbose',0,'badchans',cfg.badchans,'badlabels',cfg.badlabels);
[datatype,datafield] = ts_object_info(data);

% prepare position Info
nchan   = data.num_sensors;
pos     = zeros(nchan,3); % (x,y,z) for each channel
if isempty(cfg.layout)
  T     = data.coor_trans.device2head;
  for k = 1:nchan
    loc         = data.sensor_info(k).loc;
    if any(strmatch('grad',data.sensor_info(k).typestring)) || ...
       any(strmatch('mag' ,data.sensor_info(k).typestring))
      loc       = T*loc;
    end  
    pos(k,1:3)  = loc(1:3,4);
  end
  method = cfg.method; % gnomic, stereographic, ortographic, inverse, polar
  prj    = elproj(pos, method); % * [0 1; -1 0];
            % ELPROJ makes a azimuthal projection of a 3D electrode cloud
            %  on a plane tangent to the sphere fitted through the electrodes
            %  the projection is along the z-axis
  X = prj(:,1);   % x-coordinates
  Y = prj(:,2);   % y-coordinates  
elseif ischar(cfg.layout) && exist(cfg.layout,'file')
  [chNum,X,Y,Width,Height,Lbl,Rem] = textread(cfg.layout,'%f %f %f %f %f %q %q');
  for i=1:length(Lbl)
    if ~isempty(Rem{i})
      % this ensures that channel names with a space in them are also supported (i.e. Neuromag)
      Lbl{i} = [Lbl{i} ' ' Rem{i}];
    end
  end
%   [sel,jnk] = match_str(Lbl,{data.sensor_info.label});
%   X = X(sel);
%   Y = Y(sel);
  [sel1,sel2] = match_str({data.sensor_info.label},Lbl);
  X = X(sel2);
  Y = Y(sel2);
  if any(X(:)<0), X = X - min(X(:)); end
  if any(Y(:)<0), Y = Y - min(Y(:)); end
  data.(datafield).data  = data.(datafield).data(sel1,:);
  data.sensor_info       = data.sensor_info(sel1);
	try
  	data.(datafield).stdev = data.(datafield).stdev(sel1,:);
	catch
		data.(datafield).stdev = data.(datafield).data;
	end
%   lay.pos    = [X Y];
%   lay.width  = Width;
%   lay.height = Height;
%   lay.label  = Lbl;
elseif isstruct(cfg.layout)
  X = cfg.layout.pos(:,1);
  Y = cfg.layout.pos(:,2);
else
  error('Layout not found.');
end

cfg.normal = setdiff(1:length(X),cfg.highlight);
% Scale the data to a circle with x-axis and y-axis: -0.45 to 0.45
y = 0.9*((X-min(X))/(max(X)-min(X))-0.5); % NOTE: x becomes y and y becomes x, in griddata is also reversed which makes x x and y y again
x = 0.9*((Y-min(Y))/(max(Y)-min(Y))-0.5);
interplimits = cfg.interplimits;
  % 'electrodes' to furthest electrode
  % 'head' to edge of head
% Find limits for interpolation:
if strcmp(interplimits,'head') || strcmp(interplimits,'grid')
  xmin = min(-.5,min(x)); xmax = max(0.5,max(x));
  ymin = min(-.5,min(y)); ymax = max(0.5,max(y));
else
  xmin = max(-.5,min(x)); xmax = min(0.5,max(x));
  ymin = max(-.5,min(y)); ymax = min(0.5,max(y));
end
if ischar(cfg.zlim)
  if strcmp(cfg.zlim,'absmax')
    zmax     = max(abs(data.(datafield).data(:)));
    zmin     = -zmax;
  elseif strcmp(cfg.zlim,'maxmin')
    zmin     = min(data.(datafield).data(:));
    zmax     = max(data.(datafield).data(:));
  end
else
  zmin     = min(data.(datafield).data(:));
  zmax     = max(data.(datafield).data(:));  
end
gridscale  = cfg.gridscale;                   % resolution
interp     = cfg.interp;                      % 'linear','cubic','nearest','v4'
xi         = linspace(xmin,xmax,gridscale);   % x-axis description (row vector)
yi         = linspace(ymin,ymax,gridscale);   % y-axis description (row vector)
delta      = xi(2)-xi(1);   
cfg.labels = {data.sensor_info.label};

% Prepare data and create topoplot series
if isempty(cfg.fig)
  screensize = get(0,'ScreenSize');
  figure('NumberTitle','on','Color','w','Position',[1 1 .9*screensize(3) .9*screensize(4)]);
  hold on
elseif ~isequal(cfg.fig,0)
  figure(cfg.fig);
end
N = cfg.nframes;              % number frames
t = data.(datafield).time;    % time vector
n = floor(length(t)/N);       % number samples per frame
k = 1;                        % frame index
if N > 1
  clf
end
set(gcf,'Name',sprintf('%s',cfg.chantype));
for r = 1:cfg.nrows           % row index
  for c = 1:cfg.ncols         % column index
    if N > 1
      subplot(cfg.nrows,cfg.ncols,k)
    end
    hold on
    samp = [1:n] + (k-1)*n;   % samples for this frame
    if strcmp(cfg.avgovertime,'yes')
      dat = mean(data.(datafield).data(:,samp),2);
    else
      dat = data.(datafield).data(:,samp(1));
    end
    tt   = t(samp);
    tlim = [tt(1) tt(end)];
    cfg.comment = sprintf('t=[%3.3g %3.3g]',tlim);
    % Interpolate data; NOTE: undo the reversal of x & y
    [Xi,Yi,Zi] = griddata(y,x,dat,yi',xi,interp);
      % griddata uses meshgrid to create evenly-spaced XI & YI
    % draw the topoplot for this interval
    draw_topoplot(Xi,Yi,Zi,x,y,delta,cfg);
    hold off
    axis tight;%square
    k = k + 1;
    if isnumeric(cfg.zlim) && ~isempty(cfg.zlim)
      set(gca,'clim',cfg.zlim);
    else
      set(gca,'clim',[zmin zmax]);
    end
  end
end   
if ischar(cfg.title)
  cfg.title = strrep(cfg.title,'_','\_');
  annotation('textbox','Position',[0 .9 1 .1],'VerticalAlignment','middle','HorizontalAlignment','center',...
             'Color','k','FontSize',10,'FontWeight','bold','FitHeightToText','on','LineStyle','none','String',cfg.title);
end

function draw_topoplot(Xi,Yi,Zi,x,y,delta,cfg)
if ~cfg.iEEG_flag
  % Take data within head
  rmax       = .5;
  mask       = (sqrt(Xi.^2+Yi.^2) <= rmax);
  Zi(mask==0)= NaN;
end

% Draw topoplot on head
numcontours = cfg.numcontours;
shading     = cfg.shading;
% contour(Xi,Yi,Zi,numcontours,'k');
h = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',shading);

% calculate colormap limits
zmin = min(abs(Zi(:)));
zmax = max(abs(Zi(:)));
caxis([zmin zmax])

% draw head
if ~cfg.iEEG_flag
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
end

  % draw electrodes
  if isequal(cfg.electrodes,'yes') || isequal(cfg.electrodes,1)
    hp2 = plot(y(cfg.normal),    x(cfg.normal),    cfg.emarker,  'Color', cfg.ecolor,  'markersize', cfg.emarkersize);
  end
  if ~isempty(cfg.highlight)
    hp2 = plot(y(cfg.highlight), x(cfg.highlight), cfg.hlmarker, 'Color', cfg.hlcolor, 'markersize', cfg.hlmarkersize,'linewidth', cfg.hllinewidth);
  end
  if isequal(cfg.showlabels,'yes') || isequal(cfg.showlabels,1)
    % add labels
    for ch = 1:length(cfg.labels)
      text(y(ch), x(ch), cfg.labels{ch}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                                    'Color', cfg.ecolor, 'FontSize', cfg.fontsize);%cfg.efsize);
    end
  end
  x_COMNT =  0.6; 
  y_COMNT = -0.6;
  HorAlign = 'right';
  VerAlign = 'bottom';
% Write comment:
% if isfield(cfg, 'comment') 
%   if strcmp(cfg.commentpos, 'title')
%     title(cfg.comment, 'Fontsize', cfg.fontsize);
%   else
    text(x_COMNT, y_COMNT, cfg.comment, 'Fontsize', cfg.fontsize, 'HorizontalAlignment', HorAlign, 'VerticalAlignment', VerAlign);
%   end
% end

% colorbar
% hold off
axis off
xlim([-.6 .6]);
ylim([-.6 .6]);
% 
% 
%   if strcmp(cfg.style,'contour')
%     contour(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
%   elseif strcmp(cfg.style,'both')
%     contour(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
%     h = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',cfg.shading);
%     if exist('maskZ','var'),
%       set(h, 'AlphaData', maskZ);
%       alim([0 1]);
%       set(h, 'FaceAlpha', 'interp');
%     end
%   elseif strcmp(cfg.style,'straight')
%     h = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',cfg.shading);
%     if exist('maskZ','var'),
%       set(h, 'AlphaData', maskZ);
%       alim([0 1]);
%       set(h, 'FaceAlpha', 'interp');
%     end
%   elseif strcmp(cfg.style,'fill')
%     contourf(Xi,Yi,Zi,cfg.contournum,cfg.contcolor);
%   else
%     error('Invalid style')
%   end
