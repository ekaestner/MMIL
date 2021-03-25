function ts_plot_avg(avg_data,varargin)
%function ts_plot_avg(avg_data,[options])
%
% Purpose: plot MEG/EEG sensor waveforms for single channel
%
% Usage: ts_plot_avg(avg_data,'key1', value1,...);
%
% Required input:
%  avg_data: average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  'visible_flag': [0|1] whether to display plots (or just save to file)
%     {default = 0}
%  'tif_flag': [0|1] whether to save plot as tif file
%     {default = 0}
%  'eps_flag': [0|1] whether to save plot as eps file (vector graphics)
%     {default = 0}
%  'outstem': output file stem if tif_flag or eps_flag = 1
%     {default = 'avg_data'}
%  'conditions': vector of condition numbers (not event code)
%    used to index avg_data.averages
%    {default = [1]}
%  'condnames': cell array of condition names for legend
%    {default = []}
%  'chantype': channel type (e.g. 'grad','mag','eeg')
%    If supplied, overplot data from all channels of this type
%    Otherwise, will use single channel specified by channame or channum
%    {default = []}
%  'channame': channel name
%    supply either channame or channum (channame takes precedence)
%    {default = []}
%  'channum': channel number
%    {default = 1}
%  'legend_flag': [0|1] toggle display legend
%    {default = 1}
%  'legend_loc': legend location
%    {default = 'EastOutside'}
%  'axis_flag': [0|1] toggle display default matlab axis
%    {default = 1}
%  'vert_lines': vector of time points to make vertical lines
%    {default = [0]}
%  'horz lines': vector of y-values make horizontal lines
%    {default = [0]}
%  'scale_range': min and max values for scale (in uVolts, fT, or fT/cm)
%    {default = [] -> will use max/min}
%  'time_range': min and max time points (in msec)
%    {default = [] -> will use max/min}
%  'label_flag': [0|1] toggle display x and y axis labels and title
%    {default = 1}
%  'title': title string; if not supplied, will use channel name
%    {default = []}
%  'fontsize': font size for label text
%    {default = 12}
%  'fontname': font name for label text
%    {default = 'Helvetica'}
%  'linewidth': trace line width
%    {default = 1.5}
%  'grid_linewidth': line width for horizontal and vertical lines
%    {default = 1}
%  'color_order': cell array containing matlab color characters
%    {default = {'b' 'g' 'r' 'c' 'm' 'y' 'k'}}
%
%  created:       10/11/07   by Don Hagler
%  last modified: 03/12/12   by Don Hagler
%

%% todo: add fig_size option, no title or ylabel if empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

nconds = length(avg_data.averages);
nchans = avg_data.num_sensors;
time = avg_data.averages(1).time*1000;
ntpoints = length(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms( varargin, {...
  'visible_flag',true,[false true],...
  'tif_flag',false,[false true],...
  'eps_flag',false,[false true],...
  'outstem','avg_data',[],...
  'conditions',[1:nconds],[],...
  'condnames',[],[],...
  'chantype',[],[],...
  'channame',[],[],...
  'channum',1,[1 nchans],...
  'scale_range',[],[],...
  'time_range',[],[],...
  'linewidth',1.5,[0.1 100],...
  'grid_linewidth',1,[0.1 100],...
  'fontsize',12,[1 100],...
  'fontname','Helvetica',[],...
  'label_flag',true,[false true],...
  'title',[],[],...
  'legend_flag',true,[false true],...
  'legend_loc','EastOutside',[],...
  'axis_flag',true,[false true],...
  'vert_lines',0,[],...
  'horz_lines',0,[],...
  'color_order',{'b' 'g' 'r' 'c' 'm' 'y' 'k'},[],...
  'bgcolor','w',[]...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if isempty(parms.conditions)
  parms.conditions = [1:nconds];
end;
if ~isempty(parms.condnames)
  if ~iscell(parms.condnames), parms.condnames = {parms.condnames}; end;
  if length(parms.condnames) ~= length(parms.conditions)
    error('number of condnames (%d) does not match number of conditions (%d)',...
      length(parms.condnames),length(parms.conditions));
  end;
end;

if isempty(parms.color_order)
  error('color_order is empty');
end;
if ~iscell(parms.color_order), parms.color_order = {parms.color_order}; end;
for i=1:length(parms.color_order)
  if ~ischar(parms.color_order{i})
    error('color_order must contain color strings (e.g. ''b'', ''b*'')')
  end;
end;

if ~isempty(parms.chantype)
  parms.multichan_flag = 1;
  chantypes = {avg_data.sensor_info.typestring};
  parms.channum = ...
    find(strncmp(lower(parms.chantype),chantypes,length(parms.chantype)));
else
  parms.multichan_flag = 0;
  if ~isempty(parms.channame)
    parms.channum = find(strcmp(parms.channame,{avg_data.sensor_info.label}));
    if isempty(parms.channum)
      error('channame %s not found in avg_data.sensor_info.label',...
        parms.channame);
    end;
    chan_info = avg_data.sensor_info(parms.channum);
  else
    chan_info = avg_data.sensor_info(parms.channum);
    parms.channame = chan_info.label;
  end;
  parms.chantype = chan_info.typestring;
end;

switch parms.chantype
case {'grad1' 'grad2','grad'}
  scalefact = 10^13;
  units = 'fT/m';
case 'mag'
  scalefact = 10^15;
  units = 'fT';
case 'eeg'
  scalefact = 10^6;
  units = 'uV';
otherwise
  scalefact = 1;
  units = 'unknown';
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data

cla;
whitebg(parms.bgcolor);
hold on;

data = zeros(length(parms.conditions),ntpoints);
j=1;
legstrs = [];
for i=1:length(parms.conditions)
  if (j>length(parms.color_order)) j=1; end;
  color = parms.color_order{j}; j=j+1;
  c = parms.conditions(i);
  if c<1 | c>nconds
    error('bad condition number %d -- avg_data has %d conditions',...
      c,nconds);
  end;
  tmp_data = squeeze(avg_data.averages(c).data(parms.channum,:))*scalefact;
  if parms.multichan_flag
    plot(time,tmp_data,'LineWidth',parms.linewidth);
  else
    plot(time,tmp_data,color,'LineWidth',parms.linewidth);
  end;
  if parms.legend_flag
    if isempty(parms.condnames)
      legstrs{i} = sprintf('%d: Event Code = %d',...
        c,avg_data.averages(c).event_code);
    else
      legstrs{i} = parms.condnames{i};
    end;
  end;
end;

axis tight;
if ~isempty(parms.scale_range)
  set(gca,'YLim',parms.scale_range);
end;
if ~isempty(parms.time_range)
  set(gca,'XLim',parms.time_range);
end;
set(gca,'FontSize',parms.fontsize,'FontName',parms.fontname);
if parms.label_flag
  xlabel('Time (ms)','FontSize',parms.fontsize,'FontName',parms.fontname);
  ylabel(units,'FontSize',parms.fontsize,'FontName',parms.fontname);
  title_opts = {'FontSize',parms.fontsize,'FontName',parms.fontname,'Color','k'};
  if ~isempty(parms.title)
    title(parms.title,title_opts{:});
  elseif parms.multichan_flag
    title(['all ' parms.chantype ' channels'],title_opts{:})
  else
    title(parms.channame,title_opts{:})
  end;
end;
if parms.legend_flag
  legend(legstrs,'Location',parms.legend_loc);
end;
if ~parms.axis_flag, axis off; end;

for i=1:length(parms.vert_lines)
  x = parms.vert_lines(i);
  plot([x,x],get(gca,'YLim'),'k','LineWidth',parms.grid_linewidth);
end;
for i=1:length(parms.horz_lines)
  y = parms.horz_lines(i);
  plot(get(gca,'XLim'),[y,y],'k','LineWidth',parms.grid_linewidth);
end;

if ~parms.visible_flag
  set(gcf,'Visible','Off');
end;

if parms.tif_flag
  print('-dtiff',[parms.outstem '.tif']);
end;

if parms.eps_flag
  mmil_printeps(gcf,parms.outstem);
end;

if ~parms.visible_flag
  close(gcf);
end;

return;

