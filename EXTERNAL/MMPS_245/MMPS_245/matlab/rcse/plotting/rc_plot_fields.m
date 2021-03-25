function rc_plot_fields(avg_data,varargin)
%function rc_plot_fields(avg_data,[options])
%
% Purpose:plot MEG field on the Elekta/Neuromag Helmet Surface
%
% Required input:
%  avg_data: average data structure
%
% Optional parameters:
% 'conditions': vector of condition numbers (indices to avg_data.averages)
%   if empty, will plot all conditions
%   {default = []}
% 'cond_info': struct array of condition information
%   if empty, will place plots in a regular grid
%   {default = []}
% 'rmax': maximum eccentricity for stimuli in cond_info
%   {default = 10}
% 'time0': start time of averaged data range (msec)
%    {default = 60}
% 'time1': end time of averaged data range (msec)
%    {default = 70}
% 'badchanfile': name of text file containing bad channel labels
%    {default = []}
% 'usemags': [0|1] whether to use magnetometer data
%    {default = 0}
% 'plotsize': size of single head plot (relative to figure)
%    {default = 0.13}
% 'view_angle': 2x1 vector containing azimuth and elevation angles
%    {default = [0 20]}
% 'scale_max': max value for color scale (fT ?)
%    {default = 100}
% 'colormap': name of colormap
%    {default = 'mmil_cmap_blueblackred'}
% 'vvfpfile': file name for presaved vectorview_field_plot mat file
%    {default = 'vectorview_field_plot.mat'}
%
%  Created   05/17/06   by Don Hagler
%  Last Mod: 09/27/12   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'conditions',[],[],...
  'cond_info',[],[],...
  'rmax',10,[],...
  'time0',60,[],...
  'time1',70,[],...
  'badchanfile',[],[],...
  'usemags',false,[false true],...
  'plotsize',[0.13],[],...
  'view_angle',[0 20],[],...
  'scale_max',100,[],...
  'colormap','mmil_cmap_blueblackred',[],...
  'vvfpfile','vectorview_field_plot.mat',[],...
...
  'plotsize_init',0.01,[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(parms.vvfpfile);

% set start and end samples
sfreq = avg_data.sfreq;
t_trigger = avg_data.averages(1).time(1)*1000;
t0 = round((parms.time0 - t_trigger)*sfreq/1000);
t1 = round((parms.time1 - t_trigger)*sfreq/1000);

% MEG chan indices
chan_types = {avg_data.sensor_info.typestring};
ind_meg = find(ismember(chan_types,{'grad1','grad2','mag'}));
if parms.usemags
  ind_chans = ind_meg;
else
  ind_chans = find(ismember(chan_types,{'grad1','grad2'}));
end;

% exclude badchans
ind_badchans_avg = find(cell2mat({avg_data.sensor_info.badchan})==1);
% read badchan file
labels = {avg_data.sensor_info.label}; 
if ~isempty(parms.badchanfile)
  ind_badchans = ts_read_txt_badchans(parms.badchanfile,labels);
else
  ind_badchans = [];
end;
ind_badchans = unique([ind_badchans_avg,ind_badchans]);
ind_chans = setdiff(ind_chans,ind_badchans); % exclude bad chans

color_range = [-parms.scale_max,parms.scale_max];

if isempty(parms.conditions)
  parms.conditions = [1:length(avg_data.averages)];
end;

nconds = length(parms.conditions);

if isempty(parms.cond_info)
  nrows = ceil(sqrt(nconds));
  ncols = ceil(nconds/nrows);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create plots with initially small sizes (to avoid erasing old plots)
hold on;
for c=1:nconds
  cond = parms.conditions(c);
  B = mean(avg_data.averages(cond).data(ind_meg,t0:t1),2);
  chans = zeros(length(ind_meg),1);
  chans(ind_chans) = 1;
  if isempty(parms.cond_info)
    y = (ceil(c/nconds)-0.5)/nrows;
    x = (rem(c,ncols)-0.5)/ncols;
  else
    r = parms.cond_info(cond).ecc/(2*parms.rmax);
    th = parms.cond_info(cond).theta*pi/180;
    x = 0.5 + r*cos(th);
    y = 0.5 + r*sin(th);
  end;
  plotloc(1) = x - parms.plotsize_init/2;
  plotloc(2) = y - parms.plotsize_init/2;
  plotloc(3) = parms.plotsize_init;
  plotloc(4) = parms.plotsize_init;  
  subplot('position',plotloc);
  ts_plot_vv_field(B,'vv_field_plot',vv_field_plot,'chans',chans,...
    'crange',color_range,'view_angle',parms.view_angle,...
    'colormap',parms.colormap);
  fp_axis(c) = gca;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resize plots
for c=1:length(parms.conditions)
  cond = parms.conditions(c);
  if isempty(parms.cond_info)
    y = (ceil(c/nconds)-0.5)/nrows;
    x = (rem(c,ncols)-0.5)/ncols;
  else
    r = parms.cond_info(cond).ecc/(2*parms.rmax);
    th = parms.cond_info(cond).theta*pi/180;
    x = 0.5 + r*cos(th);
    y = 0.5 + r*sin(th);
  end;
  plotloc(1) = x - parms.plotsize/2;
  plotloc(2) = y - parms.plotsize/2;
  plotloc(3) = parms.plotsize;
  plotloc(4) = parms.plotsize;
  set(fp_axis(c),'position',plotloc);
end;

return;

