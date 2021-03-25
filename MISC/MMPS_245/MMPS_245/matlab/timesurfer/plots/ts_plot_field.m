function ts_plot_field(avg_data,varargin)
%function ts_plot_field(avg_data,[options])
%
% Purpose: plot MEG field on the Elekta/Neuromag Helmet Surface
%
% Required input:
%  avg_data: average data structure
%
% Optional parameters:
% 'condition': condition number (index to avg_data.averages)
%   {default: 1}
% 'time0': start time of averaged data range (msec)
%    {default = 60}
% 'time1': end time of averaged data range (msec)
%    {default = 70}
% 'badchanfile': name of text file containing bad channel labels
%    {default = []}
% 'usemags': [0|1] whether to use magnetometer data
%    {default = 0}
% 'view_angle': 2x1 vector containing azimuth and elevation angles
%    {default = [0 20]}
% 'scale_max': max value for color scale (fT ?)
%    {default = 100}
% 'colormap': name of colormap
%    {default = 'mmil_cmap_blueblackred'}
% 'vv_field_plot: structure loaded from vectorview_field_plot.mat
%     if empty, will load vvfpfile
%     {default = []}
% 'vvfpfile': file name for presaved vectorview_field_plot mat file
%    {default = 'vectorview_field_plot.mat'}
%
%  Created   05/17/06   by Don Hagler
%  Last Mod: 09/27/12   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'condition',1,[],...
  'time0',60,[],...
  'time1',70,[],...
  'badchanfile',[],[],...
  'usemags',false,[false true],...
  'view_angle',[0 20],[],...
  'scale_max',100,[],...
  'colormap','mmil_cmap_blueblackred',[],...
  'vv_field_plot',[],[],...
  'vvfpfile','vectorview_field_plot.mat',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.vv_field_plot)
  load(parms.vvfpfile);
  parms.vv_field_plot = vv_field_plot;
end;

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

cond = parms.condition;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = mean(avg_data.averages(cond).data(ind_meg,t0:t1),2);
chans = zeros(length(ind_meg),1);
chans(ind_chans) = 1;
ts_plot_vv_field(B,'vv_field_plot',parms.vv_field_plot,'chans',chans,...
  'crange',color_range,'view_angle',parms.view_angle,...
  'colormap',parms.colormap);

return;

