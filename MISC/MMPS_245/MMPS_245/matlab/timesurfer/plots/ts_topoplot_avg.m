function ts_topoplot_avg(avg_data,varargin);
% ts_topoplot_avg - plots eeg potentials from avg_data structure on cartoon
%   head using fieldtrips's topoplotER
%
% Usage:
%   ts_topoplot_avg(avg_data,'key1', value1,...);
%
% equired input:
%  avg_data: average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  'condition': condition number (not event code) used to index avg_data.averages
%    {default: 1}
%  'chantype': channel type (e.g. 'mag','grad1','grad2','eeg','gradpow','other')
%    'gradpow' is the hypotenuse of grad1 and grad2 pairs (power)
%    {default: 'eeg'}
%  'time0': start time of averaged data range (msec)
%    {default: 50}
%  'time1': end time of averaged data range (msec)
%    {default: 60}
%  'scale_max': max value for scale (uVolts, fT, or fT/cm)
%    if empty, uses zlim or min/max values
%    {default: []}
%  'zlim': vector of scaling values
%    if empty, uses min/max values
%    {default: []}
%  'badchanfile': name of text file containing bad channel labels
%    {default: []}
%
%  Created:  06/06/06   by Don Hagler
%  Last Mod: 08/03/15   by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'condition',1,[],...
  'chantype','eeg',{'mag','grad1','grad2','eeg','gradpow','grad','other'},...
  'time0',50,[],...
  'time1',60,[],...
  'scale_max',[],[],...
  'zlim',[],[],...
  'badchanfile',[],[],...
});

% convert to field trip format
switch parms.chantype
  case {'gradpow','grad'}
    FT_grad1_data = ts_avg2fieldtrip(avg_data,'condition',parms.condition,...
      'chantype','grad1','badchanfile',parms.badchanfile);
    FT_grad2_data = ts_avg2fieldtrip(avg_data,'condition',parms.condition,...
      'chantype','grad2','badchanfile',parms.badchanfile);
    FT_data = FT_grad1_data;
    FT_data.avg = sqrt(FT_grad1_data.avg.^2 + FT_grad2_data.avg.^2);
  case {'mag' 'grad1' 'grad2' 'eeg' 'other'}
    FT_data = ts_avg2fieldtrip(avg_data,'condition',parms.condition,...
      'chantype',parms.chantype,'badchanfile',parms.badchanfile);
  otherwise
    error('unsupported chantype (%s)',parms.chantype);
end;

if isempty(FT_data)
  error('failed to convert data');
end;

switch parms.chantype
  case {'grad1' 'grad2' 'gradpow' 'grad' 'other'}
    units = 10^-13;
  case 'mag'
    units = 10^-15;
  case 'eeg'
    units = 10^-6;
end; 

% set configuration for topoplot
cfg = [];
cfg.xlim = [parms.time0/1000 parms.time1/1000];
if ~isempty(parms.scale_max)
  cfg.zlim = [-parms.scale_max,parms.scale_max]*units;
elseif ~isempty(parms.zlim)
  cfg.zlim = parms.zlim*units;
else
  cfg.zlim = 'maxmin';
end;
cfg.fontsize = 1;
cfg.interactive = 'yes';
cfg.colorbar = 'no';
%cfg.showindex   = 'yes';
%cfg.baseline = 'yes';
%cfg.baselinetype = 'relative';
cfg.comment = [];

% plot data
ts_topoplotER(cfg,FT_data);

colormap(mmil_cmap_blueblackred);

return;

