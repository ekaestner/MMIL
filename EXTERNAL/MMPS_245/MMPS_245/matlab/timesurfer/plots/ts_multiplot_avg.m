function ts_multiplot_avg(avg_data,varargin);
% ts_multiplot_avg - plots MEG/EEG sensor waveforms
%
% Usage:
%   ts_multiplot_avg(avg_data,'key1', value1,...);
%
% equired input:
%  avg_data - average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  'conditions': vector of condition numbers (not event code)
%    used to index avg_data.averages
%    {default = 1}
%  'chantype': channel type (e.g. 'mag','grad1','grad2','eeg','grad','gradpow')
%    'grad' is both grad1 and grad2 channels
%    'gradpow' is the hypotenuse of grad1 and grad2 pairs (power)
%    {default = 'eeg'}
%  'scale_max': max value for scale (uVolts, fT, or fT/cm)
%    {default = [] -> will use max/min}
%  'linewidth': trace line width
%    {default = 1.5}
%  'badchanfile': name of text file containing bad channel labels
%    {default = []}
%  'label_flag': [0|1] whether to display labels
%    {default = 1}
%
%  created:  06/06/06   by Don Hagler
%  last mod: 08/11/11   by Don Hagler
%

%% todo: gradpow does not work if there are bad channels (mismatch)
%% todo: make a way to show a subset of channels

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'conditions',1,[],...
  'chantype','eeg',[],...
  'scale_max',[],[],...
  'linewidth',1.5,[],...
  'badchanfile',[],[],...
  'label_flag',true,[false true],...
});

FT_data = [];
for c=1:length(parms.conditions)
  % convert to field trip format
  switch parms.chantype
  case {'mag' 'grad1' 'grad2' 'eeg' 'gradpow' 'grad'},;
  otherwise
    help(mfilename);
    fprintf('\n%s: unsupported chantype (%s)\n',mfilename);
    return;
  end;
  if strcmp(parms.chantype,'gradpow')
    FT_grad1_data = ts_avg2fieldtrip(avg_data,'condition',parms.conditions(c),...
      'chantype','grad1','badchanfile',parms.badchanfile);
    FT_grad2_data = ts_avg2fieldtrip(avg_data,'condition',parms.conditions(c),...
      'chantype','grad2','badchanfile',parms.badchanfile);
    FT_data{c} = FT_grad1_data;
    FT_data{c}.avg = sqrt(FT_grad1_data.avg.^2 + FT_grad2_data.avg.^2);
  else
    FT_data{c} = ts_avg2fieldtrip(avg_data,'condition',parms.conditions(c),...
      'chantype',parms.chantype,'badchanfile',parms.badchanfile);
  end;
end;

switch parms.chantype
case {'grad1' 'grad2' 'gradpow' 'grad'}
  units = 10^-13;
case 'mag'
  units = 10^-15;
case 'eeg'
  units = 10^-6;
end; 

% set configuration for multiplot
cfg = [];
if ~isempty(parms.scale_max)
  cfg.ylim = [-parms.scale_max,parms.scale_max]*units;
end;
cfg.colorbar = 'no';
if parms.label_flag
  cfg.showlabels = 'yes';
else
  cfg.showlabels = 'no';
end;
cfg.linewidth = parms.linewidth;

% plot data
ts_multiplotER(cfg,FT_data{:});

return;

