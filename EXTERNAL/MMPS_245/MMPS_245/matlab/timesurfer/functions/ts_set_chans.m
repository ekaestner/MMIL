function [ind_chans,ind_grad,ind_mag,ind_EEG,ind_bad] = ts_set_chans(avg_data,varargin)
%function [ind_chans,ind_grad,ind_mag,ind_EEG,ind_bad] = ts_set_chans(avg_data,[options])
%
% Purpose: initialize avg_data structure for fit and err
%
% Required Input:
%   avg_data: average data structure
%
% Optional Input:
%  'badchans': vector of bad channel indices
%    {default = []}
%  'badchanfile': name of text file containing bad channel labels
%    {default = []}
%  'usegrad_flag': [0|1] use gradiometer data if available
%    {default = 1}
%  'usemag_flag': [0|1] use magnetometer data if available
%    {default = 1}
%  'useEEG_flag': [0|1] use EEG data if available
%    {default = 1}
%  'channames': cell array of channel names
%    if supplied, ignores usegrad_flag, usemag_flag, and useEEG_flag
%    {default = []}
%  'verbose': [0|1] display warnings if channel types are missing
%    {default = 1}
%
% Output:
%   ind_chans: vector of valid channel indices
%   ind_grad: valid gradiometer indices
%   ind_mag: valid magnetometer indices
%   ind_EEG: valid EEG indices
%   ind_bad: invalid channel indices
%
% Created:  09/17/13 by Don Hagler
% Last Mod: 03/06/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_chans = []; ind_grad = []; ind_mag = []; ind_EEG = []; ind_bad = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'badchans',[],[],...
  'badchanfile',[],[],...
  'usegrad_flag',true,[false true],...
  'usemag_flag',true,[false true],...
  'useEEG_flag',true,[false true],...
  'channames',[],[],...
  'verbose',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_channames = {avg_data.sensor_info.label};

% available channels
ind_grad = ...
  find(strncmp('grad',lower({avg_data.sensor_info.typestring}),length('grad')));
ind_mag  = ...
  find(strcmp('mag',lower({avg_data.sensor_info.typestring})));
ind_EEG = ...
  find(strcmp('eeg',lower({avg_data.sensor_info.typestring})));

if ~isempty(parms.channames)
  [~,ind_chans] = intersect(all_channames,parms.channames);
else
  if isempty(ind_grad)
    if parms.usegrad_flag && parms.verbose
      fprintf('%s: WARNING: no gradiometers found, setting usegrad_flag = 0\n',...
        mfilename);
    end;
    parms.usegrad_flag = 0;
  end;
  if isempty(ind_mag)
    if parms.usemag_flag && parms.verbose
      fprintf('%s: WARNING: no magnetometers found, setting usemag_flag = 0\n',...
        mfilename);
    end;
    parms.usemag_flag = 0;
  end;
  if isempty(ind_EEG)
    if parms.useEEG_flag && parms.verbose
      fprintf('%s: WARNING: no EEG channels found, setting useEEG_flag = 0\n',...
        mfilename);
    end;
    parms.useEEG_flag = 0;
  end;
  % valid chans
  if parms.usegrad_flag
    ind_chans = union(ind_chans,ind_grad);
  end;
  if parms.usemag_flag
    ind_chans = union(ind_chans,ind_mag);
  end;
  if parms.useEEG_flag
    ind_chans = union(ind_chans,ind_EEG);
  end;
end;

if isempty(ind_chans) && parms.verbose
  fprintf('%s: WARNING: no valid channels selected\n',mfilename);
end;

% bad chans
ind_bad = find(cell2mat({avg_data.sensor_info.badchan})==1);
ind_bad = union(parms.badchans,ind_bad);
sensor_labels = {avg_data.sensor_info.label}; 
if ~isempty(parms.badchanfile)
  ind_tmp = ts_read_txt_badchans(parms.badchanfile,sensor_labels);
  ind_bad = union(ind_bad,ind_tmp);
end;

% exclude bad chans
ind_chans = setdiff(ind_chans,ind_bad);

ind_grad = intersect(ind_grad,ind_chans);
ind_mag = intersect(ind_mag,ind_chans);
ind_EEG = intersect(ind_EEG,ind_chans);

