function [FT_data] = ts_data2fieldtrip(data,varargin);
% ts_data2fieldtrip - converts a {avg,epoch,timefreq}_data structure created by avg_fif_data
%                    into a FieldTrip FT_data structure, which subsequently
%                    can be used for topoplotER and other FieldTrip tools
%
% Usage:
%  [FT_data] = ts_data2fieldtrip(data,'condition',condition...);
%
% Required input:
%  data - avg_data, epoch_data, or timefreq_data structure (see avg_fif_data or ts_timefreq)
%
% Optional parameters:
%  condition - condition number (not event code) used to index data inside structure
%    {default: 1}
%  time0 - start time of data range (msec)
%    {default: first time point}
%  time1 - end time of data range (msec)
%    {default: last time point}
%  badchanfile - name of text file containing bad channel labels
%    {default: []}
%  badchans    - indices of bad channels
%
%  chantype - channel type (e.g. 'mag','grad1','grad2','eeg', or 'grad')
%    {default: 'all'} ('all' returns all non-"other" channels)
%  channels - vector of channel indices - overrides chantype
%    {default: []}
%  dimord - epoch organization -- either 'trial_chan_time' or 'chan_time'
%    {default: 'trial_chan_time'}
%
% Output:
%   FT_data - FIELDTRIP structure
%
% based on FieldTrip's eeglab2fieldtrip by Robert Oostenveld
% See http://www.ru.nl/fcdonders/fieldtrip/
%
% See also: ts_avg_fif_data, ts_timefreq
%
%  created:       04/20/06  by Don Hagler
%  last modified: 04/27/11  by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 12/01/07 - Corrected permutation of epochs to trial_chan_time

%  Revision 07/29/08 by Jason Sherfey
%  Added conversion for timefreq data with dimord='rpt_chan_freq_time'

% 03/21/09 - added support for cross-spectra
% 06/19/09 - removed factor of two in power calculation 

%  Revision 04/25/10 by Jason Sherfey
%  Add recursive calling for multiple conditions

% 04/27/11 - changed isint to mmil_isint, changed indentation (DH)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_check_nargs(nargin, 1);

parms = mmil_args2parms( varargin, ...
                         { 'channels', [], [], ...
                           'chantype', 'all', {'all', 'mag' 'grad1' 'grad2' 'eeg', 'ieeg', 'other', 'grad', 'meg'}, ...
                           'badchanfile', [], [], ....
                           'badchans',[],[],...
                           'time0', [], [], ...
                           'time1', [], [], ...
                           'condition', 1, [], ...
                           'dimord', [], {'trial_chan_time', 'chan_time', 'chan_freq_time', 'rpt_chan_freq_time'} ...
                         }, ...
                         true );

%DEFAULT_AVG_DIMORD = 'trial_chan_time';
%DEFAULT_EPOCH_DIMORD = 'chan_time';
%EPOCH_DATA = 1;
%AVERAGE_DATA = 2;
%TIMEFREQ_DATA = 3;
%FT_data = [];
  
% DATA_FIELD = ts_objecttype(data);
[datatype DATA_FIELD dataparam] = ts_object_info(data);

if ischar(parms.condition) && strcmp(parms.condition,'all')
  parms.condition = 1:length(data.(DATA_FIELD));
end

% call recursively if multiple conditions
if length(parms.condition) > 1
  allconds = parms.condition;
  FT_data  = {};
  for k = 1:length(allconds)
    parms.condition = allconds(k);
    args            = mmil_parms2args(parms);
    FT_data{k}      = ts_data2fieldtrip(data,args{:});
  end
  return
end

if strcmp(datatype,'timefreq_data') && strcmp(DATA_FIELD,'averages')
  data.timefreq = data.averages;
  data = rmfield(data,'averages');
  DATA_FIELD = 'timefreq';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set / validate dimord
if (isempty(parms.dimord))
  switch (ts_objecttype(data))
    case {'average','averages'}
      parms.dimord  = 'chan_time';
    case {'epoch','epochs'}
      parms.dimord  = 'trial_chan_time';
    case 'timefreq'
				if 		 isfield(data.timefreq(1),'power'), tfdatdims = ndims(data.timefreq(1).power);
				elseif isfield(data.timefreq(1),'cmplx'), tfdatdims = ndims(data.timefreq(1).cmplx);
				else 																			tfdatdims = 3;
				end
        if tfdatdims==4,
            parms.dimord = 'rpt_chan_freq_time';
        else   
            parms.dimord  = 'chan_freq_time';
        end;
    case {'stat','stats'}
        parms.dimord = data.stats(1).dimord;
  end;

else
  switch parms.dimord
    case 'trial_chan_time'
      if (strcmp(ts_objecttype(data), 'averages'))
        mmil_error(parms, 'dimord trial_chan_type cannot be used on averaged data.');
			elseif (strcmp(ts_objecttype(data), 'average'))
        mmil_error(parms, 'dimord trial_chan_type cannot be used on averaged data.');
      end;

    case 'chan_freq_time'
      if (~strcmp(ts_objecttype(data), 'timefreq'))
        mmil_error(parms, 'dimord chan_freq_time can only be used on timefreq data.');
      end;

    case 'rpt_chan_freq_time'
      if (~strcmp(ts_objecttype(data), 'timefreq'))
        mmil_error(parms, 'dimord rpt_chan_freq_time can only be used on timefreq data.');
      end;
  end;
end;

% validate condition
numconds = length(getfield(data,DATA_FIELD));
if length(parms.condition)~=1 | ~isnumeric(parms.condition) | ~mmil_isint(parms.condition) | parms.condition<1
  mmil_error(parms, 'condition must be a single integer >= 1');

elseif (parms.condition > numconds)
  mmil_error(parms, 'condition must be <= number of conditions (%d)',numconds);
end;

if isempty(parms.channels) 
  if isempty(parms.chantype)
    mmil_error(parms, 'Must specify channels or chantype.');
  elseif ~isstr(parms.chantype)
    mmil_error(parms, 'chantype must be a string');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defaults & Data prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate bad channels
if ~isempty(parms.badchanfile)
  badchan_i = ts_read_txt_badchans(parms.badchanfile,{data.sensor_info.label});
else
  try badchan_i = find([data.sensor_info.badchan]); 
	catch badchan_i = []; end
end;
badchan_i = sort(union(badchan_i,parms.badchans));

% Figure out channels from channel types
if isempty(parms.channels)
	if isfield(data,'sensor_info')
	  switch parms.chantype
	    case {'mag' 'grad1' 'grad2' 'eeg' 'ieeg' 'other'}
	      chans = find(strcmp(parms.chantype,{data.sensor_info.typestring}));
	    case {'grad'}
	      chans = find(strncmp(parms.chantype,{data.sensor_info.typestring},...
	        length(parms.chantype)));
	    case 'meg'
	      [a,chans] = find(ismember({data.sensor_info.typestring}, ...
	                    {'mag', 'grad1', 'grad2'}));
	    case 'all'
				try chans = setdiff(1:data.num_sensors, find(strcmp('other',{data.sensor_info.typestring}))); end;
		end
	end
	if ~exist('chans','var')
		if isfield(data.(DATA_FIELD)(1),'data')
			chans = 1:size(data.(DATA_FIELD)(1).data,1);
		elseif isfield(data.(DATA_FIELD)(1),'power')
			chans = 1:size(data.(DATA_FIELD)(1).power,1);
		end
	end 
else
  chans = parms.channels;
end;

chans = setdiff(chans,badchan_i);
if isempty(chans)
  mmil_error(parms, 'no good channels selected');
end;

% Grab the data for this condition
cond_data = getfield(data, DATA_FIELD);
cond_data = cond_data(parms.condition);

%
t_trigger = cond_data.time(1)*1000;

% set start and end samples
t0 = cond_data.time(1)*1000;
t1 = cond_data.time(end)*1000;

% Massage time0
if isempty(parms.time0)
  parms.time0 = t0;
elseif parms.time0<t0
  parms.time0=t0; 
  mmil_logstr(parms, 'WARNING: time0 < dataset''s lowest time; resetting time0=t0.');
end;

% Massage time1
if isempty(parms.time1)
  parms.time1 = t1;
elseif time1>t1
  parms.time1=t1; 
  mmil_logstr(parms, 'WARNING: time1 > dataset''s highest time; resetting time1=t1');
end;

% Massage the relation between time0 and time1
if parms.time0>parms.time1
  t0 = time0;
  t1 = time1;
  parms.time1 = t0;
  parms.time0 = t1;
  mmil_logstr(parms, 'WARNING: time0>time1; resetting to be the same.');
end;

clear('t0','t1');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data assignments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



try
	FT_data.label = {data.sensor_info(chans).label};
catch
	labels = cell(1,length(chans));
	for i=1:length(labels), labels{i} = num2str(i); end;
	FT_data.label = labels;
end
FT_data.fsample = data.sfreq;

if isempty(parms.channels) && ~strcmp(DATA_FIELD,'stats')
  switch parms.chantype
    case {'mag' 'grad1' 'grad2' 'grad', 'other', 'meg'}

      % see read_fcdc_elec for field defs

      FT_data.grad.pnt   = zeros(length(chans),3);
      FT_data.grad.label = cell(length(chans),1);
%        FT_data.grad.tra   = ones(length(chans), length(chans));

%        FT_data.grad.ori   = []; % this simply tricks FT into thinking 
                               % these 'grads' are MEG (and not EEG) grads
                               % This SHOULD be Nx3 specifying
                               % 'orientations'

    case {'eeg' 'ieeg'}
      FT_data.elec.pnt   = zeros(length(chans),3);
      FT_data.elec.label = cell(length(chans),1);
  end;

  % Apply coordinate transform to sensor data
  for c=1:length(chans)
    k=chans(c);
    try loc = data.sensor_info(k).loc; end;
    switch parms.chantype
      case {'mag' 'grad1' 'grad2' 'grad', 'meg','other'}
        % apply device2head transform
        if ~isempty(data.coor_trans.device2head)
          T = data.coor_trans.device2head;
          loc = T*loc;
        end;
        FT_data.grad.label{c} = data.sensor_info(k).label;
        FT_data.grad.pnt(c,1:3) = loc(1:3,4);
      case {'eeg' 'ieeg'}
        % locations already in "head" space
        FT_data.elec.label{c} = data.sensor_info(k).label;
        FT_data.elec.pnt(c,1:3) = loc(1:3,4);
    end;
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIMORD: not common
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since t0 and t1 need to be 1-based indices,
% we MUST add 1 to these numbers.
t0 = round((parms.time0 - t_trigger)*data.sfreq/1000) + 1; 
t1 = round((parms.time1 - t_trigger)*data.sfreq/1000) + 1;

if (t0 < 1), t0 = 1; end;
if (t1>length(cond_data.time)), t1=length(cond_data.time); end;

switch DATA_FIELD
  case 'averages'
    FT_data.avg  = cond_data.data(chans,t0:t1);
    try FT_data.var  = cond_data.stdev(chans,t0:t1); end
    FT_data.time = cond_data.time(t0:t1);
    if isfield(cond_data,'mask')
      FT_data.mask = cond_data.mask(chans,t0:t1);
    end
  case 'epochs'
    switch parms.dimord
      case 'trial_chan_time'
        ntrials = cond_data.num_trials;
        nchans = length(chans);
        tpoints = length(cond_data.time(t0:t1));

        FT_data.time  = cond_data.time(t0:t1);
        FT_data.trial = permute(cond_data.data(chans,t0:t1,:),[3 1 2]);

      case 'chan_time'
        FT_data.trial  = squeeze(num2cell(cond_data.data(chans,t0:t1,:), [1 2]));
        FT_data.time   = repmat({cond_data.time(t0:t1)}, [cond_data.num_trials 1]);
    end;

  case 'timefreq' % time / freq / powspctrm
    FT_data.time      = cond_data.time(t0:t1);
    FT_data.freq      = cond_data.frequencies;          
    if isfield(cond_data,'labelcmb')
      FT_data.labelcmb = cond_data.labelcmb;
    end
    switch parms.dimord
        case 'chan_freq_time'
          if (isfield(cond_data, 'power'))                                        % chan time freq
            FT_data.powspctrm = permute(double(cond_data.power(chans,t0:t1,:)), [1 3 2]);
          elseif (isfield(cond_data,'cmplx'))
%               FT_data.fourierspctrm = permute(cond_data.cmplx(chans,t0:t1,:),[1 3 2]);
% 							FT_data.powspctrm = permute(2*abs(double(cond_data.cmplx(chans,t0:t1,:))).^2, [1 3 2]);
            FT_data.powspctrm = permute(abs(double(cond_data.cmplx(chans,t0:t1,:))).^2, [1 3 2]);
          elseif (isfield(cond_data,'data'))
%               FT_data.powspctrm = permute(2*abs(double(cond_data.data(chans,t0:t1,:))).^2, [1 3 2]);
            FT_data.powspctrm = permute(abs(double(cond_data.data(chans,t0:t1,:))).^2, [1 3 2]);
          end;
          if isfield(cond_data,'cmplx')
%               FT_data.fourierspctrm = permute(cond_data.cmplx(chans,t0:t1,:),[1 3 2]);
          end           
          if isfield(cond_data,'mask'),  FT_data.mask      = permute(cond_data.mask(:,t0:t1,:),[1 3 2]); end
          if isfield(cond_data,'cross'), FT_data.crsspctrm = permute(cond_data.cross(:,t0:t1,:),[1 3 2]); end
          if isfield(cond_data,'coh'),   FT_data.cohspctrm = permute(cond_data.coh(:,t0:t1,:),[1 3 2]); end         
          if isfield(cond_data,'plv'),   FT_data.plvspctrm = permute(cond_data.plv(:,t0:t1,:),[1 3 2]); end      
        case 'rpt_chan_freq_time'
          if (isfield(cond_data, 'power'))
              FT_data.powspctrm = permute(double(cond_data.power(chans,t0:t1,:,:)),[4 1 3 2]);
          elseif (isfield(cond_data,'cmplx'))
% 								FT_data.powspctrm = permute(2*abs(double(cond_data.cmplx(chans,t0:t1,:,:))).^2, [4 1 3 2]);
              FT_data.powspctrm = permute(abs(double(cond_data.cmplx(chans,t0:t1,:,:))).^2, [4 1 3 2]);
          elseif (isfield(cond_data,'data'))
%                 FT_data.powspctrm = permute(2*abs(double(cond_data.data(chans,t0:t1,:,:))).^2,[4 1 3 2]);
              FT_data.powspctrm = permute(abs(double(cond_data.data(chans,t0:t1,:,:))).^2,[4 1 3 2]);
          end
          if isfield(cond_data,'cmplx')
            try
%                 FT_data.fourierspctrm = permute(cond_data.cmplx(chans,t0:t1,:,:),[4 1 3 2]);
            end
          end
          if isfield(cond_data,'cross'), FT_data.crsspctrm = permute(cond_data.cross(:,t0:t1,:,:),[4 1 3 2]); end
          if isfield(cond_data,'coh'),   FT_data.cohspctrm = permute(cond_data.coh(:,t0:t1,:,:),[4 1 3 2]); end         
          if isfield(cond_data,'plv'),   FT_data.plvspctrm = permute(cond_data.plv(:,t0:t1,:,:),[4 1 3 2]); end  
    end
  case 'stats'
    FT_data.time = cond_data.time(t0:t1);          
    cond_data = rmfield(cond_data,{'event_code','label'});
    switch parms.dimord
        case 'chan_freq_time'
          FT_data.freq = cond_data.frequencies;
          fieldlist = fieldnames(cond_data);
          for f=1:length(fieldlist)
              FT_data.(fieldlist{f}) = cond_data.(fieldlist{f});
          end
          FT_data = rmfield(FT_data,'frequencies');            
          try FT_data.prob = permute(FT_data.prob(chans,t0:t1,:),[1 3 2]); end
          try FT_data.stat = permute(FT_data.stat(chans,t0:t1,:),[1 3 2]); end        
          try FT_data.mask = permute(FT_data.mask(chans,t0:t1,:),[1 3 2]); end
          try FT_data.avg  = permute(FT_data.mask(chans,t0:t1,:),[1 3 2]); end
          try FT_data.posclusterslabelmat = permute(FT_data.posclusterslabelmat(chans,t0:t1,:),[1 3 2]); end
          try FT_data.negclusterslabelmat = permute(FT_data.negclusterslabelmat(chans,t0:t1,:),[1 3 2]); end            
        case 'chan_time'
          fieldlist = fieldnames(cond_data);
          for f=1:length(fieldlist)
              FT_data.(fieldlist{f}) = cond_data.(fieldlist{f});
          end
          try FT_data.prob = FT_data.prob(chans,t0:t1); end
          try FT_data.stat = FT_data.stat(chans,t0:t1); end
          try FT_data.mask = FT_data.mask(chans,t0:t1); end
          try FT_data.avg  = FT_data.mask; end
          try FT_data.posclusterslabelmat = FT_data.posclusterslabelmat(chans,t0:t1); end
          try FT_data.negclusterslabelmat = FT_data.negclusterslabelmat(chans,t0:t1); end              
    end
end;

FT_data.dimord = parms.dimord;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /DIMORD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  % get the full name of the function
  FT_data.cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  FT_data.cfg.version.name = st(i);
end

% add the version details of this function call to the configuration
FT_data.cfg.version.id = sprintf('$Id: %s.m,v 0.1 2007/08/08 benc Exp $',mfilename);

