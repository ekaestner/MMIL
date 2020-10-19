function [data,samp,cond] = ts_MNE_epoch_fif(varargin)
% Input continuous data (TimeSurfer epoch_data structure with one trial)
%
% Epoch based on latencies from:
%   - Eventfile (*.ev2, events_fnames)
%   - Raw data file (based on trigger channel)
%   - Vector of samples
%   - Vector of time points
%
% Read epoched data using an evntfile:
%   datafile = '/somepath/cl_screening_v2_032010.edf';
%   evntfile = '/somepath/cl_screening_v2_032010.evt';
%   prestim  = .2;
%   poststim = .6;
%   data     = ts_loadedf('datafile',datafile,'evntfile',evntfile,'prestim',prestim,'poststim',poststim);
%
% Example evntfile:
% ----------------------------------------
% 15-Apr-2010
% filename= cl_screening_v2_032010.edf
% SR= 1998.8242
% sample point	stim_num	stim_category
% 10682	0	fixation
% 20712	1	face
% 22744	2	object
% 24710	2	object
% 26710	1	face
% ----------------------------------------
% Note: SR = sampling rate (Hz).
% This would produce epoch_data with three conmod(length(varargin),2)ditions (event codes 0,1,2).
% The first condition would have one trial and the other two would have two
% each.
%
% Troubleshooting:
% If you have a problem with sopen.m, you need a newer version of biosig
% which is part of EEGLAB.  Solution: download the latest version of EEGLAB.
%
% If you have a conflict between biosig and matlab versions of str2double.m:
% make a copy of this function, uncomment the commands below (calls to
% addpath, unix, and the lines setting mv1 & mv2). Also, in those commands,
% change the paths to str2double.m.
%
% Created:  04/01/10 by Jason Sherfey
% Rcnt Mod: 08/25/10 by Jason Sherfey
% Last Mod: 09/15/12 by Don Hagler
%

% addpath(genpath('/home/jsherfey/svn/dev/packages/eeglab7_2_9_18b'));
addpath /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200
addpath /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t250
mv1      = 'mv /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double_rmbyJSS.m /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double.m';
mv2      = 'mv /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double.m /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double_rmbyJSS.m';

% if mod(length(varargin),2) ~= 0
%   % odd number of inputs
%   varargin = {'datafile',varargin{:}};
% end

parms = mmil_args2parms( varargin, ...
    { 'datafile',[],[],...
    'datatype',[],[],...
    'evntfile',[],[],...
    'evnttype',[],[],...
    'channels',[],[],...
    'stim_delay',0,[],...
    'prestim' ,.2,[],...
    'poststim',.6,[],...
    'trigchan',[],[],...
    'trigthresh',[],[],...
    'trig_minduration',[],[],...
    'trig_minseparation',.01,[],...
    'oldchanlabels',[],[],...
    'newchanlabels',[],[],...
    'continuous',0,[],...
    'valid_event_codes',[],[],...
    'write_log_flag',0,[],...
    'rootoutdir',pwd,[],...
    'second_flag',0,[],...
    'samp_vector',[],[],...
    'trigchan_flag',0,[],...
    }, ...
    false );
if isequal(parms.continuous,'yes'), parms.continuous = 1; end
if isequal(parms.continuous,'no'),  parms.continuous = 0; end
% if ~ischar(parms.evntfile) && ~iscellstr(parms.trigchan), parms.continuous = 1; end
if isempty(parms.trig_minseparation), parms.trig_minseparation = parms.stim_delay; end







% use evntfile to epoch the data
% read edf header
% unix(mv1);
% hdr = read_header(parms.datafile);
% unix(mv2);
% nsamp = hdr.nSamples * hdr.nTrials;
% nchan = hdr.nChans;




% read event file (assumes format described in help section; different
% from ev2 and events_fnames)

if isempty(parms.evnttype) && ~isempty(parms.evntfile)
    
    if strcmp(parms.evntfile(end-3:end),'.txt')
        % .txt are evntfname files
        parms.evnttype='fname';
    elseif strcmp(parms.evntfile(end-3:end),'.ev2')
        parms.evnttype='ev2';
    elseif strcmp(parms.evntfile(end-3:end),'.edf')
        parms.evnttype='edf';
    end;
    
end;

if isempty(parms.datatype)
    
    if strcmp(parms.datafile(end-3:end),'.fif')
        % .txt are evntfname files
        parms.datatype='fif';
    elseif strcmp(parms.datafile(end-3:end),'.edf')
        parms.datatype='edf';
    elseif strcmp(parms.datafile(end-3:end),'.mat')
        parms.datatype='mat';
    end;
end;

if strcmp(parms.datatype,'fif')
    unix(mv1);
hdr = read_header(parms.datafile);
unix(mv2);
nsamp = hdr.nSamples * hdr.nTrials;
nchan = hdr.nChans;
    data = ts_MNE_loadfif(parms.datafile);
elseif strcmp(parms.datatype,'edf');
    unix(mv1);
hdr = read_header(parms.datafile);
unix(mv2);
nsamp = hdr.nSamples * hdr.nTrials;
nchan = hdr.nChans;
    data = ts_loadedf(parms.datafile);
    
elseif strcmp(parms.datatype,'mat');
    load(parms.datafile);
%     data=epoch_data;
%     clear epoch_data
    nchan=length(data.sensor_info);
else
    error('this is an unsupported datatype');
end;





if isempty(parms.channels)
    chans = 1:nchan;data
else
    chans = parms.channels;
end;


if isempty(parms.evntfile) || parms.continuous || parms.trigchan_flag
  % convert continuous data to timesurfer format
%   data = ts_fieldtrip2data(dat,'epoch');
  % fix channel labels
  if iscellstr(parms.newchanlabels)
    [sel1,sel2] = match_str({data.sensor_info.label},parms.oldchanlabels);
    [data.sensor_info(sel1).label] = deal(parms.newchanlabels{:});
    clear sel1 sel2
  end
  % replace NaNs with zeros (note: FieldTrip preproc was modified to return
  % unfiltered data without NaNing the entire data set first).
  if any(isnan(data.epochs.data(:)))
    data.epochs.data(isnan(data.epochs.data)) = 0;
  end
  if ~isempty(parms.trigchan) && ~parms.continuous
    % assume these channels carry a binary code for epoching
    Fs = data.sfreq;
    % select trigger data
    trigdata     = ts_data_selection(data,'chanlabel',parms.trigchan,'verbose',0);
    % set threshold for trigger detection
    if isempty(parms.trigthresh) || ~isnumeric(parms.trigthresh)
      tmp              = trigdata.epochs.data(:);
      parms.trigthresh = mean(tmp) + std(tmp);
      clear tmp
    end
    % find threshold crossings
    t     = data.epochs.time;
    nsamp = length(t);
    ntrigchan = trigdata.num_sensors;
    [crsind{1:ntrigchan}] = deal([]);
    for k = 1:ntrigchan
      tmpdat  = trigdata.epochs.data(k,:);
      ind     = crossing(tmpdat,[],parms.trigthresh);
      % only keep left edge of pulse
      sel     = (tmpdat(ind+1)-tmpdat(ind-1)) > 0;
      ind     = ind(sel);
      % only keep crossings that occur > minduration since last crossing, if specified
      if ~isempty(parms.trig_minduration)
        sel = [1 find((t(ind(2:end))-t(ind(1:end-1))) >= parms.trig_minduration)];
        ind = ind(sel);
      end
      crsind{k} = ind;
      clear ind sel tmpdat
    end
    bincode = zeros(ntrigchan,nsamp);
    for k   = 1:length(crsind), bincode(k,crsind{k}) = 1; end
    detind  = any(bincode,1); % whether there's a detection in any trigchan
    % make sure detections for the same pulse are grouped before conversion
      delta   = round(parms.trig_minseparation*Fs/2);
      % get indices at which a detection occurs in any trigchan
      ind     = find(detind);
      % find detections with a second detection > delta steps before it
      tmp     = cellfun(@(x)any((ind>x-delta & (ind<x))),num2cell(ind));
      R_ind   = find(tmp);
      % find detections with a second detection < delta steps after it
      tmp     = cellfun(@(x)any((ind<x+delta & (ind>x))),num2cell(ind));
      L_ind   = find(tmp);
      % remove detections that are between two other detections in the same pulse
      Rsel    = ~ismember(R_ind,L_ind);
      Lsel    = ~ismember(L_ind,R_ind);
      R_ind   = R_ind(Rsel); clear Rsel
      L_ind   = L_ind(Lsel); clear Lsel
      % for each pair (L,R), set [ind(L):ind(R)] to 1 in each bincode row w/ a detection
      for k = 1:ntrigchan
        sel = cellfun(@(x,y)any(bincode(k,ind(x):ind(y))),num2cell(L_ind),num2cell(R_ind));
        sel = cellfun(@(x,y)ind(x):ind(y),num2cell(L_ind(sel)),num2cell(R_ind(sel)),'uniformoutput',false);
        sel = unique([sel{:}]);
        bincode(k,sel) = 1;
      end
    detind  = any(bincode,1); % whether there's a detection in any trigchan    
    % only keep the first pulse in contiguous detections
    samp    = find(detind(1:end-1)==0 & detind(2:end)==1) + 1;
    bincode = flipud(bincode(:,samp));
    bincode = mat2cell(bincode,ntrigchan,ones(1,length(samp)));
    % convert binary to decimal
    evcodes = cellfun(@(x)polyval(x,2),bincode);
    clear detind

    % keep valid_event_codes
    if ~isempty(parms.valid_event_codes)
      samp    = samp(ismember(evcodes,parms.valid_event_codes));
      evcodes = evcodes(ismember(evcodes,parms.valid_event_codes));
    else
      parms.valid_event_codes = evcodes;
    end
    % offset samples for nonzero stimulus delays
    samp      = samp + round(parms.stim_delay * Fs);

    % prepare event info 
    keep = ~((samp-parms.prestim*Fs)<1 | (samp+parms.poststim*Fs)>nsamp); % not outbounds
    samp = samp(keep);
    cond = evcodes(keep);
    [conds,ii]  = unique(cond);
    [dura(1:length(samp))] = deal(1);
    nevnt       = length(samp);
    ncond       = length(conds);
    begsample   = round(samp - parms.prestim*Fs );
    endsample   = round(samp + parms.poststim*Fs);
    % % CHECK
    % tmp       = zeros(1,length(t)); 
    % tmp(samp) = cond;
    % data.epochs.data(end,:) = tmp;
    % visualizer(data); clear tmp
    % % END CHECK
  end;
elseif strcmp(parms.evnttype,'ev2');
    [content,res] = mmil_readtext(parms.evntfile,'\t');
    content       = cell2mat(content);
    Fs  = hdr.Fs;
    %          T    = [-parms.prestim:1/Fs:parms.poststim];
    samp = content(:,6);
    if (parms.second_flag)
        samp=samp*Fs;
    end;
    keep = ~((samp-parms.prestim*Fs)<1 | (samp+parms.poststim*Fs)>nsamp); % not outbounds
    samp = samp(keep);
    cond = content(keep,2);
    %name = evnt(keep,1);
    [dura(1:length(samp))] = deal(1);
    [conds,ii]  = unique(cond);
    nevnt       = length(samp);
    ncond       = length(conds);
    begsample   = round(samp - parms.prestim*Fs );
    endsample   = round(samp + parms.poststim*Fs);
%Read event_info from Fnames file.
elseif strcmp(parms.evnttype,'fname');
    [content,res] = mmil_readtext(parms.evntfile,'\t');
    Fs  = hdr.Fs;
    %          T    = [-parms.prestim:1/Fs:parms.poststim];
    evnt = content(2:end,:);
    samp = cell2mat(evnt(:,2));
    keep = ~((samp-parms.prestim*Fs)<1 | (samp+parms.poststim*Fs)>nsamp); % not outbounds
    samp = samp(keep);
    cond = [evnt{keep,3}]';
    name = evnt(keep,1);
    dura = evnt(keep,4);
    [conds,ii]  = unique(cond);
    nevnt       = length(samp);
    ncond       = length(conds);
    begsample   = round(samp - parms.prestim*Fs );
    endsample   = round(samp + parms.poststim*Fs);
end;





contdata         = data.epochs.data;
data.epochs.time = [];
data.epochs.data = [];
data.epochs(1:ncond) = data.epochs;
% epoch the continuous data

if isempty(parms.channels)
    chans = 1:size(contdata,1);
else
    chans = parms.channels;
end
T     = [-parms.prestim:1/Fs:parms.poststim];
%     fdataor c = 1:ncond
%       ii  = cond==conds(c);
%       s0  = num2cell(begsample(ii));
%       sf  = num2cell(endsample(ii));
%       tmp = cellfun(@(x,y)(contdata(chans,x:y)),s0,sf,'UniformOutput',false);
%       data.epochs(c).event_code = conds(c);
%       data.epochs(c).time       = T;
%       data.epochs(c).num_trials = length(s0);
%       data.epochs(c).data       = cat(3,tmp{:}); % channel x time x trial
%       clear tmp
%     clear contdata

for c = 1:ncond
    ii  = cond==conds(c);
    s0  = num2cell(begsample(ii));
    sf  = num2cell(endsample(ii));
    tmp = cellfun(@(x,y)(contdata(chans,x:y)),s0,sf,'UniformOutput',false);
    data.epochs(c).event_code = conds(c);
    data.epochs(c).time       = T;
    data.epochs(c).num_trials = length(s0);
    data.epochs(c).data       = cat(3,tmp{:}); % channel x time x trial
    % Build the trial info field.
    data.epochs(c).trial_info={};
    data.epochs(c).trial_info.number                     = find(ii)';
    data.epochs(c).trial_info.latency                    = samp(ii)';
    data.epochs(c).trial_info.badtrial                   = zeros(1,length(s0));
    data.epochs(c).trial_info.event_code                 = cond(ii)';
    data.epochs(c).trial_info.duration                   = dura(ii)';
    [data.epochs(c).trial_info.datafile{1,1:length(s0)}] = deal(parms.datafile);
    [data.epochs(c).trial_info.events_fnames{1,1:length(s0)}] = deal(parms.evntfile);
    clear tmp
end;
clear contdata





%   % fix channel labels
%   if iscellstr(parms.newchanlabels)
%     [sel1,sel2] = match_str({data.sensor_info.label},parms.oldchanlabels);
%     [data.sensor_info(sel1).label] = deal(parms.newchanlabels{:});
%     clear sel1 sel2
%   end
%   % replace NaNs with zeros (note: FieldTrip preproc was modified to return
%   % unfiltered data without NaNing the entire data set first).
%   if any(isnan(data.epochs.data(:)))
%     data.epochs.data(isnan(data.epochs.data)) = 0;
%   end
% end
if parms.write_log_flag && exist('samp','var')
    [pathstr, name, ext] = fileparts(parms.datafile);
    logfile = sprintf('%s/%s_events.log',parms.rootoutdir,name);
    fid     = fopen(logfile,'a');
    [~,b]   = unix('echo $HOST');
    fprintf(fid,'--------------------------------------------------------------------\n');
    fprintf(fid,'Time: %s\n',datestr(now));
    fprintf(fid,'Host: %s',b);
    fprintf(fid,'Importing %s data:\n%s\n',ext,parms.datafile);
    fprintf(fid,'Number channels: %i \n',data.num_sensors);
    fprintf(fid,'Sampling frequency: %5.6g Hz\n',Fs);
    fprintf(fid,'--------------------------------------------------------------------\n');
    fprintf(fid, '%-6s\t %-10s\t %-10s\t %-12s\t %-6s\t %-20s\n','Trial','Event','Latency','Time (sec)','Type','Comment');
    for k = 1:nevnt
        fprintf(fid,'%-6g\t %-10g\t %-10g\t %-12.8g\t %-6s\t %-20s\n',...
          k,cond(k),samp(k),samp(k)/Fs,'ok','no rejection');
    end
    fprintf(fid,'--------------------------------------------------------------------\n');
    fclose(fid);
end;



function build_trial_info

