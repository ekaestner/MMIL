function [data,samp,cond] = ts_loadedf(varargin)
% Read continuous data:
%   datafile = '/somepath/cl_screening_v2_032010.edf';
%   data     = ts_loadedf('datafile',datafile);
%   save(outfile,'data','-v7.3');
% 
% Read epoched data from NYU clinical system w/ binary code in DC channels:
% Example 1:
%   datafile = '/somepath/cl_screening_v2_032010.edf';
%   trigchan = {'Pulse DC1','Pulse DC2','Pulse DC3','Pulse DC4','Pulse DC5'};
%   prestim    = .2;
%   poststim   = .6;
%   stim_delay = .06;
%   data       = ts_loadedf('datafile',datafile,'prestim',prestim,'poststim',poststim,'stim_delay',stim_delay);
%   COMMENT: If several trigger channels are listed, the function assumes
%   they carry a binary code. Specifying the same datafile without a cell
%   array of trigger channel labels will return the continuous data.
% 
% Example 2:
%   datafile  = '/home/ccarlson/data/eeg/FW_clinical_sys_100819/FW_DH.edf';
%   oldlabels = {'EEG GA64_01-Ref','EEG GA64_02-Ref','EEG GA64_05-Ref','EEG GA64_06-Ref',...
%                'EEG GA64_11-Ref','EEG GA64_12-Ref','EEG GA64_15-Ref','EEG GA64_16-Ref',...
%                'EEG GA64_17-Ref'};
%   newlabels   = {'Cz','CPz','Oz','O1','PO7','P7','O2','PO8','P8'};
%   trigchan    = {'Pulse DC1','Pulse DC2','Pulse DC3','Pulse DC4','Pulse DC5'};
%   stim_delay  = .06;
%   prestim     = .25; 
%   poststim    = .5;
%   trig_minduration   = [];
%   trig_minseparation = []; % set to 10ms by default
%   valid_event_codes  = [3 4 5 6 7];
%   data = ts_loadedf('datafile',datafile,'old_chan_labels',oldlabels,'new_chan_labels',newlabels,...
%                     'trigchan',trigchan,'stim_delay',stim_delay,'prestim',prestim,'poststim',poststim,...
%                     'trig_minduration',trig_minduration,'valid_event_codes',valid_event_codes);
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
% This would produce epoch_data with three conditions (event codes 0,1,2).
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
% Created on April-2010 by JSS
% Modified  on 25-Aug-2010 by JSS
% Modified: on 29-Mar-2012 by SRD: Fix problem calling fieldtrip read_data in the new eeglab package and to remove dependency on JSS home directory.
% Modified: on 06-April-2012 by SRD: Corrected problem with path not fully restored after the path modifications, and renamed crossing ts_crossing.
% Last Mod: on 29-Aug-2012 by SRD: Added 'no_ofd' to force read_biosig_header & read_biosig_data to not detect overflow.  Passed thru read_header & data ... to sopen.

compatible_eeglab = '/usr/pubsw/packages/eeglab/eeglab9_0_5_6b';  %%% added, SRD 032612;
compatible_fieldtrip = '/usr/pubsw/packages/MMPS/external/matlab/fieldtrip-20080624';   %%% added, SRD 032612
incompatible_eeglab_subset = '/usr/pubsw/packages/eeglab/eeglab9_0_5_6b/external/fieldtrip-partial/fileio/compat';  %%% added, SRD 040612
overflow_detect_option = '/usr/pubsw/packages/MMPS/matlab/timesurfer/external_mods/fieldtrip';   %%% added, SRD 083012 To force off overflow detect when reading edf.
% overflow_detect_option = '/home/sdeiss/matlab/svn/MMPS_test_patches';  %%% SRD 083012 added overflow_detect flag, allows disabling overflow detect when reading edf.

% addpath(genpath('/home/jsherfey/svn/dev/packages/eeglab7_2_9_18b'));  %%%%% Removed these lines after above additions, SRD 032612
% % % addpath /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200
% % % addpath /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t250
% % % mv1      = 'mv /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double_rmbyJSS.m /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double.m';
% % % mv2      = 'mv /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double.m /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double_rmbyJSS.m';

if mod(length(varargin),2) ~= 0
    % odd number of inputs
    varargin = {'datafile',varargin{:}};
end

parms = mmil_args2parms( varargin, ...
    { 'datafile',[],[],...
    'evntfile',[],[],...
    'channels',[],[],...
    'stim_delay',0,[],...
    'prestim' ,.2,[],...
    'poststim',.6,[],...
    'trigchan',[],[],...
    'trigthresh',[],[],...
    'trig_minduration',[],[],...
    'trig_minseparation',.01,[],...
    'old_chan_labels',[],[],...
    'new_chan_labels',[],[],...
    'continuous',0,[],...
    'valid_event_codes',[],[],...
    'overflow_detect',1,[],...
    }, ...
    false );
if isequal(parms.continuous,'yes'), parms.continuous = 1; end
if isequal(parms.continuous,'no'),  parms.continuous = 0; end
if ~ischar(parms.evntfile) && ~iscellstr(parms.trigchan), parms.continuous = 1; end
if isempty(parms.trig_minseparation), parms.trig_minseparation = parms.stim_delay; end
% EDF file
if isempty(parms.datafile)
    % ask user
    [filename,pathname] = uigetfile({'*.mat;*.fif;*.eeg;*.avg;*.cnt;*.vhdr;*.set'},'Pick a file.','MultiSelect','on');
    if isequal(filename,0) || isequal(pathname,0)  return; end
    if iscell(filename)
        parms.datafile = cellfun(@(x)fullfile(pathname,x),filename,'uniformoutput',false);
    else
        parms.datafile = [pathname filename];
    end
end
if iscell(parms.datafile)
    if length(parms.datafile) == 1
        parms.datafile = parms.datafile{1};
    else
        error('%s: support is limited to importing one file only.',mfilename);
    end
end

% read edf header
% % % % unix(mv1);
LePath = path; %%%%% changed, SRD 040312  Save the original path.
addpath(genpath(compatible_eeglab)); %%%%% changed, SRD 032612
rmpath(incompatible_eeglab_subset); %%%%% SRD 040612
addpath(genpath(compatible_fieldtrip));
addpath(genpath(overflow_detect_option)); %%%%% SRD 083012
hdr = read_header(parms.datafile, 'overflow_detect', parms.overflow_detect); %%%%% SRD 083012
path(LePath); %%%%% changed, SRD 040312  Restore the original path.
% % % % unix(mv2);

% read continuous data
cfg = [];
cfg.dataset    = parms.datafile;
cfg.continuous = 'yes';
if ~isempty(parms.channels)
    cfg.channel = hdr.label(parms.channels);
    parms.channels = 1:length(cfg.channel);
else
    cfg.channel = hdr.label;
end

% % % % unix(m v1);
LePath = path; %%%%% changed, SRD 040312
addpath(genpath(compatible_eeglab)); %%%%% changed, SRD 032612
rmpath(incompatible_eeglab_subset); %%%%% SRD 040612
addpath(genpath(compatible_fieldtrip));
addpath(genpath(overflow_detect_option)); %%%%% SRD 083012
cfg.overflow_detect = parms.overflow_detect;  %%%%% SRD 083012
dat  = preprocessing(cfg); %%%%% SRD 083012
path(LePath); %%%%% changed, SRD 040312
% % % % unix(mv2);

if isempty(parms.evntfile) || parms.continuous
    % convert continuous data to timesurfer format
    data = ts_fieldtrip2data(dat,'epoch');
    % fix channel labels
    if iscellstr(parms.new_chan_labels)
        [sel1,sel2] = match_str({data.sensor_info.label},parms.old_chan_labels);
        [data.sensor_info(sel1).label] = deal(parms.new_chan_labels{:});
        clear sel1 sel2
    end
    % replace NaNs with zeros (note: FieldTrip preproc was modified to return
    % unfiltered data without NaNing the entire data set first).
    if any(isnan(data.epochs.data(:)))
        data.epochs.data(isnan(data.epochs.data)) = 0;
    end
    if iscellstr(parms.trigchan) && ~parms.continuous
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
            ind     = ts_crossing(tmpdat,[],parms.trigthresh);
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
        for c = 1:ncond
            ii  = cond==conds(c);
            s0  = num2cell(begsample(ii));
            sf  = num2cell(endsample(ii));
            tmp = cellfun(@(x,y)(contdata(chans,x:y)),s0,sf,'UniformOutput',false);
            data.epochs(c).event_code = conds(c);
            data.epochs(c).time       = T;
            data.epochs(c).num_trials = length(s0);
            data.epochs(c).data       = cat(3,tmp{:}); % channel x time x trial
            clear tmp
        end
        clear contdata
    end
elseif ~exist(parms.evntfile,'file')
    error('Event file does not exist.');
else
    % use evntfile to epoch the data
    % read edf header
    % % % %   unix(mv1);
    LePath = path; %%%%% changed, SRD 040312
    addpath(genpath(compatible_eeglab)); %%%%% changed, SRD 032612
    rmpath(incompatible_eeglab_subset); %%%%% SRD 040612
    addpath(genpath(compatible_fieldtrip));
    addpath(genpath(overflow_detect_option)); %%%%% SRD 083012
    hdr = read_header(parms.datafile, 'overflow_detect', parms.overflow_detect); %%%%% SRD 083012
    path(LePath); %%%%% changed, SRD 040312
    % % % %   unix(mv2);
    nsamp = hdr.nSamples * hdr.nTrials;
    nchan = hdr.nChans;
    if isempty(parms.channels)
        chans = 1:nchan;
    else
        chans = parms.channels;
    end
    % read event file
    [content,res] = readtext(parms.evntfile,'\t');
    Fs   = content{3,1};
    Fs   = Fs(regexp(Fs,'[\d.]'));
    Fs   = str2num(Fs);
    T    = [-parms.prestim:1/Fs:parms.poststim];
    
    evnt = content(res.numberMask(:,1),:);
    samp = [evnt{:,1}]';
    keep = ~((samp-parms.prestim*Fs)<1 | (samp+parms.poststim*Fs)>nsamp); % not outbounds
    samp = samp(keep);
    cond = [evnt{keep,2}]';
    name = evnt(keep,3);
    
    [conds,ii] = unique(cond);
    names      = name(ii);
    nevnt      = length(samp);
    ncond      = length(conds);
    
    begsample = round(samp - parms.prestim*Fs );
    endsample = round(samp + parms.poststim*Fs);
    
    % initialize timesurfer structure
    tmp          = rmfield(dat,{'trial','time'});
    tmp.trial{1} = dat.trial{1}(:,1);
    tmp.time{1}  = dat.time{1}(1);
    data         = ts_fieldtrip2data(tmp,'epoch');
    data.epochs.time = [];
    data.epochs.data = [];
    data.epochs(1:ncond) = data.epochs;
    clear tmp
    
    % epoch the continuous data
    for c = 1:ncond
        ii  = cond==conds(c);
        s0  = num2cell(begsample(ii));
        sf  = num2cell(endsample(ii));
        tmp = cellfun(@(x,y)(dat.trial{1}(chans,x:y)),s0,sf,'UniformOutput',false);
        data.epochs(c).event_code = conds(c);
        data.epochs(c).time       = T;
        data.epochs(c).num_trials = length(s0);
        data.epochs(c).data       = cat(3,tmp{:}); % channel x time x trial
        clear tmp
    end
    % fix channel labels
    if iscellstr(parms.new_chan_labels)
        [sel1,sel2] = match_str({data.sensor_info.label},parms.old_chan_labels);
        [data.sensor_info(sel1).label] = deal(parms.new_chan_labels{:});
        clear sel1 sel2
    end
    % replace NaNs with zeros (note: FieldTrip preproc was modified to return
    % unfiltered data without NaNing the entire data set first).
    if any(isnan(data.epochs.data(:)))
        data.epochs.data(isnan(data.epochs.data)) = 0;
    end
end
