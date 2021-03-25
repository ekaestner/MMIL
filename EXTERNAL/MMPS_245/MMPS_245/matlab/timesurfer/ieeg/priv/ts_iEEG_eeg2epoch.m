function [epoch_data] = ts_iEEG_eeg2epoch (data_files,varargin)

% Use [epoch_data] = ts_iEEG_eeg2epoch (data_files,'badchanfile','badchanfile.txt','channamefile','chan_names.txt');
%
% Converts following file formats to the TimeSurfer
% format:
%
%  Neuroscan .eeg file
%  BrainVision .vhdr file
%
% Required Input:
%
%   data_files: full path file name (or cell array of file names) of neuroscan
%     format eeg file, BrainVision .vhdr file, EEGLAB .set file, or Matlab .mat file
%
%
% Optional Input:
%
%   badchanfile  - full path and file name of text file containing the list
%   of badchannels.  bad channels will be marked as such in the epoch data
%   structure
%   channamefile - full path and file name of a text file containing the
%   list of channel names to replace the eeg files channel names with.
%   Must be in same order as the channels in the data file.
%   timelimits    - a 2 element vector with start and end time in seconds centered
%                   around the events in the BrainVision file (default =
%                   [-2 2])
%   noise_start   - time in sec of start of noise (default = -.08 sec)
%   noise_end     - time in sec of end of noise (default = 0 sec)
%   recode_rules  - string or cell array of strings, of recoding rules
%     Each rule specifies a sequence of events and the new event code which is
%       assigned to the first event in the seqeuence
%     Examples:
%       '2->5=101' : event 2 followed directly by event 5 gets recoded as event 101
%       '[2,3,4]->5=101' : events 2, 3, or 4, followed by event 5
%       '[2]->[5,10]=101' : event 2 followed by event 5 or 10
%       '2->200:1500[5]=101' : event 2 followed by event 5 within 200 and 1500 msec
%     '->' separates events in the sequence
%     '<-' same as about but instead of the first events following it
%          specified if those events come before the others
%     '=' precedes the new event code
%     ':' is used to separate start and end times (in msec) of time window
%     Lists of event codes are placed in square brackets [], separated by commas
%     Square brackets are also required when an event code follows a time window
%
% Output:
%   epoch_data - in the standard TimeSurfer format.
%
% See also: ts_recode_events
%
% Created:  08/20/07 by Rajan Patel
% Last Mod: 03/23/11 by Don Hagler
%

%% todo: use mmil_args2parms

epoch_data = [];
T_mri2head = [];
points_sf  = 0.01;

if ~iscell(data_files), data_files = {data_files}; end;

for d=1:length(data_files)
    data_file = data_files{d};
    if ~exist(data_file,'file');
        fprintf('%s: file not found\n',data_file);
        return;
    end;
end;

try
    options = varargin;
    for index = 1:length(options)
        if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
    end;
    if ~isempty( varargin ), opt=struct(options{:});
    else opt = []; end;
catch
    fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
    return;
end;

try
    opt.badchanfile;
    badchanfile = opt.badchanfile;    
catch  
    badchanfile=[];
end;

if ~exist(badchanfile,'file')
    fprintf('No bad channels selected...\n');
    badchans = [];
else
    fid=fopen(badchanfile,'rt');
    j = 1;
    while (~feof(fid))
     badchans{j} = fgetl(fid);
     j=j+1;
    end
    fclose(fid);
end

change_names = 0;
if isfield(opt,'channamefile')
  if ~isempty(opt.channamefile)
    if exist(opt.channamefile,'file')
      [chn_path,chn_name] = fileparts(opt.channamefile);
      change_names = 1;
      cnfid = fopen(opt.channamefile,'rt');
      chan = 1;
      while ~feof(cnfid)
        chan_names{chan} = fgetl(cnfid);
        chan = chan + 1;
      end
      fclose(cnfid);
      channel_log = fullfile(chn_path,sprintf('%s_log.txt',chn_name));
      cnlog = fopen (channel_log,'w+');
      fprintf(cnlog,'      \t Old \t New\n');
      fprintf(cnlog,'Sensor\tLabel\tLabel\n');
    else
      fprintf('ERROR: Channel name file not found: %s\n',channamefile);
      return;
    end
  end
end

if isfield(opt,'noise_start') && ~isempty(opt.noise_start)
  noise_start = opt.noise_start;
else
  noise_start = -.08;
end

if isfield(opt,'noise_end') && ~isempty(opt.noise_end)
  noise_end = opt.noise_end;
else
  noise_end = 0;
end

if isfield(opt,'recode_rules') && ~isempty(opt.recode_rules)
    if ~iscell(opt.recode_rules), opt.recode_rules = {opt.recode_rules}; end
    recode_rules_obj = parse_recode_rules(opt.recode_rules);
end

if isfield(opt,'timelimits') && ~isempty(opt.timelimits)
    timelimits = opt.timelimits;
else
    timelimits = [-2 2];
end

%%
epoch_data = [];
epoch_data.num_sensors = 0;
epoch_data.sfreq = 0;
epoch_data.sensor_info = [];
epoch_data.coor_trans = [];
epoch_data.epochs=[];
epoch_data.noise.num_trials = 0;
epoch_data.noise.num_samples = 0;
epoch_data.noise.covar = [];

for d=1:length(data_files)
    data_file = data_files{d};
    [fpath,fname,filetype,dum] = fileparts(data_files{d});
    % find the appropriate EEG Lab function to call
    switch filetype
        case '.eeg'
          data = pop_loadeeg(data_file);    
        case '.vhdr'
          data = pop_loadbv(fpath,[fname filetype]);
				case '.set'
	  			data = pop_loadset('filename',[fname filetype],'filepath',fpath);
      case '.edf'
          data = pop_biosig([fname filetype]);
%           data = readedf(data_file);
        otherwise
          error('%s: File format not recognized: %s.',mfilename,data_file);
    end
    if d==1
      epoch_data.num_sensors = data.nbchan;
      epoch_data.sfreq = data.srate;
      for i=1:data.nbchan
        if change_names                                                           % was a list of new chan names given
          if length(chan_names)==epoch_data.num_sensors                           % make sure there are enough names
            epoch_data.sensor_info(i).label = chan_names{i};
            fprintf(cnlog,'%-6s\t%-5s\t%-5s\n',num2str(i),data.chanlocs(i).labels,chan_names{i}); % print to log file
          else
            error('ERROR: Number of new channel names does not match number of original channels.\n');
          end
        else
          epoch_data.sensor_info(i).label = data.chanlocs(i).labels;
        end
        %% Some default sensor values to make TimeSurfer happy
        epoch_data.sensor_info(i).typestring = 'eeg';
        epoch_data.sensor_info(i).type = 1;
        epoch_data.sensor_info(i).kind = 2;
        if strmatch(epoch_data.sensor_info(i).label,badchans,'exact')
          epoch_data.sensor_info(i).badchan = 1;
        else
          epoch_data.sensor_info(i).badchan = 0;
        end
        epoch_data.sensor_info(i).lognum = i;
        epoch_data.sensor_info(i).loc    = eye(4);
      end; % sensors loop
        epoch_data.coor_trans.device2head = eye(4);
        epoch_data.coor_trans.mri2head = T_mri2head;
        epoch_data.noise.covar         = zeros(epoch_data.num_sensors,epoch_data.num_sensors);
    end;
    
    if strcmp(filetype,'.set')
      %% EEGLAB .set file
      conds=[data.epoch.eventtype];
    elseif isfield(data.event,'epochtype')
      %% Neuroscan .eeg file
      conds=[data.event.epochtype];
    elseif isfield(data.event,'type')
      %% Brain Vision files
      tmp  = regexp({data.event.type},'\d*','match');
      keep = ~cellfun('isempty',tmp);
      tmp  = tmp(keep);
      data.event   = data.event(keep);
      data.urevent = data.urevent(keep);
      for i = 1:length(tmp)
          if ~isempty(tmp{i})
              conds(i)=mmil_cell2num(tmp{i});
          else
              conds(i)=nan;
          end
      end
      clear tmp
    end
    %% if it is continuous data - make epochs based on timelimits
    if data.trials == 1
      fprintf('%s: Epoching the data into trials with latency of %0.4f to %0.4f seconds.\n',mfilename,timelimits(1),timelimits(end));
      data.trials = 0;
      tmpdata     = data.data;
      data.data   = [];
      % epoch the data
      progress('init','textbar');
      for i = 1:length(conds)
        progress(i/length(conds));
        if ~isnan(conds(i))
            data.trials                = data.trials + 1;
            start_pnt                  = max(data.event(i).latency+(timelimits(1)*data.srate),1);
            end_pnt                    = min(data.event(i).latency+(timelimits(end)*data.srate),data.pnts);
            data.data(:,:,data.trials) = tmpdata(:,start_pnt:end_pnt);
        end
      end
      progress('close')
      % remove the empty conditions
      data.event  = data.event(~isnan(conds)); 
      % change the data to epoched format
      data.pnts   = size(data.data,2);
      data.xmin   = timelimits(1);
      data.xmax   = timelimits(end);
      fprintf('%s: Done.\n',mfilename);
    end
    % remove nonsense conditions
    conds = conds(~isnan(conds));
    % this will be our event codes
    types = unique(conds);
    badchannels = find([epoch_data.sensor_info.badchan]==1);
    %% if there is no accept information make all trials good
    if ~isfield(data.event,'accept');
        [data.event.accept] = deal(1);
    end
    fprintf('Event Codes: %s.\n',num2str(types));
    curr_epoch = 1;    
    try
      if length(conds) > data.trials
        [jnk xi xj] = unique([data.event.epoch],'first');    
        conds = conds(xi);
        try data.event = data.event(xi); end
        try data.urevent = data.urevent(xi); end
%         if issubfield(data,'event.accept')
%           acc = [data.event(xi).accept];
%           conds = conds(acc==1);
%         end        
      end
    end
    for i=1:length(types)
       index=find(conds==types(i));
       % only add a condition if it has some good trials
       if ~isempty(index) && ~isempty(find([data.event(index).accept],1))
        epoch_data.epochs(curr_epoch).event_code = types(i); 
        epoch_data.epochs(curr_epoch).num_trials = 0;
        epoch_data.epochs(curr_epoch).num_rejects.mag = 0;
        epoch_data.epochs(curr_epoch).num_rejects.grad= 0;
        epoch_data.epochs(curr_epoch).num_rejects.eeg= 0;
        epoch_data.epochs(curr_epoch).num_rejects.eog= 0;
        epoch_data.epochs(curr_epoch).num_rejects.manual= 0;
        epoch_data.epochs(curr_epoch).num_rejects.skip = 0;
        time = [data.xmin:(1/data.srate):data.xmax];
        epoch_data.epochs(curr_epoch).time = time(1:data.pnts);
        epoch_data.epochs(curr_epoch).data = [];
        noise_points = find(epoch_data.epochs(curr_epoch).time >= noise_start & epoch_data.epochs(curr_epoch).time <= noise_end);
        for j=1:length(index)
           if data.event(index(j)).accept            
            epoch_data.epochs(curr_epoch).data       = cat(3,epoch_data.epochs(curr_epoch).data,data.data(:,:,index(j)));
            epoch_data.epochs(curr_epoch).num_trials = epoch_data.epochs(curr_epoch).num_trials + 1;         
            noise = squeeze(data.data(:,noise_points,index(j)));
            num_samples = size(noise,2);
            noise = noise - mean(noise,2)*ones(1,num_samples);
            epoch_data.noise.covar = epoch_data.noise.covar + cov(noise');
            epoch_data.noise.num_samples = epoch_data.noise.num_samples + num_samples;
            epoch_data.noise.num_trials  = epoch_data.noise.num_trials + 1;
           else
            epoch_data.epochs(curr_epoch).num_rejects.eeg=epoch_data.epochs(curr_epoch).num_rejects.eeg + 1;
           end
        end  
        % remove the data in the bad channels
        epoch_data.epochs(curr_epoch).data(badchannels,:,:) = nan;
        curr_epoch = curr_epoch+1;
       else
           fprintf('WARNING: There are no valid events for event type %s, it will be excluded in the data set.\n',num2str(types(i)));
       end
    end
    if exist('recode_rules_obj','var')
        for rr = 1:length(recode_rules_obj)
            %% THIS IS BASED ON CODE BEN WROTE FOR TIMESURFER
            %% IF THERE ARE ANY QUESTIONS IT IS GOOD TO LOOK AT HIS
            %% ORIGINAL CODE FOR RECODING EVENTS
            %% I HAVE ONLY ADDED THE DIRECTIONALITY OPTION AND LACK OF TIME
            %% OPTION
            fprintf('Applying recoding rule: %s.\n',opt.recode_rules{rr});
            index = [];            
            if strcmp(func2str(recode_rules_obj(rr).event_sequence(1).method),'match_by_epoch')
                fprintf('WARNING: match_by_epoch not supported by epoched eeg files because no reference to original times exists;\n');
                fprtinf('         only the sequence of events will be considered regardless of timing.\n');
            end
            condindex = find(ismember(conds,recode_rules_obj(rr).event_sequence(1).conditions));
            for ind = 1:length(condindex)
                switch recode_rules_obj(rr).event_sequence(1).direction
                    case 'following'
                        if (condindex(ind) > 1) && ...
                            ismember(conds(condindex(ind)-1),recode_rules_obj(rr).event_sequence(2).conditions)
                           index = cat(1,index,condindex(ind));
                        end
                    case 'preceding'
                        if (condindex(ind) < length(conds)) &&  ...
                             ismember(conds(condindex(ind)+1),recode_rules_obj(rr).event_sequence(2).conditions)
                           index = cat(1,index,condindex(ind));
                        end
                    otherwise
                       error('Did not recognize the direction.'); 
                end
            end      
            if ~isempty(index) && ~isempty(find([data.event(index).accept],1))
                fprintf('Found %d matching events.\n',length(index));
                epoch_data.epochs(curr_epoch).event_code = [recode_rules_obj(rr).event_sequence.assign];
                epoch_data.epochs(curr_epoch).num_trials = 0;
                epoch_data.epochs(curr_epoch).num_rejects.mag = 0;
                epoch_data.epochs(curr_epoch).num_rejects.grad= 0;
                epoch_data.epochs(curr_epoch).num_rejects.eeg= 0;
                epoch_data.epochs(curr_epoch).num_rejects.eog= 0;
                epoch_data.epochs(curr_epoch).num_rejects.manual= 0;
                epoch_data.epochs(curr_epoch).num_rejects.skip = 0;
                time = [data.xmin:(1/data.srate):data.xmax];
                epoch_data.epochs(curr_epoch).time = time(1:data.pnts);
                epoch_data.epochs(curr_epoch).data = [];
                for j=1:length(index)
                    if data.event(index(j)).accept
                        epoch_data.epochs(curr_epoch).data = cat(3,epoch_data.epochs(curr_epoch).data,data.data(:,:,index(j)));
                        epoch_data.epochs(curr_epoch).num_trials = epoch_data.epochs(curr_epoch).num_trials + 1;
                    else
                        epoch_data.epochs(curr_epoch).num_rejects.eeg=epoch_data.epochs(curr_epoch).num_rejects.eeg + 1;
                    end
                end
                epoch_data.epochs(curr_epoch).data(badchannels,:,:) = nan;
                curr_epoch = curr_epoch+1;
            else
               fprintf('WARNING: There are no events that match the recoding rule #%d requested.\n',rr); 
            end
        end %% recoding rules
    end %% if recoding needed
end

epoch_data.noise.covar       = epoch_data.noise.covar/epoch_data.noise.num_trials;


if change_names, fclose(cnlog); end

%%%%%%%% Sub Functions %%%%%%%%
function recode_rules_obj=parse_recode_rules(recode_rules_str)
%
% recode_rules_str: must be cell array
%
% TAKEN FROM BEN'S CODE

if (~iscell(recode_rules_str))
    error('Programming error: recode_rules_str must be a cell array.');
end;
recode_rules_obj = [];
% syntax allows us to use splitstr A LOT
for i=1:length(recode_rules_str)
    recode_rules_obj(i).event_sequence = [];

    % find match and assign portions of string.
    tmp = splitstr(recode_rules_str{i},'=');
    if (length(tmp) ~= 2)
        error('Parse error when trying to find ''='' for rule ''%s''',...
            recode_rules_str{i});
    end;
    match_str = tmp{1};
    assign_str = tmp{2};
    %find transition sequence
    if strfind(match_str,'->')
      match_rules = splitstr(match_str,'->');
      direction   = 'following';
    elseif strfind(match_str,'<-')
      match_rules = splitstr(match_str,'<-');
      direction   = 'preceding';
    end
    %% todo: allow loose sequences that allow for intervening events
    %%       specified by "," rather than "->"
    for j=1:length(match_rules)
        % rule checking is in order from lowest to highest priority,
        % so that reassignment of the 'method' ends with the assignment
        % of the highest priority method.
        cur_rule = [];
        % loop through all match rules; stop when we've found an applicable one
        tmp = {}; iRule = 1;
        while (isempty(tmp))
            switch (iRule)
                % find events
                case 1
                    tmp=regexp(match_rules{j}, '^\[?([0-9,\s]+)\]?$', 'tokens');
                    if (~isempty(tmp))
                        cur_rule.conditions = double(mmil_cell2num(splitstr(tmp{1}{1},',')));
                        if j==1
                            cur_rule.mode = 'loose'; % does not have to be very next event
                        else
                            cur_rule.mode = 'strict';
                        end;
                        cur_rule.method = @match_by_sequence;
                        cur_rule.epoch = [];
                    end;
                    % find epoch
                case 2
                    if 0 % this produces the strange behavior of requiring an intervening
                        % event within the time window, but not caring which event
                        %  according to Ben's scheme, that behavior should come from a rule like this:
                        %    1->200:1500[*]=4
                        %  and 1->200:1500->[3]=4 should require no intervening events for that
                        %    time window
                        tmp=regexp(match_rules{j}, '^([0-9]*):([0-9]*)$', 'tokens');
                        if (~isempty(tmp))
                            cur_rule.conditions = [];
                            cur_rule.mode = 'relative';
                            cur_rule.method = @match_by_epoch;
                            cur_rule.epoch = [str2double(tmp{1}{1}) str2double(tmp{1}{2})];
                        end;
                    end;
                    % find epoch AND events
                case 3
                    tmp=regexp(match_rules{j},'^([0-9.]*):([0-9.]*)\[([0-9,\s]+)\]$', 'tokens');
                    if (~isempty(tmp))
                        cur_rule.conditions = mmil_cell2num(splitstr(tmp{1}{3:end},','));
                        cur_rule.mode = 'relative';
                        cur_rule.method = @match_by_epoch;
                        cur_rule.epoch = [str2double(tmp{1}{1}) str2double(tmp{1}{2})];
                    end;
                otherwise
                    tmp = 'failure!';
                    error('failed to parse rule ''%s''.',recode_rules_str{i});
            end;
            iRule = iRule+1;
        end; % loop over match rules
        clear('tmp','iRule');
        % recode only first in sequence
        cur_rule.direction  = direction;
        if j==1
            cur_rule.assign = str2num(assign_str);
        else
            cur_rule.assign = [];
        end;
        %% todo: how to recode only one event but allow it to be other than first?
        if (isempty(recode_rules_obj(i).event_sequence))
            recode_rules_obj(i).event_sequence = cur_rule;
        else
            recode_rules_obj(i).event_sequence(end+1) = cur_rule;
        end;
    end;
end;
