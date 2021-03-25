function [epoch_data] = ts_loadcnt (data_files,varargin)
%function [epoch_data] = ts_loadcnt (data_files,varargin)
%
% Use [epoch_data] = ts_loadcnt (data_files,'badchanfile','badchanfile.txt','channamefile','chan_names.txt');
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
% Output:
%   epoch_data - in the standard TimeSurfer format.
%
% Created:  02/20/09 by Jason Sherfey
% Last Mod: 06/13/14 by Don Hagler
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
    [fpath,fname,filetype] = fileparts(data_files{d});
    % find the appropriate EEG Lab function to call
    switch filetype
        case '.cnt'
          data = pop_loadcnt(data_file);
        case '.eeg'
          data = pop_loadeeg(data_file);    
        case '.vhdr'
          data = pop_loadbv(fpath,[fname filetype]);
				case '.set'
	  			data = pop_loadset('filename',[fname filetype],'filepath',fpath);
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
      tmp = regexp({data.event.type},'\d*','match');
      for i = 1:length(tmp)
          if ~isempty(tmp{i})
              conds(i)=mmil_cell2num(tmp{i});
          else
              conds(i)=nan;
          end
      end
      clear tmp
    end
    %%
    %% ALL CODE ABOVE & MOST BELOW COPIED FROM ts_iEEG_eeg2epoch()
    %%
    %% if it isn't continuous data - return
    if data.trials ~= 1
      fprintf('%s: data is not continuous. Use ts_iEEG_eeg2epoch().\n',mfilename);
      return;
    end
    badchannels = find([epoch_data.sensor_info.badchan]==1);
    epoch_data.epochs.event_code = 1;
    curr_epoch = 1;    
     index     = 1;
      epoch_data.epochs(curr_epoch).event_code        = 1; 
      epoch_data.epochs(curr_epoch).num_trials        = 1;
      epoch_data.epochs(curr_epoch).num_rejects.mag   = 0;
      epoch_data.epochs(curr_epoch).num_rejects.grad  = 0;
      epoch_data.epochs(curr_epoch).num_rejects.eeg   = 0;
      epoch_data.epochs(curr_epoch).num_rejects.eog   = 0;
      epoch_data.epochs(curr_epoch).num_rejects.manual= 0;
      epoch_data.epochs(curr_epoch).num_rejects.skip  = 0;
      time = [data.xmin:(1/data.srate):data.xmax];
      epoch_data.epochs(curr_epoch).time              = time(1:data.pnts);
      epoch_data.epochs(curr_epoch).data              = data.data;
      % remove the data in the bad channels
      epoch_data.epochs(curr_epoch).data(badchannels,:,:) = nan;
end
if change_names, fclose(cnlog); end
