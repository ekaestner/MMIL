%  epoch_data = ntools_gen_epoch_data(dat_file, dio_file, experiment, ['channels', wh_channels_to_keep], ...
%       ['max_epochs', max_epochs], ['is_continuous', 0|1], ['coords', coordinates])
%
%  Generate an epoch_data data structure from a experiment data_structure
%
%  Parameters:
%       dat_file: path to a .dat file
%       dio_file: path to a .dio.txt file
%       experiment: experiment structure
%
%  Optional keyword parameters:
%       channels: which channels to keep (default: all)
%       max_epochs: max num epochs (default: no limit)
%       is_continuous: treats all data in file as one long trial (default: no) - also see ntools_gen_cont_data 
%       coords: an nx3 array of electrode coordinates

function epoch_data = ntools_gen_epoch_data(dat_file, dio_file, experiment, varargin)
    [sfreq, num_channels, num_samples, precision, chan_names] = ntools_hdrxml_read(dat_file);
   
    max_epochs = Inf;
    wh_channels_to_keep = 1:num_channels;
    is_continuous = 0;
    
    for i = 1:length(varargin)
        if(ischar(varargin{i}))
            switch(lower(varargin{i}))
                case {'max_epochs', 'max epochs'}
                    max_epochs = varargin{i+1};
                case {'wh_channels_to_keep', 'channels'}
                    wh_channels_to_keep = varargin{i+1};
                case {'is_continuous', 'is_cont_data'}
                    is_continuous = varargin{i+1};
                case {'coords', 'coordinates', 'electrode_coords'}
                    coords = varargin{i+1};
            end
        end
    end        
    

    epoch_data = [];
    epoch_data = set_sensor_info(epoch_data, chan_names, wh_channels_to_keep);
    if(exist('coords', 'var'))
        epoch_data = ntools_gen_add_coords(epoch_data, coords, wh_channels_to_keep);
    end
    epoch_data = set_coor_trans(epoch_data);
    epoch_data = set_noise(epoch_data);
    epoch_data.sfreq = sfreq;
    epoch_data = set_epochs(epoch_data, dat_file, dio_file, experiment, sfreq, num_channels, ...
        num_samples, precision, max_epochs, wh_channels_to_keep, is_continuous);
end

function epoch_data = set_sensor_info(epoch_data, chan_names, wh_channels_to_keep)
    %NUM_NSPIKE_SENSORS = 256;

    epoch_data.num_sensors = length(wh_channels_to_keep);
    for i = wh_channels_to_keep
        epoch_data.sensor_info(i).label = chan_names{i};
        epoch_data.sensor_info(i).typestring = 'eeg';
        epoch_data.sensor_info(i).type = 1;
        epoch_data.sensor_info(i).kind = 2;
        epoch_data.sensor_info(i).badchan = 0;
        epoch_data.sensor_info(i).lognum = i;
        epoch_data.sensor_info(i).loc = eye(4);
    end
end

function epoch_data = set_coor_trans(epoch_data)
    epoch_data.coor_trans.device2head = eye(4);
    epoch_data.coor_trans.mri2head = [];
end

function epoch_data = set_noise(epoch_data)
    epoch_data.noise.num_trials = [];
    epoch_data.noise.num_samples = [];
    epoch_data.noise.covar = [];
end

function epoch_data = set_epochs(epoch_data, dat_file, dio_file, experiment, sfreq, num_channels, ...
        num_samples, precision, max_epochs, wh_channels_to_keep, is_continuous)
    NSPIKE_SAMPLING_RATE = 30000;  %  NSpike timestamps are in Nspike system sampling rate
    sfreq_to_nspike_sfreq_ratio = sfreq / NSPIKE_SAMPLING_RATE;
    samples_per_ms = sfreq / 1000;
    
    if(is_continuous)
        epoch_data.epochs(1).name = 'Whole recording';
        epoch_data.epochs(1).event_code = 0;
        epoch_data.epochs(1).num_trials = 1;
        epoch_data.epochs(1).num_rejects = struct('mag', 0, 'grad', 0, 'eeg', 0, 'eog', 0, 'manual', 0, 'skip', 0);

        epoch_data.epochs(1).time = linspace(0, num_samples-1, num_samples) / sfreq;
        data = ntools_load_samples(dat_file, precision, num_channels);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % downsample to conserve memory (silent if global not specified) -BQR 12.01.06
        global RAW_NSPIKE_DSFACTOR;
        if ~isempty(RAW_NSPIKE_DSFACTOR)
            dsfact = RAW_NSPIKE_DSFACTOR;
            data = downsample(data',dsfact)';
            epoch_data.epochs(1).time = downsample(epoch_data.epochs(1).time',dsfact)';
            epoch_data.sfreq = epoch_data.sfreq/dsfact;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        epoch_data.epochs(1).data = double(data(wh_channels_to_keep, :));
    else
        [time_codes, A, B, C, D] = ntools_dio_parse(dio_file); %#ok;

        for i = 1:length(experiment.event)
            experiment.event(i).timestamps = time_codes(ismember(C, experiment.event(i).source_code));
            
            epoch_data.epochs(i).name = experiment.event(i).name;
            epoch_data.epochs(i).event_code = experiment.event(i).source_code;
            epoch_data.epochs(i).num_trials = length(experiment.event(i).timestamps);
            epoch_data.epochs(i).num_rejects = struct('mag', 0, 'grad', 0, 'eeg', 0, 'eog', 0, 'manual', 0, 'skip', 0);

            % offsets in ms
            start_offset = experiment.event(i).start_offset;
            stop_offset = experiment.event(i).stop_offset;

            num_trials = min(length(experiment.event(i).timestamps), max_epochs);
            num_samples = round((stop_offset - start_offset) * samples_per_ms);
            epoch_data.epochs(i).time = linspace(start_offset, stop_offset, num_samples) / 1000;

            if(num_trials == 0)
                error('No trials/epochs found in recording. Use ntools_gen_cont_data for continuous data.');
            end
            
            % trial_info (added by Jason Sherfey on 10-Feb-2011)
            evt = experiment.event(i).source_code;
            lat = experiment.event(i).timestamps;
            epoch_data.epochs(i).trial_info=[];
            epoch_data.epochs(i).trial_info.number(1,:)                      = find(ismember(C,evt));
            epoch_data.epochs(i).trial_info.latency(1,:)                     = round(experiment.event(i).timestamps * sfreq_to_nspike_sfreq_ratio);
            epoch_data.epochs(i).trial_info.badtrial(1,:)                    = zeros(1,length(lat));
            epoch_data.epochs(i).trial_info.event_code(1,:)                  = evt*ones(1,length(lat));
            epoch_data.epochs(i).trial_info.duration(1,:)                    = ones(1,length(lat));
            [epoch_data.epochs(i).trial_info.datafile{1,1:length(lat)}]      = deal(dat_file);
            [epoch_data.epochs(i).trial_info.events_fnames{1,1:length(lat)}] = deal(dio_file);

            epoch_data.epochs(i).data = zeros(num_channels, num_samples, num_trials, 'double'); %double for flaws in Sig Proc Toolbox

            for j = 1:num_trials
                timestamp = experiment.event(i).timestamps(j); % in 30 kHz
                start = round(timestamp * sfreq_to_nspike_sfreq_ratio + start_offset * samples_per_ms) + 1;
                stop = round(timestamp * sfreq_to_nspike_sfreq_ratio + stop_offset * samples_per_ms);
                [data, eof] = ntools_load_samples(dat_file, precision, num_channels, start, stop);
                if(eof)
                    error('Trial %d with code %d extends beyond length of file', j, epoch_data.epochs(i).event_code);
                end
                epoch_data.epochs(i).data(:, :, j) = double(data(wh_channels_to_keep, :));
            end
        end
    end
    epoch_data.experiment = experiment;
    epoch_data.experiment.dat_file = dat_file;
    epoch_data.experiment.dio_file = dio_file;
end
