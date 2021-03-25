function outdata = ts_uncollapse(indata)
%Purpose: Bins trials in a timesurfer epoch data structure into
%         conditions as indicated in their respective trial_info.event_code
%         This function reverses the effect of ts_collapse.
%
%Example: indata = ts_uncollapse(indata)
%
%Inputs: timesurfer indata structure.
%Outputs: timesurfer indata structure with events sorted into bins of their original event codes.
%
%TODO: Make function work on timefreq data structures. 
%
%Created by BQR 05/18/12


%% sanity check
if ~isfield(indata.epochs,'trial_info')
    error('ts_collapse: Input has no trial_info. trial_info required.');
end

%% uncollapse trials across condition
events = unique(indata.epochs(1).trial_info.event_code);
for ievent = 1:length(events)
    indata.epochs(1+ievent).event_code = events(ievent);
    indata.epochs(1+ievent).num_trials = sum(indata.epochs(1).trial_info.event_code == events(ievent));
    indata.epochs(1+ievent).num_rejects = indata.epochs(1).num_rejects;
    indata.epochs(1+ievent).time = indata.epochs(1).time;
    fn = fieldnames(indata.epochs(1).trial_info);
    for ifn = 1:length(fn)
        indata.epochs(1+ievent).trial_info.(fn{ifn}) = indata.epochs(1).trial_info.(fn{ifn})(indata.epochs(1).trial_info.event_code == events(ievent));
    end
    indata.epochs(1+ievent).data = indata.epochs(1).data(:,:,indata.epochs(1).trial_info.event_code == events(ievent));
end
outdata = ts_data_selection(indata,'events',events,'keepbadchans_flag',1,'markbadchans_flag',1);
end
