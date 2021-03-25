function outdata = ts_collapse(indata)
%Purpose: Bins all events in a timesurfer epoch data structure into
%         a single condition. The event code of the new condition is equal
%         to the sum of it's component event codes. The original identities
%         of each trial are preserved in the epochs.trial_info.event_code field.
%         This function can be undone by using ts_uncollapse.
%
%Example: epoch_data = ts_collapse(epoch_data)
%
%Inputs: timesurfer epoch_data structure.
%Outputs: timesurfer epoch_data structure with all events bined into one condition.
%
%TODO: Make function work on timefreq data structures. 
%
%Created by BQR 05/18/12
%Modified by BQR 10/09/13 sort combined trials chronologically rather than by event code

[datatype,datafield,dataparam] = ts_object_info(indata);
events = [indata.(datafield).event_code];
%% sanity check
ti = 1;
if ~isfield(indata.(datafield),'trial_info')
    warning('ts_collapse: Input has no trial_info. Event identities may not be recoverable.');
    ti = 0;
end
if ti
    if ~isempty(setxor(...
        unique(cell2mat(arrayfun(@(x) [x.trial_info.event_code],indata.(datafield),'uni',0))),...
        events)) 
        warning(['ts_collapse: eventcodes in trial_info do not match those in epochs.event_code. '...
                 'Event identities may not be recoverable.']);
    end
end

%% collapse trials across condition
combinations = cellfun(@num2str,num2cell(events),'uni',0);
combinations = cellfun(@horzcat,combinations,repmat({'+'},1,length(events)),'uni',0);
combinations = [combinations{:}];
combinations = combinations(1:end-1);
outdata = ts_combine_conditions(indata,'combinations',{combinations},'neweventcodes',sum(events),'verbose',0);
outdata = ts_data_selection(outdata,'events',sum(events),'keepbadchans_flag',1,'markbadchans_flag',1);

%% sort trials chronologically
[~,sortidxs] = sort(outdata.(datafield).trial_info.latency);
fn = fieldnames(outdata.(datafield).trial_info);
for ifn = 1:numel(fn)
    outdata.(datafield).trial_info.(fn{ifn}) = ...
        outdata.(datafield).trial_info.(fn{ifn})(sortidxs);
end
outdata.(datafield).data = outdata.(datafield).data(:,:,sortidxs);
end