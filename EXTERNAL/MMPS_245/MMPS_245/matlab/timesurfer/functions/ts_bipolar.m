function outdata = ts_bipolar(indata)
%Purpose: Identifies sequences of depth electrodes by looking for a common
%         string component followed by a variable numeric component eg. 
%         {'AnTp1','AnTp2','AnTp3'}. Then a channel by channel subtraction
%         (or 1st order derivative) is performed, [X(2)-X(1)  X(3)-X(2) ...
%         X(n)-X(n-1)]. These data replace the previous data and the
%         channels are relabeled eg. {'AnTp2-1','AnTp3-2'}.Note that The
%         resulting data structure will have N less channels where N is the
%         number of depth electrode sequences. Note that only interger
%         electrode sequences numbers  0<x<10 are currently supported. 
%
%Example: epoch_data = ts_bipolar(epoch_data)
%
%Inputs:  timesurfer epoch_data structure or equivilent three diminsional
%         freqbandavg aka tfwave structure
%Outputs: The input data structure with depth data converted in to channel
%         by channel differences
%
%TODO: Make function work on timefreq_data and avg_data structures.
%TODO: Make fuction work with >1 digit electrode sequence numbers.      
%
%Created by BQR 06/07/12

%% Identify electrode sequences
labels = {indata.sensor_info.label};
numstridxs = regexp(labels, '[1-9]','start');
numstridxs = cellfun(@(x) (x-1),numstridxs,'uni',0);
changroups = labels;
for ichan = 1:numel(changroups)
    changroups{ichan} = changroups{ichan}(1:numstridxs{ichan});
end
unichangroups = unique(changroups);
grpidxs = cellfun(@(x) strmatch(x,changroups,'exact'),unichangroups,'uni',0);
grpidxs = grpidxs(~cellfun(@isempty,unichangroups)); % ignore chans without numeric component
grpidxs = grpidxs(cellfun(@numel,grpidxs)>1); % ignore single chan sequences
firstchansidxs = cellfun(@(x) x(1),grpidxs);
diffchanidxs = cell2mat(cellfun(@(x) x(2:end)',grpidxs,'uni',0));

%% Adjust data
for iepo = 1:numel(indata.epochs)
    for igrp = 1:numel(grpidxs)
        indata.epochs(iepo).data(grpidxs{igrp}(2:end),:,:) =...
            diff(indata.epochs(iepo).data(grpidxs{igrp},:,:),1,1);  
    end
    indata.epochs(iepo).data(firstchansidxs,:,:) = [];
end

%% Adjust header info
indata.num_sensors = indata.num_sensors - numel(firstchansidxs);
indata.noise = []; 
for ichan  = 1:numel(diffchanidxs)
    appendstr = sprintf('-%i',...
        str2double(indata.sensor_info(diffchanidxs(ichan)).label(end))-1);
    indata.sensor_info(diffchanidxs(ichan)).label = ...
        [indata.sensor_info(diffchanidxs(ichan)).label appendstr];
end
indata.sensor_info(firstchansidxs) = [];
outdata = indata;
end