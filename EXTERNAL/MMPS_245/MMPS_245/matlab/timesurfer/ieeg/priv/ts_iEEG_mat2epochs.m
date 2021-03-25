function [epoch_data]=ts_iEEG_mat2epochs(filename,varargin)
%ts_iEEG_mat2epochs imports .mat file to time surfer format epoch_data
%                   structure
%
%   E = ts_iEEG_mat2epochs(FILENAME) imports data with default values
%
%   E = ts_iEEG_mat2epochs(FILENAME,DIM_ORDER) specify the correct dimension
%   order.  DIM_ORDER is a 3 element row vector where each value specifies
%   which dimension is channels, samples, and trials RESPECTIVELY. default
%   is [3 2 1]
%
%   E = ts_iEEG_mat2epochs(FILENAME, SFREQ) specify the sampling frequency.
%   default is 1000 Hz.
%
%   E = ts_iEEG_mat2epochs(FILENAME, TIME_ON_OFF) specify the time in
%   seconds for the time vector.  TIME_ON_OFF is a 2 element row vector
%   where TIME_ON_OFF(1) is the onset and TIME_ON_OFF(2) is the end of the
%   epoch time window.  The size of the data matrice  is used with the vector 
%   to create the time vector.  This may result in a conflict with the SFREQ,
%   in which a warning is issued. default is [-3.5 2.5].
%
%   NOTE: Order of opitional arguments does not matter, but the filename
%   must be the first argument.
%   
%
%   For example, to import data with a sampling rate of 2000 Hz and a start
%   time of 0 and end time of 3:
%       filename='TT2_S100ALL.mat';
%       time=[0 3];
%       freq=2000;
%       epoch_data = ts_iEEG_mat2epochs(filename,time,freq)
%
%   .mat file should contain
%
%

%PRINT HELP IF NO ARGS PASSED
if nargin==0,
    help ts_iEEG_mat2epochs
    return
end

%DEFAULT VARARG VALUES
dim_order=[3 2 1];
time_on_off=[-3.5 2.5];
sfreq=1000;

%PARSE VARARG FOR NEW VALUES
if nargin>1,
    for i=1,size(varargin,2),
        if isnumeric(varargin{i}) && isvector(varargin{i}),
            if numel(varargin{i}) == 3  dim_order=varargin{i};      end;
            if numel(varargin{i}) == 2  time_on_off=varargin{i};    end;
        end
        if isnumeric(varargin{i}) && ~isvector(varargin{i}),
            sfreq=varagin{i};
        end
    end
end

if exist(filename)~=0,
    data=load(filename);
else
    error('ERROR: File specified in filename cannot be found. Check your path.')
end

conditions=fields(data);
epoch_data=[];

%SET "GLOBAL" PARAMETERS OF THE STUDY
%epoch_data
%            num_sensors ??? the number of sensors
sz_cond=size(getfield(data,conditions{1}));             %ASSUME NUM CHANNELS SAME FOR ALL CONDS!!!!
epoch_data.num_sensors=sz_cond(dim_order(1));
%This Lbl thing below just gives a name to the channels, update code later
%to be flexible to the number of channels and some cfg prefix.
Lbl=strcat(repmat('CH',epoch_data.num_sensors,1),num2str([1:epoch_data.num_sensors]','%02d'));

%            sensor_info - this is the length of the number of sensors in the data set
%                        label - the name of the channel
%                        typestring - all labels for intracranial data
%                        type - a number for the type (all 1 for intracranial)
%                        kind - a number for the kind (all 2 for intracranial)
%                        badchan - a flag to mark the channel as bad
%                        lognum - the sensor number
%                        loc - a location (4x4) matrix (all identity for intracranial)
%badchan=zeros(epoch_data.num_sensors,1);  %ASSUMING NO BAD CHANNELS
for i=1:epoch_data.num_sensors,
    epoch_data.sensor_info(i).label=Lbl(i,:);
    epoch_data.sensor_info(i).typestring='eeg';
    epoch_data.sensor_info(i).type=1;
    epoch_data.sensor_info(i).kind=2;
    epoch_data.sensor_info(i).badchan=0;%badchan(i);
    epoch_data.sensor_info(i).lognum=i;
    epoch_data.sensor_info(i).loc=eye(4);
end
%            coor_trans
%                        device2head - 4x4 matrix (all identity for intracranial)
%                        mri2head  - empty
epoch_data.coor_trans.device2head=eye(4);
epoch_data.coor_trans.mri2head=[];
%            sfreq ??? sampling frequency in Hz
epoch_data.sfreq=sfreq;
%            noise
%                        num_trials - the number of trials included in the noise calculation
%                        num_samples - number of samples included in the noise calculation
%                        cover - noise covariance matrix (num_chans x
%                        num_chans)
epoch_data.noise.num_trials=0;
epoch_data.noise.num_samples=0;
epoch_data.noise.covar=[];
%            epochs - this is the length of the number of conditions in the data set
%                        event_code - the event code for this conditions
%                        num_trials - the number of trials
%                        num_rejects - number of  trials rejected based on the following chans:
%                                    mag - 0
%                                    grad - 0
%                                    eog - 0
%                                    eeg
%                                    manual
%                                    skip
%                        time - the time vector for this condition in seconds
%                        data - the data [num_channels x num_samples x
%                        num_trials]
%LOOP THROUGH EACH CONDITION AND APPEND THE DATA TO THE STRUCTURE
for c=1:size(conditions,1);
    display(['Importing condition #' num2str(c,2)])
    %epoch_data.epochs(c).event_code=conditions{c};
    epoch_data.epochs(c).event_code=c;
    epoch_data.epochs(c).num_trials=sz_cond(dim_order(3));
    %num_rejects
    epoch_data.epochs(c).num_rejects.mag=0;
    epoch_data.epochs(c).num_rejects.grad=0;
    epoch_data.epochs(c).num_rejects.eog=0;
    epoch_data.epochs(c).num_rejects.eeg=0;
    epoch_data.epochs(c).num_rejects.manual=0;
    epoch_data.epochs(c).num_rejects.skip=0;
    %this code below can get messy with some data, i was given data from
    %-3.5 to 2.5 with only 6000 samples (there should be 6001) update this part of the code in the future;
    epoch_data.epochs(c).time=[time_on_off(1)+(1/sfreq):(1/sfreq):time_on_off(2)];
    %CHECK AND WARN IF TIME VECTOR AND SFREQ ARE NOT IN AGREEMENT
    if (numel(epoch_data.epochs(c).time)~=sz_cond(dim_order(2)))
        warning('WARNING: conflicting sampling information! SFREQ does not agree with TIME_ON_OFF w/ regards to size of time dimension of data')
        display(['numel(time) = ' num2str(numel(epoch_data.epochs(c).time)) '  numel(data_sz_time) = ' num2str(sz_cond(dim_order(2)))]);
    end
    epoch_data.epochs(c).data=permute(getfield(data,conditions{c}),dim_order);
end
