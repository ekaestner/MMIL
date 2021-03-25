function outdata = ts_load_nsX(datafile,varargin)
% outdata = ts_load_nsX(datafile,varargin)
% Imports continuous .nsX data into timesurfer (matlab) format. 
% Optional parameters:
% nsx_dsfact : 1 : [] : downsampleing factor applied when loading 
%                   (more efficent than downsampling the output)
% precision: 'double' : {'double','short'}, precision of data at import,
%            using 'double (default) is strongly recommended
%
% Created: BQR 13.02.01
% Last modified: BQR 13.22.06

%% parse parameters
parms = mmil_args2parms(varargin, { ...
  'nsx_dsfact',1,[],...
  'precision','double',{'double','short'},...
 }, ...
 false );

%% nsX import
nsXstruct = opennsX(datafile,'read','report',...
    'precision',parms.precision,...
    'skipfactor',parms.nsx_dsfact);

%% initialize output
outdata = [];
outdata.num_sensors = 0;
outdata.sensor_info = [];
outdata.coor_trans = [];
outdata.sfreq = 0;
outdata.noise.num_trials = 0;
outdata.noise.num_samples = 0;
outdata.noise.covar = [];
outdata.epochs=[];

%% assemble metadata
% sensor data
outdata.num_sensors = nsXstruct.MetaTags.ChannelCount;
for ichan = 1:outdata.num_sensors
    outdata.sensor_info(ichan).label = num2str(nsXstruct.MetaTags.ChannelID(ichan));
    outdata.sensor_info(ichan).typestring = 'eeg';
    outdata.sensor_info(ichan).type = 1;
    outdata.sensor_info(ichan).kind = 2;
    outdata.sensor_info(ichan).badchan = 0;
    outdata.sensor_info(ichan).lognum = ichan;
    outdata.sensor_info(ichan).loc = eye(4);
end
% empty/dummy values for coor and noise
outdata.coor_trans.device2head = eye(4);
outdata.coor_trans.mri2head = [];
% sfreq
outdata.sfreq = nsXstruct.MetaTags.SamplingFreq / parms.nsx_dsfact;
%noise
outdata.noise.covar = zeros(outdata.num_sensors,outdata.num_sensors);

%% assemble data
outdata.epochs(1).event_code        = 1; 
outdata.epochs(1).num_trials        = 1;
outdata.epochs(1).num_rejects.mag   = 0;
outdata.epochs(1).num_rejects.grad  = 0;
outdata.epochs(1).num_rejects.eeg   = 0;
outdata.epochs(1).num_rejects.eog   = 0;
outdata.epochs(1).num_rejects.manual= 0;
outdata.epochs(1).num_rejects.skip  = 0;
outdata.epochs(1).time = linspace(0,size(nsXstruct.Data,2)*1/outdata.sfreq,size(nsXstruct.Data,2));
outdata.epochs(1).data = nsXstruct.Data;
end



