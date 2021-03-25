function cont_data = ts_loadfif(datafile,selectflg,datafield)
if nargin < 2, selectflg = 0;        end
if nargin < 3, datafield = 'epochs'; end

hdr = ts_read_fif_header(datafile,0);
if ~(isequal(selectflg,0) || isequal(selectflg,[]))
  hdr = select_sensors(hdr,selectflg);
end
[data,hdr.nBuffs,hdr.tlast,events,sfreq] = ts_read_fif_chan(datafile,hdr.sensors.label);

cont_data                                = [];
cont_data.num_sensors                    = length(hdr.sensors.label);
cont_data.sensor_info                    = [];
cont_data.coor_trans                     = [];
cont_data.noise                          = [];
cont_data.(datafield)                    = [];
cont_data.sfreq                          = sfreq;
cont_data.coor_trans.device2head         = loadtrans(datafile);
cont_data.coor_trans.mri2head            = [];
cont_data.noise.num_trials               = 0;
cont_data.noise.num_samples              = 0;
cont_data.noise.covar                    = [];
cont_data.noise.covar                    = [];
cont_data.(datafield).event_code         = 1; 
cont_data.(datafield).num_trials         = 1;
cont_data.(datafield).num_rejects.mag    = 0;
cont_data.(datafield).num_rejects.grad   = 0;
cont_data.(datafield).num_rejects.eeg    = 0;
cont_data.(datafield).num_rejects.eog    = 0;
cont_data.(datafield).num_rejects.manual = 0;
cont_data.(datafield).num_rejects.skip   = 0;
cont_data.(datafield).time               = [0:size(data,2)-1]/hdr.sfreq + hdr.tfirst;
cont_data.(datafield).data               = single(data);
clear data
for k = 1:cont_data.num_sensors
  cont_data.sensor_info(k).label         = hdr.sensors.label{k};
  cont_data.sensor_info(k).typestring    = hdr.sensors.typestring{k};
  cont_data.sensor_info(k).type          = hdr.sensors.type(k);
  cont_data.sensor_info(k).kind          = hdr.sensors.kind(k);
  cont_data.sensor_info(k).badchan       = 0;
  cont_data.sensor_info(k).lognum        = hdr.sensors.lognum(k);
  cont_data.sensor_info(k).loc           = hdr.sensors.loc{k};
%   try cont_data.sensor_info(k).range     = hdr.sensors.range(k); end
%   try cont_data.sensor_info(k).cal       = hdr.sensors.cal(k);   end
end

function hdr = select_sensors(hdr,typeflag)
types = unique(hdr.sensors.typestring);
if ischar(typeflag) && ~isempty(find(ismember(types,typeflag)))  
  typeflag = find(ismember(types,typeflag));
  types    = types(typeflag);
elseif iscell(typeflag)
  [sel1,sel2] = match_str(types,typeflag);
  typeflag = sel1;
  types    = types(typeflag);
elseif typeflag == 1
  fprintf('Available sensor types:\n');
  for k = 1:length(types)
    fprintf('  [%g] %s\n',k,types{k});
  end
  r     = input('which sensors do you want to load (enclose numbers in []): ');
  types = types(r);
else
  types = types(typeflag);
end
sens  = [];
n     = 0; 
for k = 1:length(hdr.sensors.typestring)
  if ismember(hdr.sensors.typestring(k),types)
    n = n + 1;
    sens.label{n}       = hdr.sensors.label{k};
    sens.kind(n)        = hdr.sensors.kind(k);
    sens.lognum(n)      = hdr.sensors.lognum(k);
    sens.type(n)        = hdr.sensors.type(k);
    sens.loc{n}         = hdr.sensors.loc{k};
    sens.typestring{n}  = hdr.sensors.typestring{k};
    sens.range(n)       = hdr.sensors.range(k);
    sens.cal(n)         = hdr.sensors.cal(k);
  end
end
hdr.sensors = sens;
hdr.nChans  = n;