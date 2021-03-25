function avg_data = ts_import_neuroscan2avg(data_files,points_file,T_mri2head,points_sf);
%function avg_data = ts_import_neuroscan2avg(data_files,points_file,[T_mri2head],[points_sf]);
%
% Required Input:
%   data_files: full path file name (or cell array of file names) of neuroscan
%     format avg file
%   points_file: full path file name of neuroscan (polhemus) point file
%     this should be an ascii file converted from the .3dd format file
%
% Optional Input:
%   T_mri2head: 4x4 transformation matrix specifying registration between MRI
%     and measurements
%     {default = []}
%   points_sf: scaling factor to put points coordinates in meters
%     {default = 0.01}
%
% Created:  04/11/07 by Don Hagler
% Last Mod: 04/11/07 by Don Hagler
%

avg_data = [];

if nargin<2, help(mfilename); return; end;

if ~exist('T_mri2head','var'), T_mri2head = []; end;
if ~exist('points_sf','var') | isempty(points_sf), points_sf=0.01; end;

if ~iscell(data_files), data_files = {data_files}; end;

if ~exist(points_file,'file')
  fprintf('%s: file %s not found\n',mfilename);
  return;
end;

for d=1:length(data_files)
  data_file = data_files{d};
  if ~exist(data_file,'file');
    fprintf('%s: file %s not found\n',mfilename);
    return;
  end;
end;

points = elec_load_scan3dasc(points_file);
eeg_points = [];
for p=1:length(points.x)
  eeg_points(p).coords = points_sf*[points.x(p),points.y(p),points.z(p)];
end;
ref_coords = points_sf*points.ref(1:3);

avg_data = [];
avg_data.sensor_info = [];
avg_data.coor_trans = [];
avg_data.averages = [];
avg_data.noise = [];

for d=1:length(data_files)
  data_file = data_files{d};
  [signal, variance, chan_names, ...
     pnts, rate, xmin, xmax] = loadavg(data_file);
  if d==1
    chan_names = char(chan_names);
    nchans = size(signal,1);
    npoints = size(signal,2);
    sfreq = rate;
    avg_data.num_sensors = nchans;
    avg_data.sfreq = sfreq;
    for i=1:nchans
      avg_data.sensor_info(i).label = deblank(chan_names(i,:));
      if findstr(lower(avg_data.sensor_info(i).label),'eog')
        avg_data.sensor_info(i).typestring = 'eog';
        avg_data.sensor_info(i).type = 5;
        avg_data.sensor_info(i).kind = 202;
      else
        avg_data.sensor_info(i).typestring = 'eeg';
        avg_data.sensor_info(i).type = 1;
        avg_data.sensor_info(i).kind = 2;
      end;
      avg_data.sensor_info(i).badchan = 0;
      avg_data.sensor_info(i).lognum = i;
      avg_data.sensor_info(i).loc = eye(4);
      avg_data.sensor_info(i).loc(1:3,4) = eeg_points(i).coords;
      avg_data.sensor_info(i).loc(1:3,1) = ref_coords;
    end;
    avg_data.coor_trans.device2head = eye(4);
    avg_data.coor_trans.mri2head = T_mri2head;
  end;
  avg_data.averages(d).event_code = d;
  avg_data.averages(d).num_trials = 100;
  avg_data.averages(d).num_rejects = 0;
  time = [xmin:(1/rate):xmax];
  avg_data.averages(d).time = time(1:npoints);
  avg_data.averages(d).data = signal;
  avg_data.averages(d).stdev = sqrt(variance);
end;

