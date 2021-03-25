function avg_data = ts_iEEG_neuroscan2avg(data_files,varargin);

% avg_data = ts_iEEG_neuroscan2avg(data_files,'badchanfile','badchanfile.txt');
%
% Required Input:
%   data_files: full path file name (or cell array of file names) of neuroscan
%     format avg file
%
% Optional Input:
% 
%   badchanfile:  path/filename of a txt file containing list of bad
%   eventcodes:   list of eventcodes to assign to conditions
%   channels
%
% Removed the following Inputs:
%   points_file: full path file name of neuroscan (polhemus) point file
%     this should be an ascii file converted from the .3dd format file
%   T_mri2head: 4x4 transformation matrix specifying registration between MRI
%     and measurements
%     {default = []}
%   points_sf: scaling factor to put points coordinates in meters
%     {default = 0.01}
%
% Created:  04/11/07 by Don Hagler
% Last Modified: 08/22/07 by Rajan Patel - made specific for iEEG data
% (ignore location information)
%

avg_data = [];

if nargin<1, help(mfilename); return; end;
 
T_mri2head = [];
points_sf=0.01;

if ~iscell(data_files), data_files = {data_files}; end;

% if ~exist(points_file,'file')
%   fprintf('%s: file %s not found\n',mfilename);
%   return;
% end;

for d=1:length(data_files)
  data_file = data_files{d};
  if ~exist(data_file,'file');
    fprintf('%s: file %s not found\n',mfilename,data_file);
    return;
  end;
  clear data_file
end;

% points = elec_load_scan3dasc(points_file);
% eeg_points = [];
% for p=1:length(points.x)
%   eeg_points(p).coords = points_sf*[points.x(p),points.y(p),points.z(p)];
% end;
% ref_coords = points_sf*points.ref(1:3);

options=varargin;
if ~isempty( varargin ), opt=struct(options{:}); end;

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

if exist('opt','var') && isfield(opt,'eventcodes')
  eventcodes = opt.eventcodes;
  if length(eventcodes) ~= length(data_files)
    fprintf('Number of eventcodes not consistent with number of files.\n');
  end
else
  eventcodes = [];
end
avg_data = [];
avg_data.sensor_info = [];
avg_data.coor_trans = [];
avg_data.averages = [];
avg_data.noise.num_trials = [];
avg_data.noise.num_samples = [];
avg_data.noise.covar = [];

for d=1:length(data_files)
  data_file = data_files{d};
  [signal, variance, chan_names, pnts, rate, xmin, xmax,trials] = ts_iEEG_loadavg(data_file);
  if d==1 % set up common fields
    chan_names = char(chan_names);
    nchans = size(signal,1);
    npoints = size(signal,2);
    sfreq = rate;
    avg_data.num_sensors = nchans;
    avg_data.sfreq = sfreq;
    for i=1:nchans
      avg_data.sensor_info(i).label = deblank(chan_names(i,:));
%       if findstr(lower(avg_data.sensor_info(i).label),'eog')
%         avg_data.sensor_info(i).typestring = 'eog';
%         avg_data.sensor_info(i).type = 5;
%         avg_data.sensor_info(i).kind = 202;
%       else
        avg_data.sensor_info(i).typestring = 'eeg';
        avg_data.sensor_info(i).type = 1;
        avg_data.sensor_info(i).kind = 2;
%       end;
      
      if strmatch(avg_data.sensor_info(i).label,badchans,'exact')
       avg_data.sensor_info(i).badchan = 1;
      else
       avg_data.sensor_info(i).badchan = 0;
      end
      avg_data.sensor_info(i).lognum = i;
      avg_data.sensor_info(i).loc = eye(4);
%     avg_data.sensor_info(i).loc(1:3,4) = eeg_points(i).coords;
%     avg_data.sensor_info(i).loc(1:3,1) = ref_coords;
    end;
    avg_data.coor_trans.device2head = eye(4);
    avg_data.coor_trans.mri2head = T_mri2head;
  end;
  if isempty(eventcodes)
    avg_data.averages(d).event_code = d;
  else
    avg_data.averages(d).event_code = eventcodes(d);
  end
  avg_data.averages(d).num_trials = trials;
  avg_data.averages(d).num_rejects.mag    = 0;
  avg_data.averages(d).num_rejects.grad   = 0;
  avg_data.averages(d).num_rejects.eog    = 0;
  avg_data.averages(d).num_rejects.eeg    = 0;
  avg_data.averages(d).num_rejects.manual = 0;
  avg_data.averages(d).num_rejects.skip   = 0;
  time = [xmin:(1/rate):xmax];
  avg_data.averages(d).time = time(1:npoints);
  avg_data.averages(d).data = signal;
  if isnumeric(variance)
    avg_data.averages(d).stdev = sqrt(variance);
  else
    avg_data.averages(d).stdev = [];
  end
end;

