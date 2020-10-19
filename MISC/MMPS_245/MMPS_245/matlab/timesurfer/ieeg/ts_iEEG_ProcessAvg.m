function ts_iEEG_ProcessAvg (varargin)

% Use: ts_iEEG_ProcessAvg ('avgfile',avg_files,...);
%
% Used to convert Neuroscan .avg files into a avg_data TimeSurfer
% format with all averages suitable for further use in the analysis stream.
%
% It is best to input all the averages for a necessary subject at this
% time.
%
% Required Input:
%
%  <None> - User will be prompted for all necessary inputs
%
% Optional Input:
%
%  avgpath     - directory containing avgfiles.
%  avgfiles    - cell array of ALL averages to be put into the TimeSurfer
%                avg_data structure to be saved.
%  eventcodes  - vector of event codes to use for each file given
%  savefile    - path/name of file to save avg_data structure to.
%  badchanfile - path/name of text file containing bad channels.
%
% Created:  09/20/07 by Rajan Patel
% Rcnt Mod: 09/20/07 by Rajan Patel
% Last Mod: 09/15/12 by Don Hagler

%% Check Options

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

% set up input file
if isfield(opt,'avgfiles')
    avg_files = opt.avgfiles;
else
    [avg_files,avgpath,index]=uigetfile({'*.mat'},'Select data file','MultiSelect','on');
end
if ~iscell(avg_files), avg_files = {avg_files}; end

if isfield(opt,'avgpath')
  avgpath=opt.avgpath;
else
  avgpath = [];
end

for i = 1:length(avg_files)
  avg_files{i} = fullfile(avgpath,avg_files{i});
end

if isfield(opt,'eventcodes')
  eventcodes = opt.eventcodes;
  if length(eventcodes) ~= length(avg_files)
    fprintf('Number of event codes supplied does not match number of avg files.\n');
    return;
  end
else
  eventcodes = [];
end

if isfield(opt,'savefile'),  f_name = opt.savefile;
else
    f_name = avg_files{1};
    f_name = [f_name(1:strfind(f_name,'.cond')) '.avg.mat'];
end

[f_path, tmp_name, ext] = fileparts(f_name);
if ~exist(f_path,'dir')
    fprintf('Creating directory: %s\n',f_path);
    mkdir(f_path);
end

if isfield(opt,'badchanfile')
    badchanfile=opt.badchanfile;
    if ~exist(badchanfile,'file'), fprintf('Unable to locate bad channel file: %s\n',badchanfile); end;
else
    badchanfile=[];
end

%% Process Data
[subject, date, task, condition, data_info, data_procs] = ts_iEEG_decon_fname (avg_files);         % Find condition number in file name.


avg_data = ts_iEEG_neuroscan2avg(avg_files,'badchanfile',badchanfile);

% set event_code
for i=1:length(avg_files);
  try  
    condition_num = strtok(char(condition(i)),'cond')                                              
  catch
    condition_num = i;
  end
    if isempty(eventcodes)
      avg_data.averages(i).event_code=str2num(condition_num);                                            % Use condition number as event_code
    else
      avg_data.averages(i).event_code=eventcodes(i);
    end
end

try
    save(f_name,'avg_data');
catch
    fprintf('ERROR saving data file\n');
end;
