function ts_iEEG_makecondkey (f_name)

% ts_iEEG_condkey('/home/halgdev/projects/hh/NY07/cond_key.csv');
%
% Takes a text file containing corresponding condition numbers, event
% numbers and event names and makes a cond_key.mat file that can be used
% during plotting functions.
%
% Input: 
%
%   .csv file created from a spreadsheet where the first row contiains the
%   headers specifying the information in the column:
%
%       Event Code - the event code
%       Condition  - the condition number (place in TimeSurfer structure)
%       Name       - the name of the condition
%
%   Sample .csv file saved from Excel/Calc:
%
%    /space/monkeys/1/home/rpatel/matlab/scripts/ts_iEEG_stream/samples/condition_key_sample.csv
%
% Ouput:
%
%    cond_key.mat to the current directory containing cond_key data
%    structure with following fields:
%
%        event_code - the event code
%        cond       - the condition number
%        name       - the name of the condition
%
% Created on:       10/16/07 - Rajan Patel
% Last Modified on: 10/16/07 - Rajan Patel
%
% See also ts_iEEG_PlotAvg, ts_iEEG_PlotWavelets, ts_iEEG_Plot


[path,x] = fileparts(f_name);
data   = loadtxt (f_name,'delim',',','verbose','off');
header =  data(1,:);
data   =  data(2:end,:);

conditions = [];
events     = [];
names      = [];

for i = 1:length(header)                        % what do the columns mean
  header(i) = regexprep(header(i), '^"([^"]+)"$', '$1');
  if strcmpi(header(i),'Condition')
    conditions = i;
  elseif strcmpi(header(i),'Event Code')
    events     = i;
  elseif strcmpi(header(i),'Name')
    names      = i;
  end
end

%if isempty(conditions), fprintf('ERROR: There are no condition numbers, check the header line...\n'); return; end
if isempty(events),     fprintf('ERROR: There are no event codes, check the header line...\n');       return; end
if isempty(names),      fprintf('ERROR: There are no condition names, check the header line...\n');   return; end

for i=1:length(data)
  cond_key(i).event_code = cell2mat(data(i,events));
  if ~isempty(conditions), cond_key(i).cond       = cell2mat(data(i,conditions)); 
  else cond_key(i).cond = nan; end
  cond_key(i).name       = regexprep(data(i,names),'^"([^"]+)"$', '$1');
end

try
  save(fullfile(path,'cond_key.mat'),'cond_key');
catch
  fprintf('%s: Unable to save condition key file.\n',mfilename);
end
  