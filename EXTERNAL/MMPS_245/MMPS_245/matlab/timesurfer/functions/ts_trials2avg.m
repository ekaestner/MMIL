function outdata = ts_trials2avg(data,varargin)
% Purpose : average structures with 3D and 4D data matrices with respect to trials (epochs)
%
% Usage: outdata = ts_data_selection( indata, 'key1', val1, ... );
%
% Example: avg_data = ts_trials2avg(epoch_data);
%
% Inputs:
% 	 data: timesurfer epoch data structure 
%
% Outputs:
%    outdata: timesurfer average data structure
%
% Parameters: default : options :  description
%     events: [] :: event codes of conditions to process
%     conditions: [] :: condition number (not event code) used to index data
%     toilim: [] :: time limits ([begin end]) in seconds
%     toi: [] :: times of interest (sec)
%     keep_event_codes: 1 ::
%     stdev: 'data' :
%     verbose: 0 : 0,1 : whether to produce full output 
%     logfile: [] :: filename for log file 
%     write_avg_log: 0 : 0,1 : whether to create a logfile containing 
%                                 the number of trials in each average 
%     rootoutdir: pwd :: All output files will be saved in the rootoutdir
%                        directory (default is working directory)
%
% Created by:       Jason Sherfey 02-Feb-2009
% Last Modified by: Burke Rosen 18-Aug-2011 -changed mean to nanmean-
if nargin < 1, help(mfilename); end
data  = ts_checkdata_header(data);
parms = mmil_args2parms(varargin,...
						{'events',[],[],...
						 'conditions',[],[],...
						 'toilim',[],[],...
                         'toi',[],[],...
                         'keep_event_codes',1,[],...
                         'stdev','data',[],...
                         'verbose',1,{0,1},...
                         'logfile',     [],[],...
                         'datafield','averages',[],...
                         'write_avg_log', 0, {0,1},...
                         'rootoutdir',pwd,[],...
                                    },false);

[datatype,datafield,dataparam]   = ts_object_info(data,varargin{:});          

% extract events
if isempty(parms.events) && isempty(parms.conditions)
  conds = 1:length(data.(datafield));
elseif ~isempty(parms.conditions)
  conds = parms.conditions;
elseif ~isempty(parms.events)
  if iscell(parms.events)
    parms.events = unique([parms.events{:}]);
  end
  conds = find(ismember([data.(datafield).event_code],parms.events));
end
data.(datafield) = data.(datafield)(conds);

% abort if data is an average
if strcmp(datatype,'avg_data') || (strcmp(datatype,'timefreq_data') && ndims(data.(datafield)(1).(dataparam{1})) == 3)
  if parms.verbose, warning('%s: aborting trial averaging: data provided is already averaged over trials',mfilename); end
  outdata = data;
  return;
end

% add stdev field to indata
if ischar(parms.stdev) && any(strcmp(parms.stdev,dataparam))
  [data.(datafield).stdev] = deal([]);
end

%% Create average
outdata = rmfield(data,datafield);
for c = 1:length(data.(datafield))
  outdata.(parms.datafield)(c) = data.(datafield)(c);
  for k = 1:length(dataparam)
    parm = dataparam{k};
    ndim = ndims(data.(datafield)(c).(parm));
    if ndim < 3
      outdata.(parms.datafield)(c).(parm) = data.(datafield)(c).(parm);
      continue;
    end
    outdata.(parms.datafield)(c).(parm) = [];
    outdata.(parms.datafield)(c).(parm) = nanmean(data.(datafield)(c).(parm),ndim);    
    if ~isempty(parms.stdev) && strcmp(parms.stdev,parm)
      outdata.(parms.datafield)(c).stdev = nanstd(double(data.(datafield)(c).(parm)),0,ndim);
    end
    data.(datafield)(c).(parm) = [];
  end
end

%% create log
if parms.write_avg_log
    fid = fopen([parms.rootoutdir '/' parms.logfile],'wt');
    fprintf(fid, 'event_code\tnum_trials\n');
    for ievnt = 1:length([outdata.averages.event_code])
        event_code = outdata.averages(ievnt).event_code;
        num_trials = outdata.averages(ievnt).num_trials;
        fprintf(fid, '%i\t%i\n',event_code,num_trials);
    end
    fclose(fid);
end
    
    
