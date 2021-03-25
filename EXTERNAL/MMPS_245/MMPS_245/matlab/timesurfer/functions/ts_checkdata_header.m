function [data,hdr] = ts_checkdata_header(data,varargin)
% todo: special handling for time-freq stats

parms = mmil_args2parms(varargin,...
						{'precision','skip',{'single','double','skip'},...       
             'events',[],[],...
						 'conditions',[],[],...
             'datatype',[],[],...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...                
						},false);
          
hdr = [];         
orig_events = parms.events;
if ts_object_info(data,'verbose',0)
  if length(data) > 1
    data = ts_combine_data(data);
  end
  if ~strcmp(parms.precision,'skip'), checkprecision; end
  return; 
end
if ischar(data)
  tmp = {data};
  clear data
  data.parms.filename = tmp;
end
if iscell(data) && ischar(data{1})
  tmp = data;
  clear data
  data.parms.filename = tmp;
end
if isstruct(data) && issubfield(data,'parms.filename')
  hdr = data;
  clear data;
  % call subfunction for time-freq stats
  if ~iscell(hdr.parms.filename)
    hdr.parms.filename = {hdr.parms.filename};
  end
  n = 0;
  rmfiles = [];
  for f = 1:length(hdr.parms.filename)
    if ~exist(hdr.parms.filename{f},'file'), rmfiles = [rmfiles f]; continue; end
    mmil_logstr(parms,'%s: scanning file: %s\n',mfilename,hdr.parms.filename{f});    
%     fprintf('%s: scanning file: %s\n',mfilename,hdr.parms.filename{f});    
    S = load(hdr.parms.filename{f});
    s = fieldnames(S);
    if ~isempty(parms.datatype) && ischar(parms.datatype)
      %% TODO: ENABLE CELL ARRAY OF DATATYPES (ex. avg_data & stat_data)
      if ismember(parms.datatype,s)
        s = parms.datatype;
      else
        rmfiles = [rmfiles f];
        continue;
      end
    else
      s = s{1};
    end
    try
      if ~isempty(S.(s)) && isstruct(S.(s))
        [datatype datafield] = ts_object_info(S.(s));
				if ~isempty(parms.conditions) && isempty(parms.events)
					parms.conditions = intersect(parms.conditions,1:length(S.(s).(datafield)));
					parms.events     = [S.(s).(datafield)(parms.conditions).event_code];
      	elseif isempty(parms.conditions) && isempty(orig_events)
          parms.events     = unique([parms.events S.(s).(datafield).event_code]);
				end				
        if ~isempty(parms.events) && ~any(ismember([S.(s).(datafield).event_code],parms.events))
          rmfiles = [rmfiles f];
          continue;
        end
      end
      keep = [];      
      for cc = 1:length(S.(s).(datafield))
        if ismember(S.(s).(datafield)(cc).event_code,parms.events)
          keep = [keep cc];
        end
      end
      S.(s).(datafield) = S.(s).(datafield)(keep); 
      data(n+1) = S.(s);
      n = n + 1;
      clear keep
      mmil_logstr(parms,'%s: found %s\n',mfilename,s);
%       fprintf('%s: found %s\n',mfilename,s);
    catch
      rmfiles = [rmfiles f];
      mmil_logstr(parms,'%s: failed to load %s\n',mfilename,hdr.parms.filename{f});
%       fprintf('%s: failed to load %s\n',mfilename,hdr.parms.filename{f});
    end
    clear S s
  end
  hdr.parms.filename(rmfiles) = [];
  if n > 1
    if issubfield(data,'parms.filename')
      filenames = {};
      for c = 1:n
        filenames = {filenames{:} data(c).parms.filename{:}};
      end
      data = ts_combine_data(data); 
      data.parms.filename = filenames;
    else
      data = ts_combine_data(data); 
    end
    % make sure the event order in data matches parms.events
    if ~isempty(parms.events)
      evdata = [data.(datafield).event_code];
      evparm = parms.events;
      [a ia ib] = intersect(evparm,evdata);
      data.(datafield) = data.(datafield)(ib(ia));
    end    
  end
end
if ~strcmp(parms.precision,'skip'), checkprecision; end
if ~exist('data','var')
  mmil_logstr(parms,'%s: no data found.\n',mfilename);
%   fprintf('%s: no data found.\n',mfilename);
  data = [];
end
  function checkprecision
    [datatype datafield dataparam] = ts_object_info(data);
    for p = 1:length(dataparam)
      for cc = 1:length(data.(datafield))
        if ~strcmp(class(data.(datafield)(cc).(dataparam{p})),parms.precision)
          mmil_logstr(parms,'%s: converting data to %s precision.\n',mfilename,parms.precision);
          if strcmp(parms.precision,'single')
            if ~isa(data.(datafield)(cc).(dataparam{p}),'single')
              data.(datafield)(cc).(dataparam{p}) = single(data.(datafield)(cc).(dataparam{p}));
            end
          elseif strcmp(parms.precision,'double')
            if isa(data.(datafield)(cc).(dataparam{p}),'single')
              data.(datafield)(cc).(dataparam{p}) = double(data.(datafield)(cc).(dataparam{p}));
            end
          end
        end
      end
    end
  end
end
