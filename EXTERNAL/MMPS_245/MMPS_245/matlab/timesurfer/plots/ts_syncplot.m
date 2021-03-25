function syncplot(varargin)
% Purpose: wrap around ts_ezplot looping over rechans for synchrony plots
% Parameters:
%   datafile: string or cell array of strings listing files w/ plv_data (one
%   per condition)

if isstruct(varargin{1})
  plv_data = varargin{1};
  varargin = varargin(2:end);
  load_flag = 0;
else
  load_flag = 1;
end

parms = mmil_args2parms(varargin,...
				{  'verbose'  ,0    ,{0,1},...
           'logfid'   ,[1]  ,[],...
           'logfile'  ,[]   ,[],...
           'datafile' ,[]   ,[],...
           'zparam'   ,'plv',{'plv','coh','coherency','mscoh','lag'},...
           'prefix'   ,'proc',[],...
           'title'    ,'sync',[],...
					 'events' 	,[]    ,[],...
           'outpath'  ,'images/sync',[],...
           'refchan'  ,'all',[],...
				},false);

% remove prefix & title from varargin if present
ix = find(strcmp('prefix' ,varargin(1:2:end))); ix=2*ix-1; if ~isempty(ix), varargin([ix ix+1]) = []; end
ix = find(strcmp('title'  ,varargin(1:2:end))); ix=2*ix-1; if ~isempty(ix), varargin([ix ix+1]) = []; end
ix = find(strcmp('zparam' ,varargin(1:2:end))); ix=2*ix-1; if ~isempty(ix), varargin([ix ix+1]) = []; end
ix = find(strcmp('outpath',varargin(1:2:end))); ix=2*ix-1; if ~isempty(ix), varargin([ix ix+1]) = []; end
ix = find(strcmp('events' ,varargin(1:2:end))); ix=2*ix-1; if ~isempty(ix), varargin([ix ix+1]) = []; end

if isempty(parms.datafile) && ~exist('plv_data','var')
  [filename,pathname] = uigetfile({'*.mat'},'Load plv_data.','MultiSelect','off');
  if isequal(filename,0) || isequal(pathname,0), return; end
  parms.datafile = [pathname filename];
elseif exist('plv_data','var')
	plv_data = ts_checkdata_header(plv_data,'events',parms.events);
end
if ~iscell(parms.datafile), parms.datafile = {parms.datafile}; end

for c  = 1:length(parms.datafile)
  if load_flag, load(parms.datafile{c}); end
  if ~isempty(parms.events)
    plv_data = ts_data_selection(plv_data,'events',parms.events);
  end
  sen = {plv_data.sensor_info.label};
  cmb = plv_data.timefreq(1).labelcmb;
  
  if isempty(parms.refchan) || any(strcmpi(parms.refchan,'all'))
    allrefchans = {plv_data.sensor_info.label};
  elseif isnumeric(parms.refchan)
    allrefchans = {plv_data.sensor_info(parms.refchan).label};
  elseif ischar(parms.refchan)
    allrefchans = {parms.refchan};
  elseif iscellstr(parms.refchan)
    allrefchans = parms.refchan;
  else
    error('Reference channels not specified.');
  end
  if ~any(ismember(allrefchans,{plv_data.sensor_info.label}))
    error('The reference channel(s) were not found in the data.');
  end  
  for refchan = 1:length(sen)
    ref  = sen{refchan};  
    data = rmfield(plv_data,'timefreq');
    if ~ismember(ref,allrefchans)
      continue;
    end
    
    sel1 = strmatch(ref,cmb(:,1),'exact');
    lab1 = cmb(sel1,2);
    sel2 = strmatch(ref,cmb(:,2),'exact');
    lab2 = cmb(sel2,1);
    sel  = [sel1;sel2];
    lab  = [lab1;lab2];

    % remove duplicate entry w/o sorting (note: "unique" always sorts)
    tmp  = strmatch(ref,lab,'exact');
    if length(tmp) > 1
      tmp = tmp(2:end);
      sel(tmp) = [];
      lab(tmp) = [];
    end
    clear tmp

    [sel1 sel2]      = match_str(lab,sen);
    data.sensor_info = data.sensor_info(sel2);
    data.num_sensors = length(data.sensor_info);
    for k = 1:length(plv_data.timefreq)
      tmp              = rmfield(plv_data.timefreq(k),intersect(fieldnames(plv_data.timefreq),{'plv','lag','labelcmb','coh','mscoh','coherency'}));    
      tmp.power        = [];
      data.timefreq(k) = tmp; clear tmp
      data.timefreq(k).power = plv_data.timefreq(k).(parms.zparam)(sel,:,:);
    end
    title   = sprintf('%s %s: refchan = %s',parms.title,upper(parms.zparam),ref);
    prefix  = sprintf('%s_%s_Reference-%s',parms.prefix,upper(parms.zparam),ref);
    for k = 1:length(data.timefreq)
      this    = data.timefreq(k).event_code;
      outpath = sprintf('%s/event%g',parms.outpath,this);
      ts_ezplot(data,'events',this,'title',title,'prefix',prefix,'outpath',outpath,varargin{:});
    end
    clear data
  end

end
