function zdata = ts_zscore(data,varargin)
% Purpose: transform data by calculating z-score or performing baseline correction.
%
% Usage: zdata = ts_process_data(data,'key1',val1, ... )
%
% Example: zdata = ts_zscore(data,'blcwindow',[-.5 -.2]);
%
% Inputs:
%     data: timesurfer data struct of epoch, average, timefreq, or plv data
%
% Outputs:
%     zdata:  timesurfer data struct of the same type after transformation
%
% Parameters: default : options :  description

%     toilim: [] :: limits on the times of interest, [begin end] in seconds
%     toi: [] :: array of times to keep
%     blcwindow: [] :: [begin end] in seconds, the default is the complete trial
%     events: [] :: vector of event codes of conditions to process
%     baseline_data: [] :: a second TimeSurfer structure containing alternative data 
%                          to use for the baseline estimation
%     baselinefile: [] :: file containing baseline_data
%     baselinetype: 'zscore' : 'absolute','relchange','relative','zscore': 
%                             the kind of transformation to perform (see note below)
%     skipzero: 0: 0,1 : whether to skill conditions containing channels with zeros
%     zparam: 'all' :: fieldname containing data to transform
%     verbose: 1 : 1,0 : whether to produce fulll output
%     logfile: [] :: filename for log file (stdout will be appended)
%     logfid: 1 :: file FID for log file
%
% Note on baselinetype.
% absolute:   data - mean
% relchange:  (data - mean) / mean
% relative:   data / mean
% zscore:     (data - mean) / (standard deviation)
%   where mean and standard deviation are defined with respect to blcwindow

warning off
parms = mmil_args2parms(varargin,...
						{'tdim',2,[],...
						 'toi',[],[],...
						 'toilim',[],[],...
						 'blwindow',[],[],...
                         'blcwindow',[],[],...
                         'zbase',[],[],...
						 'events',[],[],...
						 'conditions',[],[],...
						 'baseline_data',[],[],...
						 'baselinefile',[],[],...
                         'baselinetype','zscore',{'absolute','relchange','relative','zscore'},...
                         'skipzero',0,{0,1},...
                         'zparam','all',[],...
                         'verbose',0,{0,1},...
                         'logfile',      [],[],...
                         'logfid',       [1],[], ...   
						},false);

if isempty(parms.blwindow) && ~isempty(parms.blcwindow)
  parms.blwindow = parms.blcwindow;
elseif isempty(parms.blwindow) && ~isempty(parms.zbase)
  parms.blwindow = parms.zbase;
end
  
data = ts_data_selection(data,varargin{:});
[datatype,datafield,dataparam] = ts_object_info(data,varargin{:});
tdim = parms.tdim;

if isfield(data.(datafield),'cmplx') && ~isfield(data.(datafield),'power')
  for i = 1:length(data.(datafield))
    data.(datafield)(i).power = 2*abs(data.(datafield)(i).cmplx).^2;
  end
%     data.(datafield) = rmfield(data.(datafield),'cmplx');
    dataparam{1} = 'power';
end
if ischar(parms.zparam) && ismember(parms.zparam,dataparam)
  clear dataparam
  dataparam{1} = parms.zparam;
elseif iscell(parms.zparam) && ischar(parms.zparam) && ismember(parms.zparam{1},dataparam)
  dataparam = intersect(parms.zparam,dataparam);
end
bflag = 0;
if ~isempty(parms.baselinefile) && exist(parms.baselinefile,'file')
  load(parms.baselinefile);  % baseline_data
elseif isstruct(parms.baseline_data)
  baseline_data = parms.baseline_data;
end
if exist('baseline_data','var')
  [jnk1,bfield,bparam] = ts_object_info(baseline_data,varargin{:});
%  if length(intersect(bparam,dataparam))==length(bparam)
try mmil_logstr(parms,'using baseline data combined over events [%s]',num2str(baseline_data.(bfield).event_code)); end  
% try fprintf('%s: using baseline data combined over events %s\n',mfilename,num2str(baseline_data.(bfield).event_code)); end
  baseline = baseline_data;
  [blchans jnk] = match_str({baseline.sensor_info.label},{data.sensor_info.label});
  
  bflag = 1;
%  end
end
if ~bflag
  if isempty(parms.blwindow)
    bidx = find(data.(datafield)(1).time < 0);
  elseif isnumeric(parms.blwindow)
    bidx = find(data.(datafield)(1).time >= parms.blwindow(1) & data.(datafield)(1).time <= parms.blwindow(end));
%     bidx = nearest(data.(datafield)(1).time,parms.blwindow(1)):...
%            nearest(data.(datafield)(1).time,parms.blwindow(end));
  elseif ischar(parms.blwindow) && strcmp(parms.blwindow,'all')
    bidx = 1:length(data.(datafield)(1).time);
  end
  if isempty(bidx) || length(bidx)<=1
    bidx = 1:length(data.(datafield)(1).time);
    if isempty(bidx)
      mmil_logstr(parms,'No baseline found.  Aborting z-score calculation.');
%       fprintf('%s: No baseline found.  Aborting z-score calculation.\n',mfilename);
      zdata = data;
      return;
    end
  end
end  

zdata = rmfield(data,datafield);  
if parms.verbose
  mmil_logstr(parms,'performing %s baseline correction on %g conditions.',...
      parms.baselinetype,length(data.(datafield)));  
%   fprintf('%s: performing %s baseline correction on %g conditions.\n',...
%     mfilename,parms.baselinetype,length(data.(datafield)));
end
for c = 1:length(data.(datafield))
  zdata.(datafield)(c) = data.(datafield)(c);
  for p = 1:length(dataparam)
    if strcmp(dataparam{p},'cmplx'),continue; end
    zdata.(datafield)(c).(dataparam{p}) = [];
    if bflag
      bdat = baseline.(bfield).(dataparam{p})(blchans,:,:,:);
    else
      bdat = data.(datafield)(c).(dataparam{p})(:,bidx,:,:);
    end
    repvec = ones(1,ndims(data.(datafield)(c).(dataparam{p})));
    repvec(tdim) = size(data.(datafield)(c).(dataparam{p}),tdim);
    mu    = mean(bdat,tdim);
    mu    = repmat(mu,repvec);    % average over time
    if strcmp(parms.baselinetype,'zscore')
      if ~strcmp(class(bdat),'double')
        sigma = repmat(std(double(bdat),0,tdim),repvec);
      else
        sigma = repmat(std(bdat,0,tdim),repvec);
      end
    end
    if isempty(find(mu==0)) || ~parms.skipzero
      dd           = ndims(data.(datafield)(c).(dataparam{p}));
      if ndims(mu) == dd - 1
        tmpvec      = ones(1,dd);
        tmpvec(end) = size(data.(datafield)(c).(dataparam{p}),dd);
        mu    = repmat(mu,tmpvec);
        if exist('sigma','var'), sigma = repmat(sigma,tmpvec); end
      end
      if strcmp(parms.baselinetype,'zscore')
        zdat  = (data.(datafield)(c).(dataparam{p}) - mu) ./ sigma;
      elseif strcmp(parms.baselinetype,'absolute')
        zdat  = (data.(datafield)(c).(dataparam{p}) - mu);
      elseif strcmp(parms.baselinetype,'relative')
        zdat  = (data.(datafield)(c).(dataparam{p})) ./ mu;
      elseif strcmp(parms.baselinetype,'relchange')
        zdat  = (data.(datafield)(c).(dataparam{p}) - mu) ./ mu;
      else
        mmil_logstr(parms,'baseline correction type not recognized.');
%         fprintf('baseline correction type not recognized.\n');
        clear mu sigma zdat bdat;
        zdata.(datafield)(c).(dataparam{p}) = data.(datafield)(c).(dataparam{p});
        continue;
      end
      zdat(abs(zdat) == inf) = nan;
      zdata.(datafield)(c).(dataparam{p}) = zdat;
    else
      mmil_logstr(parms,'zero baseline in at least one channel. skipping correction for condition %g',data.(datafield)(c).event_code);
%       fprintf('zero baseline in at least one channel. skipping correction for condition %g\n',data.(datafield)(c).event_code);
      zdata.(datafield)(c).(dataparam{p}) = data.(datafield)(c).(dataparam{p});
    end
    clear mu sigma zdat bdat;
  end
end
warning on