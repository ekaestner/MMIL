function zdata = ts_blc(data,varargin)

parms = mmil_args2parms(varargin,...
						{'tdim',2,[],...
						 'toi',[],[],...
						 'toilim',[],[],...
						 'blwindow',[],[],...
             'blcwindow',[],[],...
             'baseline_start',-inf,[],...
             'baseline_end',0,[],...
             'zbase',[],[],...
						 'events',[],[],...
						 'conditions',[],[],...
						 'baseline_data',[],[],...
						 'baselinefile',[],[],...
             'baseline','no',{'absolute','relchange','relative','zscore','no'},...             
             'baselinetype','absolute',{'absolute','relchange','relative','zscore','no'},...
             'skipzero',0,{0,1},...
             'verbose',1,{0,1},...
						},false);

% backwards compatibility
if ~strcmp(parms.baseline,'no')
  parms.baselinetype = parms.baseline;
end
if isempty(parms.blwindow) && ~isempty(parms.blcwindow)
  parms.blwindow = parms.blcwindow;
elseif isempty(parms.blwindow) && ~isempty(parms.zbase)
  parms.blwindow = parms.zbase;
else
  parms.blwindow = [parms.baseline_start parms.baseline_end];
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
elseif strcmp(datatype,'timefreq_data')
  clear dataparam;
  dataparam{1} = 'power';
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
  baseline = baseline_data;
  [blchans jnk] = match_str({baseline.sensor_info.label},{data.sensor_info.label});
  
  bflag = 1;
%  end
end
if ~bflag
  if isempty(parms.blwindow)
    bidx = find(data.(datafield)(1).time < 0);
  elseif isnumeric(parms.blwindow)
    bidx = find(data.(datafield)(1).time > parms.blwindow(1) & data.(datafield)(1).time < parms.blwindow(end));
%     bidx = nearest(data.(datafield)(1).time,parms.blwindow(1)):...
%            nearest(data.(datafield)(1).time,parms.blwindow(end));
  elseif ischar(parms.blwindow) && strcmp(parms.blwindow,'all')
    bidx = 1:length(data.(datafield)(1).time);
  end
  if isempty(bidx)
    fprintf('%s: No baseline found.  Aborting z-score calculation.\n',mfilename);
    zdata = data;
    return;
  end
end  

zdata = rmfield(data,datafield);  
if parms.verbose
  fprintf('%s: performing %s baseline correction on %g conditions.\n',...
    mfilename,parms.baselinetype,length(data.(datafield)));
end
for c = 1:length(data.(datafield))
  zdata.(datafield)(c) = data.(datafield)(c);
  for p = 1:length(dataparam)
    zdata.(datafield)(c).(dataparam{p}) = [];
    if bflag
      bdat = baseline.(bfield).(dataparam{p})(blchans,:,:,:);
    else
      bdat = data.(datafield)(c).(dataparam{p})(:,bidx,:,:);
    end
    if tdim ~= 1
      zdat         = [];
      repvec       = ones(1,ndims(bdat));
      repvec(tdim) = size(data.(datafield)(c).(dataparam{p}),tdim);
      for k = 1:size(bdat,1)
        if data.sensor_info(k).badchan, continue; end
        mu = repmat(nanmean(bdat(k,:,:),tdim),repvec);
        if ~strcmp(class(bdat),'double')
          sigma = repmat(nanstd(double(bdat(k,:,:)),0,tdim),repvec);
        else
          sigma = repmat(nanstd(bdat(k,:,:),0,tdim),repvec);
        end
        if strcmp(parms.baselinetype,'zscore')
          zdat  = cat(1,zdat,(data.(datafield)(c).(dataparam{p})(k,:,:,:) - mu) ./ sigma);
        elseif strcmp(parms.baselinetype,'absolute')
          zdat  = cat(1,zdat,(data.(datafield)(c).(dataparam{p})(k,:,:,:) - mu));
        elseif strcmp(parms.baselinetype,'relative')
          zdat  = cat(1,zdat,(data.(datafield)(c).(dataparam{p})(k,:,:,:)) ./ mu);
        elseif strcmp(parms.baselinetype,'relchange')
          zdat  = cat(1,zdat,(data.(datafield)(c).(dataparam{p})(k,:,:,:) - mu) ./ mu);
        else        
          fprintf('baseline correction type not recognized.\n');
          clear mu sigma zdat bdat;
          zdata.(datafield)(c).(dataparam{p}) = data.(datafield)(c).(dataparam{p});
          continue;
        end    
        clear mu sigma
      end
      zdat(abs(zdat) == inf) = nan;
      zdata.(datafield)(c).(dataparam{p}) = zdat;
      clear zdat
    else
      fprintf('Warning: %s does not support averaging along the first dimension\nCalling ts_zscore...\n',mfilename);
      zdata = ts_zscore(data,varargin{:});
      return;
    end
    clear bdat
  end
end
