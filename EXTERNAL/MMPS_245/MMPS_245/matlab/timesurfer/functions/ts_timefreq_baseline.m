function outdata = ts_timefreq_baseline(datafile,varargin)

% toi & toilim units must match what is in the data structure

parms = mmil_args2parms(varargin,...
						{'events',[],[],...
						 'conditions',[],[],...
						 'toilim',[],[],...
             'toi',[],[],...
             'blcwindow',[],[],...
             'foi',[],[],...
             'baseline_start',-inf,[],...
             'baseline_end',0,[],...
             'keep_event_codes',1,[],...
             'out','yes',{'yes','avg_data','epoch_data','timefreq_data','timefreq'},...
             'reject_data',[],[],...
             'rejectfile',[],[],...
						},false);

if ~iscell(datafile), datafile = {datafile}; end

if ischar(parms.rejectfile) && exist(parms.rejectfile,'file')
  fprintf('%s: loading reject data: %s\n',mfilename,parms.rejectfile);
  load(parms.rejectfile);
  parms.reject_data = reject_data;
  clear reject_data;
end

% get toilims          
if isempty(parms.toilim)
  if ~isempty(parms.toi)
    parms.toilim = [min(parms.toi) max(parms.toi)];
  elseif ~isempty(parms.blcwindow)
    parms.toilim = [parms.blcwindow(1) parms.blcwindow(end)];
  else
    parms.toilim = [parms.baseline_start parms.baseline_end];
  end
end
tmin = parms.toilim(1);
tmax = parms.toilim(2);

base = [];
code = [];
nall = 0;
n    = 1;
fprintf('%s: calculating composite baseline from %g files\n',mfilename,length(datafile));
for f = 1:length(datafile)
  fprintf('%s: processing file %g of %g\n',mfilename,f,length(datafile));
  thisfile = datafile{f};
  if ~exist(thisfile,'file') || ~findstr(thisfile,'.mat'), continue; end
  S = load(thisfile);
  if ~isfield(S,'timefreq_data')
    clear S;
    continue;
  end
  data = S.timefreq_data;
  clear S;
  if ~isfield(data,'timefreq')
    clear data
    continue;
  end
  if ~isfield(data.timefreq,'power')
    clear data
    continue;
  end
  if ~isempty(parms.reject_data) && isstruct(parms.reject_data)
    data = ts_data_selection(data,'reject_data',parms.reject_data,'keepbadchans',1);
  end  
  nc   = length(data.timefreq);
  nall = nall + nc;
  for c = 1:nc
    t    = data.timefreq(c).time;
    % get indices to toi
    tix  = find(t >= tmin & t <= tmax);
    dat  = data.timefreq(c).power(:,tix,:,:);
    data.timefreq(c).power = [];
    % average over trials if necessary
    if ndims(dat) == 4, dat = nanmean(dat,4); end    
    if n == 1
      % 1st condition
      nch = size(dat,1);
      nt  = size(dat,2);
      nfq = size(dat,3);
      code = data.timefreq(c).event_code;
      ntrl = data.timefreq(c).num_trials;
      base = dat*ntrl;
      outdata = rmfield(data,'timefreq');
      outdata.timefreq.time = t(tix);
      outdata.timefreq.num_rejects = data.timefreq(c).num_rejects;
      outdata.timefreq.frequencies = data.timefreq(c).frequencies;
      fieldorder = setdiff(fieldnames(data.timefreq),{'cmplx','cross','plv','coh'});
    else
      % make sure subsequent dimensions are consistent
      if size(dat,1)~=nch || size(dat,2)~=nt || size(dat,3)~=nfq
        fprintf('%s: skipping event %g, inconsistent dimensions.\n',mfilename,data.timefreq(c).event_code);
        clear dat
        continue;
      end
      code = [code data.timefreq(c).event_code];
      ntrl = [ntrl data.timefreq(c).num_trials];
      base = dat*ntrl(end) + base;
    end
    clear dat
    n = n + 1;
  end
  clear data nc
end
if isempty(base)
  fprintf('%s: failed to find any data\n',mfilename);
  return;
end
outdata.timefreq.power      = base / sum(ntrl);
outdata.timefreq.event_code = code;
outdata.timefreq.num_trials = sum(ntrl);
outdata.timefreq = orderfields(outdata.timefreq,fieldorder);
if issubfield(outdata,'parms.filename')
  outdata = rmsubfield(outdata,'parms.filename');
end
fprintf('%s: finished averaging baseline data over %g of %g conditions found in %g files\n',mfilename,n-1,nall,length(datafile));
