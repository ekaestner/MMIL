function plv_data = ts_PLV(timefreq_data,varargin)
% data (complex): chan x time x freq x trials
% timefreq_data should be an average over trials containing a parms
% structure with its filename.  The same directory must have a
% timefreq_data structure for each channel containing all trials.  The name
% of those structures must match the name of the average file with _chan###
% appended to it. (ex. timefreq_cond3.mat ==> timefreq_cond3_chan001.mat)
% foi: vector frequencies over which to calc PLV
% events: vector of event codes for which to calc PLV (default = first condition)
% band: [0|1] whether to compute an average PLV over a frequency band
% channelcmb: either
%   (1) cell array of vectors of channels indices for which to calc PLV
%       ex. {[1 2] [1 3]}
%   (2) 2D cell array of channel labels of channels for which to calc PLV
%       ex. [{'EEG1' 'EEG2'}; {'EEG1' 'EEG3'}]
%   (3) 'all' - to compute for all channel pairs
% output: {'plv','coh','coherency','mscoh','lag'}
tic
parms = mmil_args2parms(varargin,...
						{'foi',[],[],...
             'foilim',[],[],...
             'toilim',[],[],...
						 'band',0,[],...
             'events',[],[],...
						 'channelcmb','all',[],...
             'output','plv',[],...
             'loadflag',0,{0,1},...
             'load_flag',[],[],...
             'combine_flag',0,{0,1},...
             'neweventcodes',[],[],...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...    
             'dsfact',1,[],...
						},false);
if ~iscell(parms.output), parms.output={parms.output}; end
if ~iscell(parms.events), parms.events={parms.events}; end
evcodes = unique([parms.events{:}]);
if ~isempty(parms.load_flag) && isnumeric(parms.load_flag), parms.loadflag = parms.load_flag; end
if parms.loadflag
%   if issubfield(timefreq_data,'parms.filename') && iscell(timefreq_data.parms.filename)
%     timefreq_data.parms.filename = timefreq_data.parms.filename(1);
%   end
  [timefreq_data,hdr]  = ts_checkdata_header(timefreq_data,'events',evcodes);
  fnames = hdr.parms.filename;
  if ~iscell(fnames), fnames = {fnames}; end
end
timefreq_data = backwardcompatible(timefreq_data,parms);
[datatype datafield] = ts_object_info(timefreq_data);
if strcmp(datafield,'epochs')
  tfrflag = 0;
  % mscohere parameters
  % nfft = ...
  fs       = data.sfreq;
  window   = hanning(nfft); % windowing function
  noverlap = nfft/2;        % # of samples by which sections overlap (default = 50% overlap)  
  if mod(nfft,2)            % odd nfft
    nfreq = (nfft+1)/2;
  else
    nfreq = nfft/2+1;
  end
elseif strcmp(datafield,'timefreq')
  tfrflag = 1;
end

% events
if isempty(evcodes)
  evcodes      = [timefreq_data.(datafield).event_code];
  parms.events = {evcodes};
end

% frequencies
if tfrflag, f = timefreq_data.(datafield)(1).frequencies; else f = 1; end
if ~isempty(parms.foi)
  [jk1 fix jk2] = intersect(f,parms.foi);
  fix           = sort(fix);
elseif ~isempty(parms.foilim)
  fix           = find(parms.foilim(1)<=f & f<=parms.foilim(2));
else
  fix           = 1:length(f);
end
f     = f(fix);
nfreq = length(f);

% time
if parms.dsfact ~= 1 && ~isempty(parms.dsfact)
  if ~isempty(parms.toilim)
    timefreq_data = ts_data_selection(timefreq_data,'toilim',parms.toilim);
  end    
  for k = 1:length(timefreq_data.(datafield))
    timefreq_data.(datafield)(k).time = downsample(timefreq_data.(datafield)(k).time,parms.dsfact);
  end
  timefreq_data.sfreq = timefreq_data.sfreq / parms.dsfact;
end
ntime = length(timefreq_data.(datafield)(1).time);

sens  = timefreq_data.sensor_info;

if parms.loadflag
%   [path fname] = fileparts(timefreq_data.parms.filename);
%   fname = [path '/' fname];
  timefreq_data.(datafield) = rmfield(timefreq_data.(datafield),...
                           intersect(fieldnames(timefreq_data.(datafield)),{'power','data','cross','cmplx'}));
  plv_data = timefreq_data;
  plv_data.(datafield) = plv_data.(datafield)(1);
  clear timefreq_data
else
  % events
  [evcodes ix jx]           = intersect(evcodes,[timefreq_data.(datafield).event_code]);
  timefreq_data.(datafield) = timefreq_data.(datafield)(jx);
  
  if isfield(timefreq_data.(datafield),'power')
    timefreq_data.(datafield) = rmfield(timefreq_data.(datafield),'power');
  end
  plv_data = timefreq_data;
  plv_data.timefreq = rmfield(plv_data.timefreq,...
                      intersect(fieldnames(plv_data.timefreq),{'power','data','cross','cmplx'}));
end

cmb = [];
% channel pairs for which to calc PLV
if isnumeric(parms.channelcmb) && size(parms.channelcmb,2)==2
  cmb = parms.channelcmb;
elseif iscell(parms.channelcmb) && isnumeric(parms.channelcmb{1,1})
  cmb = reshape([parms.channelcmb{:}],2,size(parms.channelcmb,1))';
elseif iscell(parms.channelcmb) && ischar(parms.channelcmb{1,1})
  lbl = {parms.channelcmb{:}};
  [jnk sel] = match_str(lbl,{sens.label}); % note: {sens(sel).label} == lbl
  cmb = reshape(sel,size(parms.channelcmb,1),2);
elseif ischar(parms.channelcmb) && ismember(parms.channelcmb,{sens.label})
  [chx jnk] = match_str({sens.label},parms.channelcmb);
  [chs jnk] = match_str({sens.typestring},sens(chx).typestring);
  cmb = [chx*ones(length(chs),1) chs];
else
  % default to all
  cmb = selectall(sens);
end
if ~isnumeric(cmb) || size(cmb,2)~=2
  error('invalid specification of channel pairs.');
end
% cmb is now a matrix of unique channel index pairs
numcmb = size(cmb,1);

if ismember('plv',      parms.output), plvflg = 1; else plvflg = 0; end   % phase locking value
if ismember('coherency',parms.output), coyflg = 1; else coyflg = 0; end   % coherency (complex)
if ismember('coh',      parms.output), cohflg = 1; else cohflg = 0; end   % coherence |complex|
if ismember('mscoh',    parms.output), mscflg = 1; else mscflg = 0; end   % magnitude squared coherence (real) - analogous to Pearson correlation r^2
if ismember('lag',      parms.output), lagflg = 1; else lagflg = 0; end   % phase lag

if ~parms.combine_flag
  parms.events = num2cell(evcodes);
else
  if isempty(parms.neweventcodes) || length(parms.neweventcodes)~=length(parms.events)
    maxval = max([evcodes [plv_data.(datafield).event_code]]);
    parms.neweventcodes = [maxval+1:maxval+length(parms.events)];
  end
end

% calculate plv
for e = 1:length(parms.events)
  events = parms.events{e};
  if coyflg, coy = zeros(numcmb,ntime,nfreq); end
  if cohflg, coh = zeros(numcmb,ntime,nfreq); end
  if mscflg, msc = zeros(numcmb,ntime,nfreq); end  
  if lagflg, lag = zeros(numcmb,ntime,nfreq); end
  if plvflg
    if parms.band
      plv = zeros(numcmb,ntime);
    else
      plv = zeros(numcmb,ntime,nfreq); 
    end
  end
  plv_data.timefreq(e).labelcmb    = [];
  plv_data.timefreq(e).frequencies = f;
  for k=1:numcmb
    i = cmb(k,1);
    j = cmb(k,2);
    mmil_logstr(parms,'%s: calculating %s matrix %i of %i (%s,%s) for event(s): [%s] (%gs elapsed)\n',mfilename,[parms.output{:}],k,numcmb,sens(i).label,sens(j).label,num2str(events),toc);    
    if parms.verbose
      fprintf('%s: calculating %s matrix %i of %i (%s,%s) for event(s): [%s] (%gs elapsed)\n',mfilename,[parms.output{:}],k,numcmb,sens(i).label,sens(j).label,num2str(events),toc);
    end
    % use mscohere if given epoch_data
    % ...
    % continue;
    if parms.loadflag
      fname = fnames(ismember(evcodes,events));
      data  = getdata(fname,events,i,j,parms.toilim,parms.dsfact);
    else
      data = ts_data_selection(timefreq_data,'channels',[i j],'removebadchans',1,'events',events);
      if length(events) > 1
        data = {data.timefreq.cmplx};
        data = cat(4,data{:});
      else
        data = data.timefreq.cmplx;
      end
      if ~isa(data,'double')
        data = double(data);
      end
      if i == j
        data(2,:,:,:) = data(1,:,:,:);
      end
      if parms.dsfact ~=1 && ~isempty(parms.dsfact)
        data = data(:,1:parms.dsfact:end,:,:);
      end      
    end
    ntrl = size(data,4);
    if plvflg
      phase1 = angle(squeeze(data(1,:,fix,:)));
      phase2 = angle(squeeze(data(2,:,fix,:)));
      if parms.band
        plv(k,:) = getPLV3D(phase1,phase2,1)';
      else
        plv(k,:,:) = getPLV3D(phase1,phase2,0);
      end
      clear phase1 phase2
    end
    if cohflg || coyflg || mscflg || lagflg
      x   = squeeze(data(1,:,fix,:)); % time x freq x trials
      y   = squeeze(data(2,:,fix,:));
      sxy = sum(x.*conj(y),3) / size(x,3); % time x freq - cross-spectrum analogous to covariance
      pxx = sum(x.*conj(x),3) / size(x,3);
      pyy = sum(y.*conj(y),3) / size(y,3);
      if coyflg, coy(k,:,:) = sxy ./ sqrt(pxx .* pyy);      end   % coherency (complex)
      if cohflg, coh(k,:,:) = abs(sxy ./ sqrt(pxx .* pyy)); end
      if mscflg, msc(k,:,:) = abs(sxy).^2 ./ (pxx.*pyy);    end % magnitude squared coherence (real) - analogous to Pearson correlation r^2
      if lagflg, lag(k,:,:) = atan(imag(sxy) ./ real(sxy)); end % phase difference, [-pi,pi]
    end
    plv_data.timefreq(e).labelcmb = [plv_data.timefreq(e).labelcmb; {sens(i).label sens(j).label}];
    clear data;
  end
  if parms.combine_flag
    plv_data.timefreq(e).event_code = parms.neweventcodes(e);
  end  
  if plvflg, plv_data.timefreq(e).plv      = single(plv); end
  if coyflg, plv_data.timefreq(e).cohcmplx = single(coy); end
  if cohflg, plv_data.timefreq(e).coh      = single(coh); end
  if mscflg, plv_data.timefreq(e).mscoh    = single(msc); end   
  if lagflg, plv_data.timefreq(e).lag      = single(lag); end
	plv_data.timefreq(e).num_trials          = ntrl;
	clear plv coy coh msc lag
end
toc

% subfunctions
function plv = getPLV3D(phase1,phase2,fband)
% Calculate plv for phases 1 and 2
% phase1, phase2: time x freq x trials
ntrials = size(phase1,3);
phase = phase1-phase2;				% calc phase difference
phase = exp(j*phase);					% make complex
phase = sum(phase,3)/ntrials; % make complex circular mean
if fband
	phase = mean(phase,2);				% average over freq band
end
plv 	= abs(phase);						% get modulus --> plv

function timefreq_data = backwardcompatible(timefreq_data,parms)
if isfield(timefreq_data,'opt')
  timefreq_data.parms = timefreq_data.opt;
  timefreq_data = rmfield(timefreq_data,'opt');
end
if ~issubfield(timefreq_data,'parms.filename') && parms.loadflag
  error('timefreq_data filename not found.');
elseif parms.loadflag
  if iscell(timefreq_data.parms.filename)
    timefreq_data.parms.filename = timefreq_data.parms.filename{1}; 
  end
%   load(timefreq_data.parms.filename);
%   if issubfield(timefreq_data,'parms.filename') && iscell(timefreq_data.parms.filename)
%     timefreq_data.parms.filename = timefreq_data.parms.filename{1};
%   end
end

function data = getdata(fname,e,i,j,toilim,dsfact)
for f = 1:length(fname)
  file   = fname{f}; [filepath,filename] = fileparts(file);
  prefix = [filepath '/' filename];
  
  S = dir(filepath);
  files = {S.name};
  files = files(~[S.isdir]);
  ind   = cellfun(@(x)(strfind(x,filename)),files,'UniformOutput',false);
  files = files(~cellfun(@isempty,ind));
  
  % first channel
  ind   = cellfun(@(x)(strfind(x,sprintf('chan%03i',i))),files,'UniformOutput',false);
  ind   = find(~cellfun(@isempty,ind));
  if numel(ind)==1
    infile = fullfile(filepath,files{ind});
  else
    infile = sprintf('%s_chan%03i.mat',prefix,i);
  end
  load(infile); % timefreq_data
%   if ~isempty(toilim)
%     timefreq_data = ts_data_selection(timefreq_data,'toilim',toilim);
%   end
  tmp = {timefreq_data.timefreq(ismember([timefreq_data.timefreq.event_code],e)).cmplx};
  tmp = cat(ndims(tmp{1}),tmp{:});
  if ~isempty(toilim)
    tmp = tmp(:,nearest(timefreq_data.timefreq(1).time,toilim(1)):nearest(timefreq_data.timefreq(1).time,toilim(2)),:,:);
  end
  clear timefreq_data
  if dsfact ~=1 && ~isempty(dsfact)
    tmp = tmp(:,1:dsfact:end,:,:);
  end
  if ~isa(tmp,'double'), tmp = double(tmp); end
  if f == 1
    dat1 = tmp;
  else
    dat1 = cat(ndims(tmp),dat1,tmp);
  end
  clear tmp
  % second channel
  ind   = cellfun(@(x)(strfind(x,sprintf('chan%03i',j))),files,'UniformOutput',false);
  ind   = find(~cellfun(@isempty,ind));
  if numel(ind)==1
    infile = fullfile(filepath,files{ind});
  else
    infile = sprintf('%s_chan%03i.mat',prefix,j);
  end
  load(infile); % timefreq_data
%   if ~isempty(toilim)
%     timefreq_data = ts_data_selection(timefreq_data,'toilim',toilim);
%   end
  tmp = {timefreq_data.timefreq(ismember([timefreq_data.timefreq.event_code],e)).cmplx};
  tmp = cat(ndims(tmp{1}),tmp{:});
  if ~isempty(toilim)
    tmp = tmp(:,nearest(timefreq_data.timefreq(1).time,toilim(1)):nearest(timefreq_data.timefreq(1).time,toilim(2)),:,:);
  end
  clear timefreq_data
  if dsfact ~=1 && ~isempty(dsfact)
    tmp = tmp(:,1:dsfact:end,:,:);
  end  
  if ~isa(tmp,'double'), tmp = double(tmp); end
  if f == 1
    dat2 = tmp;
  else
    dat2 = cat(ndims(tmp),dat2,tmp);
  end
  clear tmp
end
data = cat(1,dat1,dat2);

% function data = getdata(fname,e,i,j)
% load(sprintf('%s_chan%03i.mat',fname,i));
% data = double(timefreq_data.timefreq([timefreq_data.timefreq.event_code]==e).cmplx);
% load(sprintf('%s_chan%03i.mat',fname,j));
% data = cat(1,data,double(timefreq_data.timefreq([timefreq_data.timefreq.event_code]==e).cmplx));

function cmb = selectall(sens)
nchan = length(sens);
c1  = [repmat(1:nchan,[nchan 1])']'; c1 = c1(:);
c2  = repmat([1:nchan]',[nchan 1]);
cmb = [c1 c2];
cup = cmb(c2> c1,:);  % upper triangular off-diagonal matrix
cdn = cmb(c2< c1,:);  % lower triangular off-diagonal matrix
cdg = cmb(c2==c1,:);  % diagonal
[jnk cupidx] = setdiff(cup,fliplr(cdn),'rows');
[jnk cdnidx] = setdiff(fliplr(cdn),cup,'rows');
[jnk ix tss] = intersect(cup,fliplr(cdn),'rows');
cupidx = [cupidx ix];
cmb = sortrows([cup(cupidx,:); cdg; cdn(cdnidx,:)]);