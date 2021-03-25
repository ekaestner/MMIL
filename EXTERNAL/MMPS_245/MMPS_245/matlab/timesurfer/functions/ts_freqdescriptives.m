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

parms = mmil_args2parms(varargin,...
						{'foi',[],[],...
						 'band',0,[],...
             'events',[],[],...
						 'channelcmb','all',[],...
             'verbose',1,{0,1},...
             'output','plv',{'plv','coh','coherency','mscoh'},...
             'loadflag',0,{0,1},...
						},false);

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
if isempty(parms.events), parms.events = timefreq_data.(datafield)(1).event_code; end
[events ix jx]         = intersect(parms.events,[timefreq_data.(datafield).event_code]);
timefreq_data.(datafield) = timefreq_data.(datafield)(jx);

sens  = timefreq_data.sensor_info;
ntime = length(timefreq_data.(datafield)(1).time);
if tfrflag
  nfreq = length(timefreq_data.(datafield)(1).frequencies);
end

if parms.loadflag
  [path fname] = fileparts(timefreq_data.parms.filename);
  fname = [path '/' fname];
  timefreq_data.(datafield) = rmfield(timefreq_data.(datafield),...
                           intersect(fieldnames(timefreq_data.(datafield)),{'power','data','cross','cmplx'}));
  plv_data = timefreq_data;                                              
  clear timefreq_data
else
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

% calculate plv
for e = 1:length(events)
  if coyflg, coy = zeros(numcmb,ntime,nfreq); end
  if cohflg, coh = zeros(numcmb,ntime,nfreq); end
  if mscflg, msc = zeros(numcmb,ntime,nfreq); end  
  if plvflg
    if parms.band
      plv = zeros(numcmb,ntime);
    else
      plv = zeros(numcmb,ntime,nfreq); 
    end
  end
  plv_data.timefreq(e).labelcmb = [];
  for k=1:numcmb
    i = cmb(k,1);
    j = cmb(k,2);
    if parms.verbose
      fprintf('%s: calculating %s matrix %i of %i (%s,%s) for event %i\n',mfilename,parms.output,k,numcmb,sens(i).label,sens(j).label,events(e));
    end
    % use mscohere if given epoch_data
    % ...
    % continue;
    if parms.loadflag
      data = getdata(fname,events(e),i,j);
    else
      data = ts_data_selection(timefreq_data,'channels',[i j],'removebadchans',1);
      data = double(data.timefreq.cmplx);
      if i == j
        data(2,:,:,:) = data(1,:,:,:);
      end
    end
    ntrl = size(data,4);
    if plvflg
      phase1 = angle(squeeze(data(1,:,:,:)));
      phase2 = angle(squeeze(data(2,:,:,:)));
      clear data;
      if parms.band
        plv(k,:) = getPLV3D(phase1(:,parms.foi,:),phase2(:,parms.foi,:),parms.band)';
      else
        plv(k,:,:) = getPLV3D(phase1,phase2,parms.band);
      end
      clear phase1 phase2
    end
    if cohflg || coyflg || mscflg
      x   = squeeze(data(1,:,:,:)); % time x freq x trials
      y   = squeeze(data(2,:,:,:));
      sxy = sum(x.*conj(y),3) / size(x,3); % time x freq - cross-spectrum analogous to covariance
      pxx = sum(x.*conj(x),3) / size(x,3);
      pyy = sum(y.*conj(y),3) / size(y,3);
      if coyflg, coy(k,:,:) = sxy ./ sqrt(pxx .* pyy);      end   % coherency (complex)
      if cohflg, coh(k,:,:) = abs(sxy ./ sqrt(pxx .* pyy)); end
      if mscflg, msc(k,:,:) = abs(sxy).^2 ./ (pxx.*pyy);    end % magnitude squared coherence (real) - analogous to Pearson correlation r^2
    end
    plv_data.timefreq(e).labelcmb = [plv_data.timefreq(e).labelcmb; {sens(i).label sens(j).label}];
  end
  try save(sprintf('/home/halgdev/projects/jsherfey/plv/plv%g.mat',events(e)),'plv_data','plv'); end
  if plvflg, plv_data.timefreq(e).plv      = plv; end
  if coyflg, plv_data.timefreq(e).cohcmplx = coy; end
  if cohflg, plv_data.timefreq(e).coh      = coh; end
  if mscflg, plv_data.timefreq(e).mscoh    = msc; end    
	plv_data.timefreq(e).num_trials          = ntrl;
	clear plv coy coh msc
end

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

function data = backwardcompatible(data,parms)
if isfield(data,'opt')
  data.parms = data.opt;
  data = rmfield(data,'opt');
end
if ~issubfield(data,'parms.filename') && parms.loadflag
  error('timefreq_data filename not found.');
elseif parms.loadflag
  if iscell(data.parms.filename)
    data.parms.filename = data.parms.filename{1}; 
  end
end

function data = getdata(fname,e,i,j)
load(sprintf('%s_chan%03i.mat',fname,i));
data = double(timefreq_data.timefreq([timefreq_data.timefreq.event_code]==e).cmplx);
load(sprintf('%s_chan%03i.mat',fname,j));
data = cat(1,data,double(timefreq_data.timefreq([timefreq_data.timefreq.event_code]==e).cmplx));

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