function syncdata = ts_sync(varargin)
% Purpose: calculate bivariate synchrony metrics (ex. plv, coherence)
% 
% Inputs:
%  data: either -
%   - timefreq_data, [cmplx] = chan x time x freq x trial
%   - cell array of file names with timefreq_data trial spectra
%   - struct with data.parms.filename containing a list of file names
%
% Undocumented options: trialdim

% check inputs & load data if necessary
data = [];
if mod(nargin,2)      % odd number of inputs
  data = varargin{1};
  if nargin > 1
    varargin = varargin(2:end);
  else
    varargin = {};
  end
  if ~(issubfield(data,'timefreq.cmplx') && ndims(data.timefreq(1).cmplx)==4)
    try     % load TF trial spectra
      data = ts_load_timefreq_data(data,varargin{:});
    catch   % load average TF spectra
      data = ts_checkdata_header(data,varargin{:});
    end
  end
end
% => [tfr] = chan x time x freq x trial  ||  chan x time x freq

% Set up parameters
parms = mmil_args2parms(varargin,...
						{'foi',[],[],...
             'foilim',[],[],...
             'toi',[],[],...
             'toilim',[],[],...
             'events',[],[],...
             'chanlabels',[],[],...
             'chantype',[],[],...
						 'channelcmb','all',[],...
             'output','plv',[],...
             'exc_freqs',[],[],...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...    
             'trialdim',4,[],...
             'rejectfile',[],[],...
             'reject_data',[],[],...
             'prescale_mag'            ,10^15,[],...  % to fT
             'prescale_grad'           ,10^13,[],... % to fT/cm
             'prescale_eeg'            ,1,[],...   % to uV
             'prescale_eog'            ,1,[],...   % to uV
             'refchan','all',[],...
						},false);
          
% Backwards compatibility
parms = backcompatible(parms,varargin{:});     

% Check parameters for consistency
if ~iscell(parms.output), parms.output = {parms.output}; end
  
% Select data to process
data  = ts_checkdata_header(data,'precision','single','events',parms.events,'verbose',0);
data  = ts_data_selection(data,'chanlabels',parms.chanlabels,'chantype',parms.chantype,'toi',parms.toi,'toilim',parms.toilim,...
  'events',parms.events,'foi',parms.foi,'foilim',parms.foilim,'reject_data',parms.reject_data,'rejectfile',parms.rejectfile); 

% Get info on data
[datatype datafield dataparam] = ts_object_info(data);
sens  = data.sensor_info;
T     = data.(datafield)(1).time;
F     = data.(datafield)(1).frequencies;
Fix   = find(~ismember(F,parms.exc_freqs));
ncond          = length(data.(datafield));
nchan          = data.num_sensors;
ntime          = length(T);
nfreq          = length(Fix);
ntrials        = [data.(datafield).num_trials];
parms.channels = 1:data.num_sensors;
trialdim       = parms.trialdim;

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

% Reference channel specification
% if isempty(parms.refchan) || any(strcmpi(parms.refchan,'all'))
%   allrefchans = {data.sensor_info.label};
% elseif isnumeric(parms.refchan)
%   allrefchans = {data.sensor_info(parms.refchan).label};
% elseif ischar(parms.refchan)
%   allrefchans = {parms.refchan};
% elseif iscellstr(parms.refchan)
%   allrefchans = parms.refchan;
% else
%   error('Reference channels not specified.');
% end
% if ~any(ismember(allrefchans,{data.sensor_info.label}))
%   error('The reference channel(s) were not found in the data.');
% end
% [sel1 sel2] = match_str({data.sensor_info.label},allrefchans);
% cmb   = cmb(ismember(cmb(:,1),sel1) | ismember(cmb(:,2),sel1),:);
% nchan = length(sel1);

% cmb is now a matrix of unique channel index pairs
numcmb = size(cmb,1);

if ismember('plv',      parms.output), plvflg = 1; else plvflg = 0; end   % phase locking value
if ismember('coherency',parms.output), coyflg = 1; else coyflg = 0; end   % coherency (complex)
if ismember('coh',      parms.output), cohflg = 1; else cohflg = 0; end   % coherence |complex|
if ismember('mscoh',    parms.output), mscflg = 1; else mscflg = 0; end   % magnitude squared coherence (real) - analogous to Pearson correlation r^2
if ismember('lag',      parms.output), lagflg = 1; else lagflg = 0; end   % phase lag

if ~issubfield(data,'timefreq.cmplx')
  error('complex spectra must be provided for calculating synchrony metrics.');
else
  syncdata  = rmsubfield(data,'timefreq.cmplx');
  if isfield(syncdata.timefreq,'power'), syncdata.timefreq = rmfield(syncdata.timefreq,'power'); end
end

% prescale MEG data to avoid rounding spectra < eps to zero
typestring = {data.sensor_info.typestring};
mag_i   = strcmp ('mag' ,typestring);
grad_i  = strncmp('grad',typestring,4);
eeg_i   = strcmp ('eeg' ,typestring);
eog_i   = strcmp ('eog' ,typestring);
for c = 1:ncond
  data.timefreq(c).cmplx(mag_i,:,:,:)  = data.timefreq(c).cmplx(mag_i,:,:,:) * parms.prescale_mag;
  data.timefreq(c).cmplx(grad_i,:,:,:) = data.timefreq(c).cmplx(grad_i,:,:,:)* parms.prescale_grad;
  data.timefreq(c).cmplx(eeg_i,:,:,:)  = data.timefreq(c).cmplx(eeg_i,:,:,:) * parms.prescale_eeg;
  data.timefreq(c).cmplx(eog_i,:,:,:)  = data.timefreq(c).cmplx(eog_i,:,:,:) * parms.prescale_eog;
end;

% calculate synchrony metrics from observed data
% preallocate synchrony matrices
for c = 1:ncond
  if plvflg, syncdata.timefreq(c).plv   = zeros(numcmb,ntime,nfreq,'single'); end
  if coyflg, syncdata.timefreq(c).coy   = zeros(numcmb,ntime,nfreq,'single'); end
  if cohflg, syncdata.timefreq(c).coh   = zeros(numcmb,ntime,nfreq,'single'); end
  if mscflg, syncdata.timefreq(c).mscoh = zeros(numcmb,ntime,nfreq,'single'); end  
  if lagflg, syncdata.timefreq(c).lag   = zeros(numcmb,ntime,nfreq,'single'); end  
end  
% loop over conditions and calculate metrics
for c = 1:ncond
  % for each channel, calculate metrics for all it's remaining pairs
  for k = 1:nchan
    if parms.verbose, fprintf('channel %g of %g (%s)\n',k,nchan,[parms.output{:}]); end
    ind = cmb(:,1)==k;
    sel = cmb(ind,:); % select all remaining channel pairs for this channel
    x   = double(data.timefreq(c).cmplx(sel(:,1),:,Fix,:)); % TF spectrum, sensor 1
    y   = double(data.timefreq(c).cmplx(sel(:,2),:,Fix,:)); % TF spectrum, sensor 2
    if cohflg || coyflg || mscflg || lagflg
      sxy = sum(x.*conj(y),trialdim) / size(x,trialdim); % time x freq - cross-spectrum analogous to covariance
      pxx = sum(x.*conj(x),trialdim) / size(x,trialdim);
      pyy = sum(y.*conj(y),trialdim) / size(y,trialdim);
      if coyflg, syncdata.timefreq(c).coy(ind,:,:)   = single(sxy ./ sqrt(pxx .* pyy));      end   % coherency (complex)
      if cohflg, syncdata.timefreq(c).coh(ind,:,:)   = single(abs(sxy ./ sqrt(pxx .* pyy))); end
      if mscflg, syncdata.timefreq(c).mscoh(ind,:,:) = single(abs(sxy).^2 ./ (pxx.*pyy));    end % magnitude squared coherence (real) - analogous to Pearson correlation r^2
      if lagflg, syncdata.timefreq(c).lag(ind,:,:)   = single(atan(imag(sxy) ./ real(sxy))); end % phase difference, [-pi,pi]
      clear sxy pxx pyy
    end
    if plvflg
      X   = x.*conj(y);                               % cross-spectrum
      X   = X ./ abs(X);                              % normalized cross-spectrum (= exp(i*(relative phase)))      
    syncdata.timefreq(c).plv(ind,:,:) = single(abs(sum(X,trialdim) / size(X,trialdim)));
    clear X
    end
    clear x y sel ind
  end
  syncdata.timefreq(c).labelcmb = [{data.sensor_info(cmb(:,1)).label}' {data.sensor_info(cmb(:,2)).label}'];
end
fprintf('success!\n');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parms = backcompatible(parms,varargin)

opt = mmil_args2parms(varargin,{...
        'event_codes',[],[],...
        'blc',[],[],...
        },false);
  
if isempty(parms.events) && ~isempty(opt.event_codes)
  parms.events = opt.event_codes;
end
if ~isempty(opt.blc)
  if isequal(opt.blc,false)
    parms.blc = 'no';
  elseif isequal(opt.blc,true)
    parms.blc = 'yes';
  end
end
