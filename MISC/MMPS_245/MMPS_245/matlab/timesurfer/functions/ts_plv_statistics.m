function [stat_data,stats,pvals] = ts_plv_statistics(varargin)


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
else
  data = ts_load_timefreq_data(varargin{:});
end
% => [tfr] = chan x time x freq x trial  ||  chan x time x freq

% Set up parameters
parms = mmil_args2parms(varargin,...
						{'foi',[],[],...
             'foilim',[],[],...
             'toi',[],[],...
             'toilim',[],[],...
             'freqbins',[],[],...
             'timebins',[],[],...
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
             'numrandomization',500,[],...
             'method','montecarlo',{'montecarlo','ttest'},...
						},false);
          
% Backwards compatibility
parms = backcompatible(parms,varargin{:});     

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

if ncond > 2
  error('This function only supports comparing two conditions at a time.');
end
% check time and frequency bins
if isempty(parms.freqbins), parms.freqbins{1} = [F(1) F(end)]; end
if isempty(parms.timebins), parms.timebins{1} = [T(1) T(end)]; end
if isnumeric(parms.freqbins), parms.freqbins = {parms.freqbins}; end
if isnumeric(parms.timebins), parms.timebins = {parms.timebins}; end
nbins_freq = numel(parms.freqbins);
nbins_time = numel(parms.timebins);

cmb = [];
% channel pairs for which to calc PLV stats
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

if strcmpi(parms.method,'montecarlo')
  output_nfreq = 1;
  output_ntime = length(T);
elseif strcmpi(parms.method,'ttest')
  output_nfreq = length(parms.freqbins);
  output_ntime = length(parms.timebins);
else
  error('blah');
end

% preallocate variables
stats = zeros(numcmb,output_ntime,output_nfreq,'single');
pvals = inf  (numcmb,output_ntime,output_nfreq,'single');
for c = 1:ncond
  testdata{c} = zeros(numcmb,ntime,nfreq,ntrials(c),'single');
end

% for each channel, calculate metrics for all it's remaining pairs
for k = 1:nchan
  if parms.verbose, fprintf('channel %g of %g\n',k,nchan); end
  ind = cmb(:,1)==k;
  sel = cmb(ind,:); % select all remaining channel pairs for this channel
  % loop over conditions and calculate mean relative phase over both conds
  MeanRelativePhase = zeros(length(find(ind)),ntime,nfreq,'double');
  for c = 1:ncond    
    x   = double(data.timefreq(c).cmplx(sel(:,1),:,Fix,:)); % TF spectrum, sensor 1
    y   = double(data.timefreq(c).cmplx(sel(:,2),:,Fix,:)); % TF spectrum, sensor 2
    X   = x.*conj(y);                                       % cross-spectrum
    MeanRelativePhase = sum(cat(4,MeanRelativePhase,X./abs(X)),4);
  end
  clear x y X
  MeanRelativePhase = angle(MeanRelativePhase/sum(ntrials));
  % loop over conditions and calculate pre-stat metric for each condition
  for c = 1:ncond
    x   = double(data.timefreq(c).cmplx(sel(:,1),:,Fix,:)); % TF spectrum, sensor 1
    y   = double(data.timefreq(c).cmplx(sel(:,2),:,Fix,:)); % TF spectrum, sensor 2    
    theta = angle(x) - angle(y);
    clear x y
    for trl = 1:ntrials(c)
      testdata{c}(ind,:,:,trl) = single(cos(MeanRelativePhase - theta(:,:,:,trl)));
    end
  end
  clear theta
end
labelcmb = [{data.sensor_info(cmb(:,1)).label}' {data.sensor_info(cmb(:,2)).label}'];
clear MeanRelativePhase
tic
switch lower(parms.method)
  case {'montecarlo','mc','mcstats','nonparametric'}
    % prepare test data structure
    epoch_data = rmfield(data,datafield);
    epoch_data.epochs = rmfield(data.(datafield),intersect(fieldnames(data.(datafield)),{'power','cmplx','frequencies'}));
    % prepare channel labels
    epoch_data.sensor_info(1:numcmb) = epoch_data.sensor_info(1);
    for k = 1:numcmb
      epoch_data.sensor_info(k).label = sprintf('%s-%s',labelcmb{k,1},labelcmb{k,2});
    end
    epoch_data.num_sensors = numcmb;
    clear data
    fprintf('Calculating Monte Carlo stats on %g channel pairs, %g freqs and %g times\n',numcmb,nbins_freq,ntime);
    % loop over frequency bins
    foi = F(Fix);
    for k = 1:length(parms.freqbins)
      fix = foi>=parms.freqbins{k}(1) & foi<=parms.freqbins{k}(2);
      if isempty(find(fix))
        fprintf('skipping freqbin %g of %g\n',k,nbins_freq);
        continue
      else
%         fprintf('freqbin %g of %g\n',k,nbins_freq);      
      end
      for c = 1:ncond
        epoch_data.epochs(c).data = squeeze(mean(testdata{c}(:,:,fix,:),3)); % chan pair x time x trials
      end
%       stat_data = ts_statistics_wrapper(epoch_data,'events',[epoch_data.epochs.event_code]);
      stat_data = ts_statistics(epoch_data,'events',[epoch_data.epochs.event_code],'numrandomization',parms.numrandomization);
      stats(:,:,k) = stat_data.stats(1).stat;
      pvals(:,:,k) = stat_data.stats(1).prob;
    end % end loop over frequency bins
    clear testdata
  case {'ttest','parametric'}
    clear data
    fprintf('Calculating t-stats on %g channel pairs, %g freq bins and %g time bins\n',numcmb,nbins_freq,nbins_time);
    % loop over frequency bins
    foi = F(Fix);
    for k = 1:length(parms.freqbins)
      fix = foi>=parms.freqbins{k}(1) & foi<=parms.freqbins{k}(2);
      if isempty(find(fix))
        fprintf('skipping freqbin %g of %g\n',k,nbins_freq);
        continue
      else
        fprintf('freqbin %g of %g\n',k,nbins_freq);      
      end
      % loop over time bins
      for j = 1:length(parms.timebins)
%         tix = T>=parms.timebins{j}(1) & T<=parms.timebins{j}(2);
        tix = nearest(T,parms.timebins{j}(1)):nearest(T,parms.timebins{j}(2));
        if (parms.timebins{j}(1)<T(1)) || (parms.timebins{j}(2)>T(end)) || isempty(tix)
          fprintf('skipping freqbin %g of %g: timebin %g of %g\n',k,nbins_freq,j,nbins_time);
          continue
        end
        for c = 1:ncond
          tfoi_data{c} = squeeze(mean(mean(testdata{c}(:,tix,fix,:),2),3)); % chan pair x time x trials
        end
        % loop over channel pairs
        for c = 1:numcmb
          x1 = tfoi_data{1}(c,:);
          x2 = tfoi_data{2}(c,:);
          [H,P,CI,STATS] = ttest2(x1,x2,.05,'both','equal');
          pvals(c,j,k) = P;
          stats(c,j,k) = STATS.tstat;
        end % end loop over channel pairs
      end % end loop over time bins
    end % end loop over freq bins
    stat_data = [];
    clear testdata
  otherwise
    error('unknown method requested.');
end
toc
% [stats|pvals] = cmb x ntime (x Fbin)
% [stats|pvals] = cmb x Tbin   x Fbin


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