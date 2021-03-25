function prior = rc_extract_prior(fname_prior,time,areas,contrasts,avg_flag)
%function prior = rc_extract_prior(fname_prior,[time],[areas],[contrasts],[avg_flag])
%
% Required Input:
%   fname_prior: full path name of mat file containing
%     'results' struct that can be from RCSE or Group RCSE
%     May also directly supply results struct
%
% Optional Input:
%   time: vector of times for output waveform S_ref
%     if empty, S_ref will be a copy of S
%     {default = []}
%   areas: vector of of visual areas included in source estimates
%     if empty, will use all areas in prior
%     {default = []}
%   contrasts: vector of contrast levels included in source estimates
%     if empty, will use all contrast levels in prior
%     {default = []}
%   avg_flag: [0|1|2] whether to use averaged prior
%     0: no averaging -- requires that prior and fitted waveform have
%          same number of visual areas and contrast levels
%     1: average prior across contrast levels, but not visual areas
%     2: average prior across visual areas, but not contrast levels
%     3: average prior across contrast levels and visual areas
%     {default = 0}
%
% Output:
%   prior: struct containing these fields:
%     S: original prior source estimates 
%     S_ref: source estimates resampled to reference time series
%     S_sem: original prior standard error of mean source estimates
%     S_sem_ref: SEM resampled to reference time series
%     time: copy of original time in prior results
%     time_ref: time vector for S_ref
%
%   NOTE: S_sem and S_sem_ref will be empty unless S_sem found in prior
%     e.g. if resamp_flag=1 when MEG_MMIL_Group_RCSE was run
%     
% Created:  03/09/11 by Don Hagler
% Last Mod: 03/05/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('time','var'), time=[]; end;
if ~exist('areas','var'), areas=[]; end;
if ~exist('contrasts','var'), contrasts=[]; end;
if ~exist('avg_flag','var') || isempty(avg_flag), avg_flag=0; end;
prior = [];

if isstruct(fname_prior)
  results = fname_prior;
elseif ~exist(fname_prior,'file')
  error('prior file %s not found',fname_prior);
else
  results = [];
  load(fname_prior);
end;
if isempty(results)
  error('prior file %s has missing or empty results struct',fname_prior);
end;

% get information about prior waveforms
prior.S = results.S;
prior.S_sem = mmil_getfield(results,'S_sem',[]);
prior.time = results.time;
prior.ntpoints = length(prior.time);
prior.areas = results.areas;
prior.nareas = results.nareas;
prior.contrasts = results.contrasts;
prior.ncontrasts = results.ncontrasts;
prior.nwforms = prior.nareas*prior.ncontrasts;

% set information about prior_ref waveforms (to be used as reference)
if isempty(time)
  prior.time_ref = prior.time;
else
  prior.time_ref = time;
end;
prior.ntpoints_ref = length(prior.time_ref);
if isempty(areas)
  prior.areas_ref = prior.areas;
  prior.nareas_ref = prior.nareas;
else
  prior.areas_ref = areas;
  prior.nareas_ref = length(areas);
end;
if isempty(contrasts)
  prior.contrasts_ref = prior.contrasts;
  prior.ncontrasts_ref = prior.ncontrasts;
else
  prior.contrasts_ref = contrasts;
  prior.ncontrasts_ref = length(contrasts);
end;
prior.nwforms_ref = prior.nareas_ref * prior.ncontrasts_ref;

% check that input areas and contrasts are found in prior
ind_areas = find(ismember(prior.areas_ref,prior.areas));
ind_contrasts = find(ismember(prior.contrasts,prior.contrasts_ref));
if length(ind_areas)<prior.nareas_ref && ~ismember(avg_flag,[2,3])
  error('input areas not found in prior');
end;
if length(ind_contrasts)<prior.ncontrasts_ref && ~ismember(avg_flag,[1,3])
  error('input contrast levels not found in prior');
end;

% reduce prior to selected areas and contrast levels
if ~isempty(ind_areas) && prior.nareas_ref ~= prior.nareas
  prior.S = prior.S(:,ind_areas,:);
  prior.areas = prior.areas_ref;
  prior.nareas = prior.nareas_ref;
  if ~isempty(prior.S_sem)
    prior.S_sem = prior.S_sem(:,ind_areas,:);
  end;
end;
if ~isempty(ind_contrasts) && prior.ncontrasts_ref ~= prior.ncontrasts
  prior.S = prior.S(:,:,ind_contrasts);
  prior.contrasts = prior.contrasts_ref;
  prior.ncontrasts = prior.ncontrasts_ref;
  if ~isempty(prior.S_sem)
    prior.S_sem = prior.S_sem(:,:,ind_contrasts);
  end;
end;

% average across contrasts or areas or both
%  to avoid biasing between-area latency differences and
%  between-contrast amplitude differences
size_S = [prior.ntpoints,prior.nareas_ref,prior.ncontrasts_ref];
switch avg_flag
  case 0 % no averaging -- just check that number of wforms match
    if prior.nareas_ref ~= prior.nareas
      error('number of areas in source estimate matrix (%d) does not match prior (%d) from %s',...
        prior.nareas_ref,prior.nareas,fname_prior); % maybe redundant
    elseif prior.ncontrasts_ref ~= prior.ncontrasts
      error('number of contrast levels in source estimate matrix (%d) does not match prior (%d) from %s',...
        prior.ncontrasts_ref,prior.ncontrasts,fname_prior); % maybe redundant
    end;
  case 1 % average prior across contrast levels, but not visual areas
    if prior.nareas_ref ~= prior.nareas
      error('number of areas in source estimate matrix (%d) does not match prior (%d) from %s',...
        prior.nareas_ref,prior.nareas,fname_prior); % maybe redundant
    end;
    tmpS = mean(prior.S,3);
    % set waveform for each contrast to average
    prior.S = zeros(size_S);
    for c=1:prior.ncontrasts_ref
      prior.S(:,:,c) = tmpS;
    end;
    if ~isempty(prior.S_sem)
      tmpS = mean(prior.S_sem,3);
      % set waveform for each contrast to average
      prior.S_sem = zeros(size_S_sem);
      for c=1:prior.ncontrasts_ref
        prior.S_sem(:,:,c) = tmpS;
      end;
    end;
  case 2 % average prior across visual areas, but not contrast levels
    if prior.ncontrasts_ref ~= prior.ncontrasts
      error('number of contrast levels in source estimate matrix (%d) does not match prior (%d) from %s',...
        prior.ncontrasts_ref,prior.ncontrasts,fname_prior); % maybe redundant
    end;
    tmpS = mean(prior.S,2);
    % set waveform for each area to average
    prior.S = zeros(size_S);
    for s=1:prior.nareas_ref
      prior.S(:,s,:) = tmpS;
    end;
    if ~isempty(prior.S_sem)
      tmpS = mean(prior.S_sem,2);
      % set waveform for each area to average
      prior.S_sem = zeros(size_S);
      for s=1:prior.nareas_ref
        prior.S_sem(:,s,:) = tmpS;
      end;
    end;
  case 3 % average prior across contrast levels and visual areas
    tmpS = mean(reshape(prior.S,...
      [prior.ntpoints,prior.nareas*prior.ncontrasts]),2);
    % set waveform for each area and contrast to average
    prior.S = zeros(size_S);
    for s=1:prior.nareas_ref
      for c=1:prior.ncontrasts_ref
        prior.S(:,s,c) = tmpS;
      end;
    end;
    if ~isempty(prior.S_sem)
      tmpS = mean(reshape(prior.S_sem,...
        [prior.ntpoints,prior.nareas*prior.ncontrasts]),2);
      % set waveform for each area and contrast to average
      prior.S_sem = zeros(size_S);
      for s=1:prior.nareas_ref
        for c=1:prior.ncontrasts_ref
          prior.S_sem(:,s,c) = tmpS;
        end;
      end;
    end;    
end;

% resample prior waveforms to match time samples
size_S_ref = [prior.ntpoints_ref,prior.nareas_ref,prior.ncontrasts_ref];
prior.S_ref = zeros(size_S_ref);
ind_pre = find((prior.time(1) - prior.time_ref)>0);
ind_post = find((prior.time(end) - prior.time_ref)<0);
for s=1:prior.nwforms_ref
  prior.S_ref(:,s) = spline(prior.time,prior.S(:,s),prior.time_ref);
  prior.S_ref(ind_pre,s) = 0;
  prior.S_ref(ind_post,s) = 0;
end;

if ~isempty(prior.S_sem)
  prior.S_sem_ref = zeros(size_S_ref);
  ind_pre = find((prior.time(1) - prior.time_ref)>0);
  ind_post = find((prior.time(end) - prior.time_ref)<0);
  for s=1:prior.nwforms_ref
    prior.S_sem_ref(:,s) = spline(prior.time,prior.S_sem(:,s),prior.time_ref);
    prior.S_sem_ref(ind_pre,s) = 0;
    prior.S_sem_ref(ind_post,s) = 0;
  end;
else
  prior.S_sem_ref = [];
end;

clear results
