function synth_data = rc_synth_sensors_from_RCSE(varargin)
%function synth_data = rc_synth_sensors_from_RCSE([options])
%
% Purpose: simulate sensor data from RCSE source waveforms
%
% Usage:
%  synth_data = rc_synth_sensors_from_RCSE('key1', value1,...);
%
% Optional Inupt:
%  'prefix': prefix of RCSE output files
%    {default = 'RCSE'}
%  'forward_prefix': prefix of RCSE output files different forward
%    If empty, use prefix
%    {default = []}
%  'rootdir': directory containing matfiles dir
%    {default = pwd}
%  'forward_rootdir':  directory containing matfiles dir for different forward
%    If empty, use rootdir
%    {default = []}
%  'fiterr_flag': which type of sensor data
%    0: RCSE fitted data plus noise
%    1: residual error of RCSE fit
%    2: RCSE fit projected into sensor data plus noise
%    3: data minus RCSE fit projected into sensor data
%     {default = 0}
%  'areas': when fiterr_flag=0, use only these area numbers
%    if empty, use all areas
%    {default = []}
%  'sourcefact': amplitude of modeled sources
%    {default = 1}
%  'baselineflag': use repeated baseline as noise added to synth_data
%    {default = 1}
%  'sources': matrix of source time courses
%     if empty, will use RCSE estimated source time courses
%     size must be [ntpoints,narea,ncontrasts]
%     only used if fiterr_flag = 0
%    {default = []}
%
% Output:
%   synth_data: synthesized data in avg_data structure
%
% NOTE: currently will not work if indy_locs_flag=1 or
%   loose_flag=1
%
% Created:  03/15/10 by Don Hagler
% Last Mod: 04/26/13 by Don Hagler
%

%% todo: option? to not resample to time_forward, instead, change
%%       synth_data.averages.data, stdev, and time
%%       need to also change forward.inv_scale_matrix

synth_data = [];

%if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'prefix','RCSE',[],...
  'forward_prefix',[],[],...
  'rootdir',pwd,[],...
  'forward_rootdir',[],[],...
  'fiterr_flag',0,[0:3],...
  'areas',[],[],...
  'sourcefact',1,[-Inf,Inf],...
  'baselineflag',true,[false true],...
  'sources',[],[],...
});
if isempty(parms.forward_prefix)
  parms.forward_prefix = parms.prefix;
end;
if isempty(parms.forward_rootdir)
  parms.forward_rootdir = parms.rootdir;
end;
orig_parms = parms;

% load parms
matfile=sprintf('%s/matfiles/%s_parms.mat',parms.rootdir,parms.prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);
fields = fieldnames(orig_parms);
for f=1:length(fields)
  parms = setfield(parms,fields{f},getfield(orig_parms,fields{f}));
end;

if ~strcmp(parms.forward_rootdir,parms.rootdir) ||...
   ~strcmp(parms.forward_prefix,parms.prefix)
  % load parms for forward
  matfile=sprintf('%s/matfiles/%s_parms.mat',...
    parms.forward_rootdir,parms.forward_prefix);
  if ~exist(matfile,'file')
    error('file %s not found',matfile);
  end;
  tmp = load(matfile);
  forward_parms = tmp.parms;
else
  forward_parms = parms;
end;

% load other mat files
flist = {'ret_forward','forward_prep','avg_data','results'};
if parms.fiterr_flag==1
  rootflags = [1 1 1 1];
else
  rootflags = [1 1 1 0];
end;
for f=1:length(flist)
  if rootflags(f)
    rootdir = parms.forward_rootdir;
    prefix = parms.forward_prefix;
  else
    rootdir = parms.rootdir;
    prefix = parms.prefix;
  end;
  matfile=sprintf('%s/matfiles/%s_%s.mat',rootdir,prefix,flist{f});
  if ~exist(matfile,'file')
    error('file %s not found',matfile);
  end;
  load(matfile);
end;

synth_data = avg_data;
num_conds = length(synth_data.averages);
% clear data matrices
for cond=1:num_conds
  synth_data.averages(cond).data = zeros(size(synth_data.averages(cond).data),'single');
end;
if isempty(parms.areas)
  parms.areas = [1:results.retmap.num_areas];
end;
parms.areas = intersect(parms.use_areas,parms.areas);

time_source = results.time;
time_forward = avg_data.averages(1).time;
num_sources = size(results.S,2);
num_sensors = length(forward_parms.goodchans);
num_tpoints = length(time_source);
num_tpoints_forward = length(time_forward);

switch parms.fiterr_flag
  case 0 % synthesize sensor data from retinotopy fit
    if ~isempty(parms.sources)
      time_source_orig = time_source;
      num_tpoints_orig = num_tpoints;
      num_sources_orig = num_sources;
      num_contrasts_orig = parms.ncontrasts;
      [num_tpoints,num_sources,num_contrasts] = size(parms.sources);
      if num_contrasts>1 && num_contrasts~=num_contrasts_orig
        error('number of elements in third dimension of sources should be %d',...
          num_contrasts_orig);
      end; 
      if num_sources>1 && num_sources~=num_sources_orig
        error('number of elements in second dimension of sources should be %d',...
          num_sources_orig);
      end;
      if num_tpoints>num_tpoints_orig
        dt = time_source_orig(2) - time_source_orig(1);
        nextra = num_tpoints - num_tpoints_orig;
        tlast = time_source_orig(end);
        tfirst_new = tlast + dt;
        tlast_new = tlast + nextra*dt;
        time_source = [time_source linspace(tfirst_new,tlast_new,nextra)];
      elseif num_tpoints<num_tpoints_orig
        time_source = time_source_orig(1:num_tpoints);
      end;
      sources = parms.sources;
      if num_sources==1 && num_contrasts==1
        sources = repmat(sources,[1,num_sources_orig,num_contrasts_orig]);
      elseif num_sources==1
        sources = repmat(sources,[1,num_sources_orig,1]);
      elseif num_contrasts==1
        sources = repmat(sources,[1,1,num_contrasts_orig]);
      end;
      [num_tpoints,num_sources,num_contrasts] = size(sources);
    else
      sources = results.S;
    end;
  
    % calculate fitted sensor waveforms for all conditions
    if forward_parms.ncontrasts~=parms.ncontrasts ||...
       forward_parms.nconds~=parms.nconds
      % average over multiple source contrast levels
      S = single(0);
      for i=1:parms.ncontrasts
        % select waveforms for this contrast
        tmpS = zeros(num_tpoints,num_sources);
        tmpS(:,parms.areas) = squeeze(parms.sourcefact*sources(:,parms.areas,i));
        S = S + tmpS;
      end;
      S = S/parms.ncontrasts;
      % use spline to adjust for difference in time points
      S = resample_sources(S,time_source,time_forward);
      % use average for multiple forward contrast levels
      Yfit = zeros(num_tpoints_forward,forward_parms.nconds*num_sensors);
      for i=1:forward_parms.ncontrasts
        % calculate fit
        tmpYfit = (retforward.F*S')';
        % insert result into Yfit
        j = 1 + (i-1)*length(forward_parms.unique_location_conds)*num_sensors;
        k = j + length(forward_parms.unique_location_conds)*num_sensors - 1;
        Yfit(:,j:k) = tmpYfit;
      end;
    else
      % loop over multiple contrast levels
      Yfit = zeros(num_tpoints_forward,parms.nconds*num_sensors);
      for i=1:parms.ncontrasts
        % select waveforms for this contrast
        S = zeros(num_tpoints,num_sources);
        S(:,parms.areas) = parms.sourcefact*sources(:,parms.areas,i);
        % use spline to adjust for difference in time points
        S = resample_sources(S,time_source,time_forward);
        % calculate fit
        tmpYfit = (retforward.F*S')';
        % insert result into Yfit
        j = 1 + (i-1)*length(parms.unique_location_conds)*num_sensors;
        k = j + length(parms.unique_location_conds)*num_sensors - 1;
        Yfit(:,j:k) = tmpYfit;
      end;
    end;

    % place in struct
    cond_nums = [forward_parms.cond_info.cond_number];
    for cond=1:num_conds
      matsize = size(synth_data.averages(cond).data);
      synth_data.averages(cond).data = zeros(matsize,'single');
      synth_data.averages(cond).stdev = zeros(matsize,'single');
      c = find(cond==cond_nums);
      if ~isempty(c)
        % get fitted data for this condition
        k1=1+length(forward_parms.goodchans)*(c-1);
        k2=k1+length(forward_parms.goodchans)-1;
        synth_data.averages(cond).data(forward_parms.goodchans,:)=Yfit(:,k1:k2)';
        % scale back to original units
        synth_data.averages(cond).data = ...
          synth_data.averages(cond).data.*forward.inv_scale_matrix;
      end;
    end;
  case 1 % use residual error from retinotopy fit
    for c=1:length(results.retmap.orig_cond_order)
      cond = results.retmap.orig_cond_order(c);
      k1=1+length(forward_parms.goodchans)*(c-1);
      k2=k1+length(forward_parms.goodchans)-1;
      synth_data.averages(cond).data(forward_parms.goodchans,:)=results.E(:,k1:k2)';
      synth_data.averages(cond).data = ...
        synth_data.averages(cond).data.*forward.inv_scale_matrix;
    end;
  case {2,3} % calculate projection of source waveforms into sensor data
    % calculate average source waveform across multiple contrast levels
    if forward_parms.ncontrasts~=parms.ncontrasts
      % average over multiple source contrast levels
      S = single(0);
      for i=1:parms.ncontrasts
        S = S + results.S(:,parms.areas,i);
      end;
      S = S/parms.ncontrasts;
      % use spline to adjust for difference in time points
      S = resample_sources(S,time_source,time_forward);
      if parms.fiterr_flag==2 % projection
        P = S*pinv(S);
      else % data minus projection
        P = eye(num_tpoints) - S*pinv(S);
      end;
      % project out of sensor data for each condition with this contrast
      for cond=1:length(synth_data.averages)
        tmp_data = synth_data.averages(cond).data(forward_parms.goodchans,:)';
        tmp_data = P*tmp_data;
        synth_data.averages(cond).data(forward_parms.goodchans,:) = tmp_data';
      end;
    else
      % loop over all conditions
      contrasts = [parms.orig_cond_info.contrast];
      uniq_contrasts = unique(contrasts(contrasts>0));
      for cond=1:num_conds
        contrast = contrasts(cond);
        if contrast==0, continue; end;
        ind_cont = find(contrast == uniq_contrasts);
        S = results.S(:,parms.areas,ind_cont);
        % use spline to adjust for difference in time points
        S = resample_sources(S,time_source,time_forward);
        if parms.fiterr_flag==2 % projection
          P = S*pinv(S);
        else % data minus projection
          P = eye(num_tpoints) - S*pinv(S);
        end;
        tmp_data = avg_data.averages(cond).data(forward_parms.goodchans,:)';
        tmp_data = P*tmp_data;
        synth_data.averages(cond).data(forward_parms.goodchans,:) = tmp_data';
      end;
    end;
end;

% add noise to all goodchans for all conditions
if parms.baselineflag && ismember(parms.fiterr_flag,[0,2])
  for cond=1:num_conds
    time = avg_data.averages(cond).time;
    ind_baseline = find(time<0);
    baseline_noise = ...
      avg_data.averages(cond).data(forward_parms.goodchans,ind_baseline);
    num_repeats = ceil(length(time)/length(ind_baseline));
    full_noise = zeros(size(avg_data.averages(cond).data),'single');
    j = 0;
    for i=1:num_repeats
      if ~mod(i,2) % even, keep same
        tmp_noise = baseline_noise;
      else % off, flip time
        tmp_noise = flipdim(baseline_noise,2);
      end;
      if length(ind_baseline) > (length(time) - j)
        num_samples = length(time) - j;
      else
        num_samples = length(ind_baseline);
      end;
      full_noise(forward_parms.goodchans,j+1:j+num_samples) = ...
        tmp_noise(:,1:num_samples);
      j = j + num_samples;
    end;
    synth_data.averages(cond).data = ...
      synth_data.averages(cond).data + full_noise;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S_res = resample_sources(S,time,time_res)
  if length(time)~=length(time_res) || any(time ~= time_res)
    num_tpoints = length(time_res);
    num_sources = size(S,2);
    dt = time_res(2) - time_res(1);
    S_res = zeros(num_tpoints,num_sources);
    ind_pre = find((time(1) - time_res)>dt);
    ind_post = find((time(end) - time_res)<-dt);
    for s=1:num_sources
      S_res(:,s) = spline(time,S(:,s),time_res);
      S_res(ind_pre,s) = 0;
      S_res(ind_post,s) = 0;
    end;
  else
    S_res = S;
  end;
return;

