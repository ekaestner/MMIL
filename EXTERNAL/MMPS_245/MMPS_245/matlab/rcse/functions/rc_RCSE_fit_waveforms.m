function [parms,results,retforward,inverse]=...
  rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward,retforward)
%function [parms,results,retforward,inverse]=...
%  rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward,[retforward])
%
% Purpose: calculate RCSE waveforms
%
% Required Input:
%   parms: RCSE parms struct
%   retmap: retinotopy mapping struct
%   avg_data: average data struct
%   forward: forward solution struct
%
% Optional Input:
%   retforward: retinotopy constrained forward struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/05/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,4), return; end;

% construct retintotopy constrained forward solution
if ~exist('retforward','var') || isempty(retforward)
  fprintf('%s: constructing retinotopy mapping matrix...\n',mfilename);
  % use unique_location_conds only (e.g. multiple contrast levels)
  tmp_event_codes = []; tmp_conditions = [];
  if ~isempty(parms.event_codes)
    tmp_event_codes = parms.event_codes(parms.unique_location_conds);
  end;
  if ~isempty(parms.conditions)
    tmp_conditions = parms.conditions(parms.unique_location_conds);
  end;
  %% todo: should not be necessary if M already exists
  %%        and if M was updated during rc_choose_retdips

  retforward.retmap = rc_construct_ret_mapping(retmap,...
    'use_areas',parms.use_areas,'conditions',tmp_conditions,...
    'event_codes',tmp_event_codes,parms.ret_mapping_opts{:});
  if isempty(retforward.retmap)
    error('failed to create retinotopy mapping matrix');
  end;

  % construct retinotopy forward matrix
  fprintf('%s: constructing retinotopy-constrained forward matrix...\n',mfilename);
  tmp_uniq_verts_lh = retforward.retmap.uniq_verts_lh;
  tmp_uniq_verts_rh = retforward.retmap.uniq_verts_rh;
  retforward.retmap.orig_uniq_verts_lh = forward.lh_verts;
  retforward.retmap.orig_uniq_verts_rh = forward.rh_verts;
  tic;
  retforward.F = rc_construct_ret_forward(...
    forward.G_norm,forward.G_tang1,forward.G_tang2,retforward.retmap,...
    parms.num_dips_lh,parms.num_dips_rh,parms.indy_locs_flag,parms.loose_flag,...
    parms.near_nbr_weight);
  toc;
  if isempty(retforward.F)
    error('failed to calculate retinotopy forward matrix');
  end;
  retforward.retmap.uniq_verts_lh = tmp_uniq_verts_lh;
  retforward.retmap.uniq_verts_rh = tmp_uniq_verts_rh;
  clear tmp_uniq_verts_lh tmp_uniq_verts_rh;
end;
[num_measurements,num_sources] = size(retforward.F);
num_sensors = length(parms.goodchans);

if parms.condF_thresh>0
  if isempty(parms.condF)
    num_locs = length(parms.unique_location_conds);
    parms.condF = zeros(num_locs,1);
    for i=1:num_locs
      j = (i-1)*num_sensors + 1;
      k = j + num_sensors - 1;
      tmpF = retforward.F(j:k,:);
      parms.condF(i) = cond(tmpF'*tmpF);
    end;
  end;
  if isempty(parms.cond_weights)
    parms.cond_weights = 1.0*(parms.condF>parms.condF_thresh);    
  elseif parms.reweight_flag
    parms.cond_weights = parms.cond_weights.*(parms.condF>parms.condF_thresh);    
  end;
  fprintf('%s: %d stimulus locations with condition # <= %0.1f\n',...
    mfilename,length(find(parms.cond_weights)),parms.condF_thresh);
end;

% store original data and forward solution
if ~isempty(parms.cond_weights)
  orig_retforward = retforward;
  orig_avg_data = avg_data;
end;

fprintf('%s: initializing output...\n',mfilename);
% initialize output
time = avg_data.averages(1).time;
ntpoints = length(time);
nmeas = parms.nconds*num_sensors;
results.time = time;
results.Y = zeros(ntpoints,nmeas);
results.Yfit = zeros(ntpoints,nmeas);
results.S = zeros(ntpoints,num_sources,parms.ncontrasts);

% extract prior waveform from mat file and resample to desired time
%   and number of waveforms
%   (e.g. to go from multiple contrasts to a single waveform)
if isempty(parms.fname_prior)
  parms.prior = [];
elseif ~isfield(parms,'prior') || isempty(parms.prior)
  fprintf('%s: extracting prior...\n',mfilename);
  tic;
  parms.prior = rc_extract_prior(parms.fname_prior,...
    time,parms.use_areas,parms.unique_contrasts,parms.prior_avg_flag);
  toc;
end;

% calculate inverse operator
switch parms.inverse_type
  case 0 % unregularized pseudo-inverse
    fprintf('%s: calculating simple inverse...\n',mfilename);
    inverse.W = rc_calc_simple_inverse(retforward.F);
  case 1 % regularized pseudo-inverse with identity matrix for noise covar
    fprintf('%s: calculating inverse with no source covar...\n',mfilename);
    inverse.W = rc_calc_ret_inverse_operator_noprior(retforward.F,parms.SNR);
  case 2 % regularized pseudo-inverse with noise and source covar
    % construct P and R source covariance matrices
    if isfield(parms,'P') && isfield(parms,'R')
      inverse.P = parms.P;
      inverse.R = parms.R;
    else
      fprintf('%s: constructing source covariance matrix...\n',mfilename);
      [inverse.P,inverse.R] = rc_calc_ret_sourcecov(...
        retforward.retmap,...
        'indy_locs_flag',parms.indy_locs_flag,...
        'theta_smfact',parms.indy_smfact,...
        'ecc_smfact',parms.ecc_smfact,...
        'upperlower_smfact',parms.upperlower_smfact,...
        'hemi_smfact',parms.hemi_smfact,...
        'loose_flag',parms.loose_flag,...
        'loose_tang_weight',parms.loose_tang_weight,...
        'visual_area_weight',parms.visual_area_weight,...
        'ret_dips_weight',parms.ret_dips_weight,...
        'nonret_dips_weight',parms.nonret_dips_weight);
      parms.R = inverse.R;
      parms.P = inverse.P;
    end;

    % construct C noise covariance matrix
    if isfield(parms,'C')
      inverse.C = parms.C;
    else
      fprintf('%s: constructing sensor noise covariance matrix...\n',mfilename);
      inverse.C = sparse(num_measurements,num_measurements);
      for i=1:retforward.retmap.num_locs
        j = 1 + (i-1)*num_sensors;
        k = j + num_sensors - 1;
        inverse.C(j:k,j:k) = parms.noisecovar;
      end;
      parms.C = inverse.C;
    end;

    fprintf('%s: calculating inverse with source covar...\n',mfilename);
    inverse.W = rc_calc_ret_inverse_operator(retforward.F,parms.SNR,...
      inverse.R,inverse.C);
end;

if ~isempty(parms.cond_weights)
  % adjust forward and data
  fprintf('%s: updating forward and data...\n',mfilename);
  retforward = orig_retforward;
  avg_data = orig_avg_data;
  % scale data
  for j=1:retforward.retmap.orig_num_locs
    i = parms.same_location_conds(j);
    c = retforward.retmap.orig_cond_order(j);
    avg_data.averages(c).data = ...
      parms.cond_weights(i) * avg_data.averages(c).data;
  end;

  % scale forward matrix
  for i=1:retforward.retmap.num_locs
    j = 1 + (i-1)*num_sensors;
    k = j + num_sensors - 1;
    retforward.F(j:k,:) = parms.cond_weights(i) * retforward.F(j:k,:);
  end;
end;

% loop over multiple contrast levels
for i=1:parms.ncontrasts
  % apply inverse to data
  fprintf('%s: applying inverse...\n',mfilename);

  % get conditions for this contrast level
  if ~isempty(parms.cond_info)
    ind_contrast = find(parms.contrasts == parms.unique_contrasts(i));
    cond_numbers = [parms.cond_info.cond_number];
    cond_order = cond_numbers(ind_contrast);
  else
    cond_order = retforward.retmap.cond_order;
  end;

  % load data
  tmpY = rc_load_ret_avgt(avg_data,cond_order,...
    parms.goodchans,forward.scale_matrix);
  if isempty(tmpY)
    error('failed to load data');
  end;

  % calculate source estimates
  tmpS = (inverse.W*tmpY')';

  % load original (unweighted) data
  if ~isempty(parms.cond_weights)
    tmpY = rc_load_ret_avgt(orig_avg_data,cond_order,...
      parms.goodchans,forward.scale_matrix);
    if isempty(tmpY)
      error('failed to load data');
    end;
    % revert to original (unweighted) forward solution
    retforward = orig_retforward;
  end;

  % calculate fit
  tmpYfit = (retforward.F*tmpS')';

  % store results
  j = 1 + (i-1)*length(parms.unique_location_conds)*num_sensors;
  k = j + length(parms.unique_location_conds)*num_sensors - 1;

  results.Y(:,j:k) = tmpY;
  results.Yfit(:,j:k) = tmpYfit;
  results.S(:,:,i) = tmpS;
end;

% calculate residual error
[results.E,results.norm_var_E,...
 results.Yeeg,results.Yfiteeg,results.Eeeg,...
 results.var_Eeeg,results.var_Yeeg,results.norm_var_Eeeg,...
 results.Ygrad,results.Yfitgrad,results.Egrad,results.var_Egrad,...
 results.var_Ygrad,results.norm_var_Egrad,...
 results.Ymag,results.Yfitmag,results.Emag,...
 results.var_Emag,results.var_Ymag,results.norm_var_Emag]...
  = rc_RCSE_calc_fit_err(results.Y,results.Yfit,parms);

% calculate overall error / cost
if isempty(parms.prior) || parms.prior_weight<1
  if parms.fit_range_flag
    % average normalized error over specified time range
    parms.total_error = mean(results.norm_var_E(parms.fit_t0:parms.fit_t1));
  else
    % average normalized error across entire time range
    parms.total_error = mean(results.norm_var_E);
  end;
  % increase cost if polarity in specified time range is positive (or negative)
  if parms.polarity_penalty
    tmp = results.S(parms.fit_t0:parms.fit_t1,:);
    tmp = reshape(tmp,[size(tmp,1),size(tmp,2)*size(tmp,3)]);
    penalty = 0;
    for s=1:size(tmp,2)
      tmp2 = mean(squeeze(tmp(:,s)));
      if (tmp2>0 && parms.polarity_penalty>0) ||...
         (tmp2<0 && parms.polarity_penalty<0)
        penalty = penalty + tmp2*parms.polarity_penalty;
      end;    
    end;
    parms.total_error = parms.total_error*(1+penalty);
  end;
end;

if ~isempty(parms.prior)
  if parms.prior_zscore_flag
    S_sem_ref = parms.prior.S_sem_ref;
  else
    S_sem_ref = [];
  end;
  if parms.prior_mindiff_flag==1
    % absolute difference between prior and current estimate
    prior_error = ...
      rc_calc_ref_absdiff(results.S,parms.prior.S_ref,...
        'S_sem_ref',S_sem_ref,...
        'sem_min',parms.prior_sem_min,...
        'srange',[parms.corr_t0,parms.corr_t1]);
  elseif parms.prior_mindiff_flag==2
    % optimally scale S_ref then calculate absolute difference
    prior_error = ...
      rc_calc_ref_absdiff(results.S,parms.prior.S_ref,...
        'S_sem_ref',S_sem_ref,...
        'sem_min',parms.prior_sem_min,...
        'srange',[parms.corr_t0,parms.corr_t1],...
        'scale_flag',1,...
        'indy_wform_flag',parms.prior_indy_wform_flag);
  else
    % correlation with prior as cost function
    prior_error = 1.0 -...
      rc_calc_ref_corr(results.S,parms.parms.prior.S_ref,[parms.corr_t0,parms.corr_t1]);
  end;
  parms.total_error = parms.prior_weight * prior_error +...
                  (1-parms.prior_weight) * parms.total_error;
end;

results.retmap = retforward.retmap;

results.areas = results.retmap.use_areas;
results.nareas = results.retmap.num_areas;
results.contrasts = [results.retmap.orig_cond_info.contrast];
results.ncontrasts = length(unique(results.contrasts(results.contrasts>0)));
