function [parms,results] = rc_RCSE_hybrid_offset_search(parms,avg_data,forward)
%function [parms,results] = rc_RCSE_hybrid_offset_search(parms,avg_data,forward)
%
% Purpose: calculate source estimates excluding locations with larger error
%   perform multi-start random offset search to optimize dipole offsets
%   for good locations; then optimize dipoles for bad locations, using estimates
%   from good locations as prior
%
% Required Input:
%   parms: RCSE parms struct
%   avg_data: average data struct
%   forward: forward solution struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 05/18/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

orig_parms = parms;

% initialize retmap
retmap = ...
  rc_RCSE_retmap_from_retfit(parms,parms.r_offset,parms.th_offset,[],forward);

% calculate source estimates
[f0_parms,f0_results]=...
  rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
f0_parms.best_retmap = f0_results.retmap;

% calculate average abs(error) for each condition
num_sensors = length(parms.goodchans);
if parms.fit_range_flag
  fit_range = [parms.fit_t0:parms.fit_t1];
else
  fit_range = [1:size(f0_results.E,1)];
end;
errvec = [];
for i=1:length(parms.cond_info)
  % calculate average error for this condition
  j = 1+(i-1)*num_sensors;
  k = j+num_sensors-1;
  err = f0_results.E(fit_range,j:k);
  err = mean(abs(err(:)));
  errvec = [errvec err];
end

% identify conditions with best and worst err
nbins = 100;
[N,X] = hist(errvec,nbins);
P = cumsum(N); Pnorm = P/P(end);
[tmp,ind_min] = min(abs(Pnorm-parms.hybrid_cutoff));
errmin = X(ind_min);
ind_good = find(errvec>errmin);
ind_bad = find(errvec<=errmin);

% for stim locs with best correlation recalculate fRCSE
%   with random offset search if offset_niters>0
parms = orig_parms;
parms.event_codes = parms.event_codes(ind_good);
parms.conditions = parms.conditions(ind_good);
parms = rc_RCSE_check_conds(parms,avg_data);
parms.offset_niters = parms.hybrid_pre_niters;
if parms.offset_niters>0  % multi-start rand offset search
  fprintf('%s: searching for best offsets with good stim locs...\n',mfilename);
  [f1_parms,f1_results] = ...
    rc_RCSE_mstart_offset_search(parms,avg_data,forward);    
else
  if parms.num_r_offsets>1 || parms.num_th_offsets>1
    fprintf('%s: searching for best offsets with good stim locs...\n',mfilename);
    % if r_offset or th_offset are vectors, do exhaustive search
    f1_parms = rc_RCSE_full_offset_search(parms,avg_data,forward);
    parms.retmap = f1_parms.best_retmap;
  else
    fprintf('%s: recalculating source estimates with good stim locs...\n',mfilename);
    parms.retmap = ...
      rc_RCSE_retmap_from_retfit(parms,parms.r_offset,parms.th_offset,[],forward);
  end;
  [f1_parms,f1_results]=...
    rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
  if ~isfield(f1_parms,'best_retmap') | isempty(f1_parms.best_retmap)
    f1_parms.best_retmap = f1_results.retmap;
  end;
end;
%  f1_wforms = rc_RCSE_extract_source_wforms(f1_parms,avg_data,f1_results);

% for stim locs with worst correlation, search for best offsets
%   to give max corr with new fRCSE (one area x stim loc at a time?)
%% todo: include all stim locs?
%% todo: one stim loc at a time?
parms = orig_parms;
parms.event_codes = parms.event_codes(ind_bad);
parms.conditions = parms.conditions(ind_bad);
parms = rc_RCSE_check_conds(parms,avg_data);
parms.prior = rc_extract_prior(f1_results);
if parms.offset_niters>0  % multi-start rand offset search
  fprintf('%s: searching for best offsets with bad stim locs...\n',mfilename);
  [f2_parms,f2_results] = ...
    rc_RCSE_mstart_offset_search(parms,avg_data,forward);    
else
  if parms.num_r_offsets>1 || parms.num_th_offsets>1
    fprintf('%s: searching for best offsets with bad stim locs...\n',mfilename);
    % if r_offset or th_offset are vectors, do exhaustive search
    f2_parms = rc_RCSE_full_offset_search(parms,avg_data,forward);
    parms.retmap = f2_parms.best_retmap;
  else
    fprintf('%s: recalculating source estimates with bad stim locs...\n',mfilename);
    parms.retmap = ...
      rc_RCSE_retmap_from_retfit(parms,parms.r_offset,parms.th_offset,[],forward);
  end;
  [f2_parms,f2_results]=...
    rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
  if ~isfield(f2_parms,'best_retmap') | isempty(f2_parms.best_retmap)
    f2_parms.best_retmap = f2_parms.retmap;
  end;
end;
%  f2_wforms = rc_RCSE_extract_source_wforms(f2_parms,avg_data,f2_results);

% combine offsets for all conditions
cond_info = [];
fields = fieldnames(f1_parms.best_cond_info);
for i=1:length(fields)
  cond_info.(fields{i}) = [];
end;
cond_info(ind_good) = f1_parms.best_cond_info;
cond_info(ind_bad) = f2_parms.best_cond_info;

% recalculate results
parms = orig_parms;
retmap = rc_RCSE_retmap_from_retfit(parms,0,0,cond_info,forward);
[parms,results]=...
  rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
parms.best_retmap = results.retmap;

