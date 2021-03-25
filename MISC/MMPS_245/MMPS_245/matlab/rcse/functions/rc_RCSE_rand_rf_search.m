function [parms,results] = rc_RCSE_rand_rf_search(parms,avg_data,forward)
%function [parms,results] = rc_RCSE_rand_rf_search(parms,avg_data,forward)
%
% Purpose: optimize receptive field sizes through random search
%
% Required Input:
%   parms: RCSE parms struct
%   avg_data: average data struct
%   forward: forward solution struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 12/08/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

num_better_fits = 0;
best_results = [];
if ~isfield(parms,'tmp_best_rf_sizes')
  parms.tmp_best_rf_sizes = parms.best_rf_sizes;
end;
if isempty(parms.tmp_best_rf_sizes)
  parms.tmp_best_rf_sizes = parms.init_rf_sizes;
end;
if ~isfield(parms,'tmp_best_rf_slopes')
  parms.tmp_best_rf_slopes = parms.best_rf_slopes;
end;
if isempty(parms.tmp_best_rf_slopes)
  parms.tmp_best_rf_slopes = parms.init_rf_slopes;
end;
if ~isfield(parms,'tmp_min_error')
  parms.tmp_min_error = parms.min_error;
end;

% choose offsets for retfit to retmap
if ~isfield(parms,'tmp_best_r_offsets')
  tmp_r_offsets = zeros(parms.nareas,parms.nconds,parms.npatches);
else
  tmp_r_offsets = parms.tmp_best_r_offsets;
end;    
if ~isfield(parms,'tmp_best_th_offsets')
  tmp_th_offsets = zeros(parms.nareas,parms.nconds,parms.npatches);
else
  tmp_th_offsets = parms.tmp_best_th_offsets;
end;
tmp_cond_info = rc_RCSE_set_cond_info_offsets(parms,tmp_r_offsets,tmp_th_offsets);

% find best rf sizes
for iter=1:parms.rf_niters
  if iter>1
    parms.total_rf_niters = parms.total_rf_niters + 1;
    % randomly adjust field sizes
    tmp_rf_sizes = parms.tmp_best_rf_sizes +...
                   parms.rf_size_step*randn(size(parms.tmp_best_rf_sizes));
    % randomly adjust rf ecc slope
    tmp_rf_slopes = parms.tmp_best_rf_slopes +...
                   parms.rf_slope_step*randn(size(parms.tmp_best_rf_slopes));
  else
    % use previous best estimate
    tmp_rf_sizes = parms.tmp_best_rf_sizes;
    tmp_rf_slopes = parms.tmp_best_rf_slopes;
  end;
  % check bounds
  tmp_min_size = parms.rf_min_sizes(1);
  tmp_max_size = parms.rf_max_sizes(1);
  tmp_min_slope = parms.rf_min_slopes(1);
  tmp_max_slope = parms.rf_max_slopes(1);
  for i=1:length(parms.tmp_best_rf_sizes)
    % assume that higher areas have larger receptive fields
    tmp_rf_sizes(i) = rc_check_bounds(tmp_rf_sizes(i),...
      [tmp_min_size,tmp_max_size]);
    if i<length(parms.tmp_best_rf_sizes)
      tmp_min_size = max([parms.rf_min_sizes(i+1),tmp_rf_sizes(i)]);
      tmp_max_size = parms.rf_max_sizes(i+1);
    end;
    % assume that higher areas have larger slope
    tmp_rf_slopes(i) = rc_check_bounds(tmp_rf_slopes(i),...
      [tmp_min_slope,tmp_max_slope]);
    if i<length(parms.tmp_best_rf_slopes)
      tmp_min_slope = max([parms.rf_min_slopes(i+1),tmp_rf_slopes(i)]);
      tmp_max_slope = parms.rf_max_slopes(i+1);
    end;
  end;
  parms.rf_sizes = tmp_rf_sizes;
  parms.rf_slopes = tmp_rf_slopes;
  % need to re-compute vf2ctx
  if ~isempty(parms.retfit_results) && parms.vf2ctx_flag
    if parms.verbose
      fprintf('%s: calculating mapping between visual field and cortex...\n',...
        mfilename);
      tic
    end;
    args = mmil_parms2args(parms,parms.vf2ctx_tags);
    [parms.vf2ctx,parms.retfit_results] = ...
      rc_calc_vf2ctx(parms.retfit_results,args{:});
    if parms.verbose,toc; end;
  end;
  % create retmap
  retmap = rc_RCSE_retmap_from_retfit(parms,0,0,tmp_cond_info,forward);
  % calculate goodness of fit
  [parms,results,retforward,inverse]=...
    rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
  if isempty(best_results), best_results = results; end;
  fprintf('%s: rand rf sizes search iteration %d: error = %f\n',...
    mfilename,iter,parms.total_error);
  fprintf('     rf sizes: [%s], rf slopes: [%s]\n',...
    sprintf('%0.2f ',tmp_rf_sizes),...
    sprintf('%0.2f ',tmp_rf_slopes));
  % if error is reduced, use current dipoles as benchmark
  if ~isnan(parms.total_error) && ...
     ((parms.total_error<parms.tmp_min_error) ||...
     isempty(parms.tmp_best_retmap))
    best_results = results;
    parms.tmp_min_error = parms.total_error;
    parms.tmp_best_retmap = retforward.retmap;
    parms.tmp_best_rf_sizes = tmp_rf_sizes;
    parms.tmp_best_rf_slopes = tmp_rf_slopes;
    parms.tmp_best_vf2ctx = parms.vf2ctx;
    num_better_fits = num_better_fits + 1;
    fprintf('%s: ### rand rf sizes search better fit\n',mfilename);
    fprintf('     best rf sizes: [%s], rf slopes: [%s]\n',...
      sprintf('%0.2f ',parms.tmp_best_rf_sizes),...
      sprintf('%0.2f ',parms.tmp_best_rf_slopes));
  end;
  if parms.plotflag
    fprintf('%s: plotting...\n',mfilename);
    figure(1);
    rc_plot_sources('parms',parms,'results',results,...
      'ylim',parms.plot_ylim,'units',parms.source_units,...
      'area_names',parms.area_names,'area_colors',parms.area_colors,...
      'labelflag',0);
    figure(2);
    rc_plot_fitvar_results(parms,results,0);
    drawnow;
  end;
end;
results = best_results;
fprintf('%s: ### number of rand rf sizes search better fits = %d\n',mfilename,num_better_fits);
fprintf('     best rf sizes: [%s], rf slopes: [%s]\n',...
  sprintf('%0.2f ',parms.tmp_best_rf_sizes),...
  sprintf('%0.2f ',parms.tmp_best_rf_slopes));

parms.rf_sizes = parms.tmp_best_rf_sizes;
parms.rf_slopes = parms.tmp_best_rf_slopes;
if ~(parms.mstart_flag && parms.mstart_rfsearch_flag)
  parms.best_rf_sizes = parms.tmp_best_rf_sizes;
  parms.best_rf_slopes = parms.tmp_best_rf_slopes;
  parms.vf2ctx = parms.tmp_best_vf2ctx;
  parms.best_retmap = parms.tmp_best_retmap;
  parms.min_error = parms.tmp_min_error;
end;
