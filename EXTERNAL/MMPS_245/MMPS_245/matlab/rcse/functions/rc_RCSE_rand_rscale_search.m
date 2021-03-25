function [parms,results] = rc_RCSE_rand_rscale_search(parms,avg_data,forward)
%function [parms,results] = rc_RCSE_rand_rscale_search(parms,avg_data,forward)
%
% Purpose: optimize eccentricity scaling through random search
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
if ~isfield(parms,'tmp_best_rscale')
  parms.tmp_best_rscale = parms.best_rscale;
end;
if ~isfield(parms,'tmp_min_error')
  parms.tmp_min_error = parms.min_error;
end;
parms.orig_r_max = parms.r_max;

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

% find best rscale
for iter=1:parms.rscale_niters
  if iter>1
    parms.total_rscale_niters = parms.total_rscale_niters + 1;
    % randomly adjust rscale
    tmp_rscale = parms.tmp_best_rscale + parms.rscale_step*randn;
  else
    % use previous best estimate
    tmp_rscale = parms.tmp_best_rscale;
  end;
  % check bounds
  tmp_rscale = rc_check_bounds(tmp_rscale,parms.rscale_range);
  parms.r_max = parms.orig_r_max * tmp_rscale;
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
  fprintf('%s: rand rscale search iteration %d: error = %f\n',...
    mfilename,iter,parms.total_error);
  fprintf('     rscale: %0.3f, r_max:  %0.3f\n',tmp_rscale,parms.r_max);
  % if error is reduced, use current dipoles as benchmark
  if ~isnan(parms.total_error) && ...
     ((parms.total_error<parms.tmp_min_error) ||...
     isempty(parms.tmp_best_retmap))
    best_results = results;
    parms.tmp_min_error = parms.total_error;
    parms.tmp_best_retmap = retforward.retmap;
    parms.tmp_best_rscale = tmp_rscale;
    num_better_fits = num_better_fits + 1;
    fprintf('%s: ### rand rscale search better fit\n',mfilename);
    fprintf('     best rscale: %0.3f, best r_max:  %0.3f\n',...
      tmp_rscale,parms.r_max);
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
fprintf('%s: ### number of rand rscale search better fits = %d\n',mfilename,num_better_fits);
fprintf('     best rscale: %0.3f, best r_max:  %0.3f\n',...
  tmp_rscale,parms.r_max);
parms.best_rscale = parms.tmp_best_rscale;
parms.r_max = parms.orig_r_max * parms.best_rscale;
parms.min_error = parms.tmp_min_error;

