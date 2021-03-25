function parms = rc_RCSE_rand_nbrhood_search(parms,avg_data,forward)
%function parms = rc_RCSE_rand_nbrhood_search(parms,avg_data,forward)
%
% Purpose: optimize dipoles by choosing single vertex within patch
%   for each stimulus - random search
%
% Required Input:
%   parms: RCSE parms struct
%   avg_data: average data struct
%   forward: forward solution struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 10/24/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

if ~isfield(parms.retmap,'M')
  % construct ret mapping structure for all areas
  use_areas = [1:length(parms.retmap.areas)];
  % use unique_location_conds only (e.g. multiple contrast levels)
  tmp_event_codes = []; tmp_conditions = [];
  if ~isempty(parms.event_codes)
    tmp_event_codes = parms.event_codes(parms.unique_location_conds);
  end;
  if ~isempty(parms.conditions)
    tmp_conditions = parms.conditions(parms.unique_location_conds);
  end;
  retmap = rc_construct_ret_mapping(parms.retmap,...
    'use_areas',parms.use_areas,'conditions',tmp_conditions,...
    'event_codes',tmp_event_codes,parms.ret_mapping_opts{:});
  if isempty(retmap)
    error('failed to create retinotopy mapping matrix');
  end
else
  retmap = parms.retmap;
end;
num_better_fits = 0;
for iter=1:parms.nbrhd_niters
  parms.total_nbrhd_niters = parms.total_nbrhd_niters + 1;
  % randomly choose dipoles within neighborhoods for each stim loc and visual area
  tmp_retmap = rc_choose_retdips(retmap);
  [parms,results,retforward,inverse]=...
    rc_RCSE_fit_waveforms(parms,tmp_retmap,avg_data,forward);
  fprintf('%s: nbrhd search iteration %d: error = %f\n',mfilename,iter,parms.total_error);
  % if error is reduced, use current dipoles as benchmark
  if ~isnan(parms.total_error) & ...
     ((parms.total_error<parms.min_error) | isempty(parms.best_retmap))
    parms.min_error = parms.total_error;
    parms.best_retmap = retforward.retmap;
    num_better_fits = num_better_fits + 1;
    fprintf('%s: ### nbrhd search better fit\n',mfilename);
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
fprintf('%s: ### number of nbrhd search better fits = %d\n',mfilename,num_better_fits);
fprintf('%s: ### number of total nbrhd iters = %d\n',...
  mfilename,parms.total_nbrhd_niters);
fprintf('%s: ### nbrhd search min error = %0.4f\n',mfilename,parms.min_error);
