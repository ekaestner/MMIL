function parms = rc_RCSE_full_nbrhood_search(parms,avg_data,forward)
%function parms = rc_RCSE_full_nbrhood_search(parms,avg_data,forward)
%
% Purpose: optimize dipoles by choosing single vertex within patch
%   for each stimulus - exhaustive search
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
parms.best_retmap = retmap;
for iter=1:parms.nbrhd_niters
  parms.total_nbrhd_niters = parms.total_nbrhd_niters + 1;
  % loop over areas and stimulus locations
  for a=parms.use_areas  
    for c=1:parms.nconds
      % exhaustive search within patch

      % left hemisphere dipoles for this location
      v_lh = retmap.orig_areas(a).verts(c).v_lh;
      w_lh = retmap.orig_areas(a).verts(c).w_lh;
      % right hemisphere dipoles for this location
      v_rh = retmap.orig_areas(a).verts(c).v_rh;
      w_rh = retmap.orig_areas(a).verts(c).w_rh;

      % threshold relative to max value
      if parms.w_thresh_patch>0
        % normalize to max across both hemispheres
        w = rc_norm_weights([w_lh;w_rh],1);
        w_lh = w(1:length(w_lh));
        w_rh = w(length(w_lh)+1:end);
        % apply threshold
        [w_lh,v_lh]=rc_thresh_weights(w_lh,v_lh,parms.w_thresh_patch);
        [w_rh,v_rh]=rc_thresh_weights(w_rh,v_rh,parms.w_thresh_patch);
      end;

      n_lh = length(v_lh);
      n_rh = length(v_rh);
      fprintf('%s: searching through %d vertices for cond %d...\n',...
        mfilename,n_lh+n_rh,c);
      parms.tmp_min_error = 10^29;
      for i=1:n_lh
        tmp_retmap = parms.best_retmap; % current best estimate
        % choose from patch
        tmp_retmap.orig_areas(a).verts(c).v_lh = v_lh(i);
        tmp_retmap.orig_areas(a).verts(c).w_lh = 1;
        [parms,results,retforward,inverse]=...
          rc_RCSE_fit_waveforms(parms,tmp_retmap,avg_data,forward);
        % if error is reduced, use current dipoles as benchmark
        fprintf('%s: nbrhd search iteration %d, area %d, condition %d, lh vertex %d: error = %f\n',...
          mfilename,iter,a,c,i,parms.total_error);
        if ~isnan(parms.total_error) & ...
           ((parms.total_error<parms.tmp_min_error) )
          parms.tmp_min_error = parms.total_error;
          parms.best_retmap = retforward.retmap;
          num_better_fits = num_better_fits + 1;
          fprintf('%s: ### nbrhd search better fit\n',mfilename);
        end;
      end;
      parms.tmp_min_error = 10^29;
      for i=1:n_rh
        tmp_retmap = parms.best_retmap; % current best estimate
        tmp_retmap.orig_areas(a).verts(c).v_rh = v_rh(i);
        tmp_retmap.orig_areas(a).verts(c).w_rh = 1;
        [parms,results,retforward,inverse]=...
          rc_RCSE_fit_waveforms(parms,tmp_retmap,avg_data,forward);
        % if error is reduced, use current dipoles as benchmark
        fprintf('%s: nbrhd search iteration %d, area %d, condition %d, rh vertex %d: error = %f\n',...
          mfilename,iter,a,c,i,parms.total_error);
        if ~isnan(parms.total_error) & ...
           ((parms.total_error<parms.tmp_min_error) )
          parms.tmp_min_error = parms.total_error;
          parms.best_retmap = retforward.retmap;
          num_better_fits = num_better_fits + 1;
          fprintf('%s: ### nbrhd search better fit\n',mfilename);
        end;
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
  end;
  parms.simple_inverse_flag = 0;
  [parms,results,retforward,inverse]=...
    rc_RCSE_fit_waveforms(parms,parms.best_retmap,avg_data,forward);
  parms.min_error = parms.total_error;
  parms.best_retmap = retforward.retmap;
  fprintf('%s: nbrhd search iteration %d: error = %f\n',...
    mfilename,iter,parms.total_error);
  if parms.plotflag
    fprintf('%s: plotting...\n',mfilename);
    figure(1);
    plot_sources('parms',parms,'results',results,...
      'ylim',parms.plot_ylim,'units',parms.source_units,...
      'labelflag',0);
    figure(2);
    plot_fitvar_results(parms,results,0);
    drawnow;
  end;
end;
fprintf('%s: ### number of nbrhd search better fits = %d\n',mfilename,num_better_fits);
fprintf('%s: ### number of total nbrhd iters = %d\n',...
  mfilename,parms.total_nbrhd_niters);
fprintf('%s: ### nbrhd search min error = %0.4f\n',mfilename,parms.min_error);
