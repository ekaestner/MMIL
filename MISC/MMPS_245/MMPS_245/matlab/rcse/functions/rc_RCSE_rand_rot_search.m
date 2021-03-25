function [parms,results] = rc_RCSE_rand_rot_search(parms,avg_data,forward)
%function [parms,results] = rc_RCSE_rand_rot_search(parms,avg_data,forward)
%
% Purpose: optimize coordinate system rotation through random search
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

num_better_fits = 0;
best_results = [];
parms.tmp_best_rot_forward = [];
if ~isfield(parms,'tmp_best_rot')
  parms.tmp_best_rot = parms.best_rot;
end;
if ~isfield(parms,'tmp_min_error')
  parms.tmp_min_error = parms.min_error;
end;

% create retmap from retfit
if ~isempty(parms.retfit_results)
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
  retmap = rc_RCSE_retmap_from_retfit(parms,0,0,tmp_cond_info,forward);
end;

for iter=1:parms.rot_niters
  if iter>1
    parms.total_rot_niters = parms.total_rot_niters + 1;
    % select random rotation
    tmp_rot = parms.tmp_best_rot +...
              parms.rot_step*randn(size(parms.tmp_best_rot));
  else
    % use previous best estimate (0 if first time)
    tmp_rot = parms.tmp_best_rot;
  end;
  tmp_rot = rc_check_bounds(tmp_rot,[-parms.rot_max,parms.rot_max]);
  tmp_rot = rc_check_wrap(tmp_rot,[-180,180]);
  % rotate dipoles
  tmp_trans = rc_rotate_trans(parms.trans,tmp_rot);
  tmp_forward = forward;
  [tmp_forward.G_norm,tmp_forward.G_tang1,tmp_forward.G_tang2]=...
    ts_gain_xyz2norm(tmp_forward.G_xyz,...
    forward.lh_dip_info,forward.rh_dip_info,...
    forward.lh_dec_dips,forward.rh_dec_dips,tmp_trans);

  % fit waveforms with adjusted forward solution
  [parms,results,retforward,inverse]=...
    rc_RCSE_fit_waveforms(parms,retmap,avg_data,tmp_forward);
  if isempty(best_results), best_results = results; end;
  fprintf('%s: rotation search iteration %d: error = %f\n',...
    mfilename,iter,parms.total_error);
  fprintf('     rotation: [%s] degrees\n',sprintf('%0.2f ',tmp_rot));
  % if error is reduced, use current dipoles as benchmark
  if ~isnan(parms.total_error) && ...
     ((parms.total_error<parms.tmp_min_error) ||...
     isempty(parms.tmp_best_retmap) ||...
     isempty(parms.tmp_best_rot_forward))
    best_results = results;
    parms.tmp_min_error = parms.total_error;
    parms.tmp_best_rot = tmp_rot;
    parms.tmp_best_rot_forward = tmp_forward;
    num_better_fits = num_better_fits + 1;
    fprintf('%s: ### rotation search better fit\n',mfilename);
    fprintf('     best rotation: [%s] degrees\n',sprintf('%0.2f ',parms.tmp_best_rot));
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
fprintf('%s: ### number of rotation search better fits = %d\n',mfilename,num_better_fits);
fprintf('     best rotation: [%s] degrees\n',sprintf('%0.2f ',parms.tmp_best_rot));

if ~(parms.mstart_flag && parms.mstart_rotsearch_flag)
  parms.best_rot_forward = parms.tmp_best_rot_forward;
  parms.best_rot = parms.tmp_best_rot;
  parms.min_error = parms.tmp_min_error;
end;
