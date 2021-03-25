function [parms,results] = rc_RCSE_mstart_offset_search(parms,avg_data,forward)
%function [parms,results] = rc_RCSE_mstart_offset_search(parms,avg_data,forward)
%
% Purpose: multi-start dipole offset optimization
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

parms.best_r_offsets = parms.init_r_offsets;
parms.best_th_offsets = parms.init_th_offsets;
parms.best_cond_info = parms.cond_info;
best_results = [];
num_better_mstart_fits = 0;
if parms.mstart_average_flag
  sourcewf_sum = single(0);
  sourcewf_sumsq = single(0);
end;
tmp_forward = forward; % for rotation search
for ms_iter=1:parms.offset_mstarts
  parms.ms_iter = ms_iter;
  parms = rc_RCSE_init_mstart(parms);

  % search for best offset (to account for errors in map fit)
  if parms.fmincon_flag
    [parms,results] = rc_RCSE_fmincon_offset_search(parms,avg_data,tmp_forward);
  else
    [parms,results] = rc_RCSE_rand_offset_search(parms,avg_data,tmp_forward);
  end;
  
  % search for best rf sizes
  if parms.rf_niters && parms.offset_mstarts>1 &&...
      parms.mstart_rfsearch_flag &&...
      ~(ms_iter==parms.offset_mstarts && parms.offset_niters_last>0)
    [parms,results] = rc_RCSE_rand_rf_search(parms,avg_data,tmp_forward);
  end;

  % search for best rotation (to account for errors in registration)
  if parms.rot_niters && parms.offset_mstarts>1 &&...
      parms.mstart_rotsearch_flag &&...
      ~(ms_iter==parms.offset_mstarts && parms.offset_niters_last>0)
    [parms,results] = rc_RCSE_rand_rot_search(parms,avg_data,tmp_forward);
    tmp_forward = parms.tmp_best_rot_forward;
  end;

  if isempty(best_results), best_results = results; end;
  if parms.mstart_flag
    fprintf('%s: multistart iteration %d: error = %f\n',...
      mfilename,parms.ms_iter,parms.tmp_min_error);
  end;
  % if error is reduced, use current dipoles as benchmark
  if parms.tmp_min_error<parms.min_error | isempty(parms.best_retmap)
    best_results = results;
    parms.min_error = parms.tmp_min_error;
    parms.best_retmap = parms.tmp_best_retmap;
    parms.best_r_offsets = parms.tmp_best_r_offsets;
    parms.best_th_offsets = parms.tmp_best_th_offsets;
    parms.best_cond_info = parms.tmp_best_cond_info;
    if parms.rf_niters && ...
       parms.mstart_flag && parms.mstart_rfsearch_flag
      parms.best_rf_sizes = parms.tmp_best_rf_sizes;
      parms.best_rf_slopes = parms.tmp_best_rf_slopes;
      parms.rf_sizes = parms.tmp_best_rf_sizes;
      parms.rf_slopes = parms.tmp_best_rf_slopes;
      parms.vf2ctx = parms.tmp_best_vf2ctx;
    end;
    if parms.rot_niters && ...
       parms.mstart_flag && parms.mstart_rotsearch_flag
      parms.best_rot = parms.tmp_best_rot;
      parms.best_rot_forward = parms.tmp_best_rot_forward;
    end;
    if parms.offset_mstarts>1
      num_better_mstart_fits = num_better_mstart_fits + 1;
      fprintf('%s: ### better multistart fit\n',mfilename);
    end;
  end;

  % add current waveforms to running total for later averaging
  if parms.mstart_average_flag
    sourcewf_sum = sourcewf_sum + results.S;
    sourcewf_sumsq = sourcewf_sumsq + results.S.^2;
  end;
end;
if parms.mstart_average_flag
  N = parms.offset_mstarts;
  if N<2
    fprintf('%s: WARNING: less than two iterations in mstart average\n',...
      mfilename);
  else
    fprintf('%s: %d iterations included in mstart average\n',mfilename,N);
    avgS = sourcewf_sum / (eps+N);
    stdS = sqrt((N*sourcewf_sumsq - sourcewf_sum.^2)./(eps+N*(N-1)));
    best_results.S = avgS;
    best_results.S_stdv = stdS;
  end;
end;
results = best_results;

for a=1:parms.nareas
  if ~ismember(a,parms.use_areas), continue; end;
  fprintf('%s: %s best r_offsets:\n  ',mfilename,parms.area_names{a});
  fprintf('%0.4f ',squeeze(parms.best_r_offsets(a,:,:)));
  fprintf('\n');
  fprintf('%s: %s best th_offsets:\n  ',mfilename,parms.area_names{a});
  fprintf('%0.4f ',squeeze(parms.best_th_offsets(a,:,:)));
  fprintf('\n');
end;
if parms.rf_niters && ...
   parms.offset_mstarts>1 && parms.mstart_rfsearch_flag
  fprintf('%s: best rf sizes: [%s], rf slopes = [%s]\n',mfilename,...
    sprintf('%0.2f ',parms.best_rf_sizes),...
    sprintf('%0.2f ',parms.best_rf_slopes));
  fprintf('%s: ### number of total rf sizes iters = %d\n',...
    mfilename,parms.total_rf_niters);
end;
if parms.rot_niters && ...
   parms.offset_mstarts>1 && parms.mstart_rotsearch_flag
  fprintf('%s: best rotation: [%s] degrees\n',...
    mfilename,sprintf('%0.2f ',parms.best_rot));
  fprintf('%s: ### number of total rot iters = %d\n',...
    mfilename,parms.total_rot_niters);
end;
fprintf('%s: ### number of total offset iters = %d\n',...
  mfilename,parms.total_offset_niters);
if parms.total_offset_mstarts>0
  fprintf('%s: ### number of better multi-start fits = %d\n',...
    mfilename,num_better_mstart_fits);
  fprintf('%s: ### number of total offset multi-starts = %d\n',...
    mfilename,parms.total_offset_mstarts);
end;
fprintf('%s: ### min error = %0.4f\n',mfilename,parms.min_error);

