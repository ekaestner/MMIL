function [parms,results] = rc_RCSE_rand_offset_search(parms,avg_data,forward)
%function [parms,results] = rc_RCSE_rand_offset_search(parms,avg_data,forward)
%
% Purpose: optimize dipole offset through random search
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
for iter=1:parms.tmp_offset_niters
  if iter>1
    parms.total_offset_niters = parms.total_offset_niters + 1;
    % select random offset steps for each area and condition
    r_steps = parms.r_step*randn(parms.nareas,parms.nconds,parms.npatches);
    th_steps = parms.th_step*randn(parms.nareas,parms.nconds,parms.npatches);
    % force some conditions to have same steps (e.g. within quad)
    for i=1:parms.nconds
      j=parms.offset_group_conds(i);
      if i==j, continue; end;
      r_steps(:,i,:) = r_steps(:,j,:);
      th_steps(:,i,:) = th_steps(:,j,:);
    end;
    if parms.offset_group_patches_flag % force patches to have same offsets
      for p=2:parms.npatches
        r_steps(:,:,p) = r_steps(:,:,1);
        th_steps(:,:,p) = th_steps(:,:,1);
      end;
    end;
    if parms.offset_group_areas_flag % force areas to have same offsets
      for a=2:parms.nareas
        r_steps(a,:,:) = r_steps(1,:,:);
        th_steps(a,:,:) = th_steps(1,:,:);
      end;
    end;
    for k=1:length(parms.offset_const_areas)
      % keep offsets constant (set to zero)
      a = parms.offset_const_areas(k);
      r_steps(a,:,:) = 0;
      th_steps(a,:,:) = 0;
    end;
  else
    r_steps = zeros(parms.nareas,parms.nconds,parms.npatches);
    th_steps = zeros(parms.nareas,parms.nconds,parms.npatches);
  end;
  tmp_r_offsets = parms.tmp_best_r_offsets + r_steps;
  tmp_th_offsets = parms.tmp_best_th_offsets + th_steps;
  tmp_cond_info = rc_RCSE_set_cond_info_offsets(parms,tmp_r_offsets,tmp_th_offsets);
  retmap = rc_RCSE_retmap_from_retfit(parms,0,0,tmp_cond_info,forward);
  [parms,results,retforward,inverse]=...
    rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
  if isempty(best_results), best_results = results; end;
  fprintf('%s: rand offset search iteration %d: error = %f\n',...
    mfilename,iter,parms.total_error);
  % if error is reduced, use current dipoles as benchmark
  if ~isnan(parms.total_error) & ...
     ((parms.total_error<parms.tmp_min_error) |...
     isempty(parms.tmp_best_retmap))
    best_results = results;
    parms.tmp_min_error = parms.total_error;
    parms.tmp_best_retmap = retforward.retmap;
    parms.tmp_best_r_offsets = tmp_r_offsets;
    parms.tmp_best_th_offsets = tmp_th_offsets;
    parms.tmp_best_cond_info = tmp_cond_info;
    num_better_fits = num_better_fits + 1;
    fprintf('%s: ### rand offset search better fit\n',mfilename);
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

fprintf('%s: ### number of rand offset search better fits = %d\n',mfilename,num_better_fits);
