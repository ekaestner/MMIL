function parms = rc_RCSE_init_mstart(parms)
%function parms = rc_RCSE_init_mstart(parms)
%
% Purpose: initialize parameters for multi-start dipole offset optimization
%
% Required Input:
%   parms: RCSE parms struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/04/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

% initialize tmp_best
fprintf('%s: multistart initialization...\n',mfilename);
parms.tmp_offset_niters = parms.offset_niters;
if (parms.ms_iter==1 && parms.total_offset_mstarts==0) ||...
  ~parms.mstart_randstart_flag || ...
  (parms.ms_iter==parms.offset_mstarts && parms.offset_niters_last>0)
  fprintf('%s: using current best offsets...\n',mfilename);
  % use current best offsets
  parms.tmp_min_error = parms.min_error;
  parms.tmp_best_r_offsets = parms.best_r_offsets;
  parms.tmp_best_th_offsets = parms.best_th_offsets;
  parms.tmp_best_retmap = parms.best_retmap;
  % use different number of iterations if this is the last mstart
  %% todo: could also use a smaller r_step, th_step for this search
  %%       or use offset_group_flag = 0
  if parms.mstart_flag &&...
     parms.offset_niters_last>0 && parms.ms_iter==parms.offset_mstarts
    parms.tmp_offset_niters = parms.offset_niters_last;
  end;
  parms.tmp_best_rot = parms.best_rot;
  parms.tmp_best_rf_sizes = parms.best_rf_sizes;
  parms.tmp_best_rf_slopes = parms.best_rf_slopes;
else
  fprintf('%s: randomly selecting offsets...\n',mfilename);
  parms.tmp_min_error = 10^29;
  % select random offsets for each area and condition
  parms.tmp_best_r_offsets = parms.r_offset_range(1) +...
    (parms.r_offset_range(2)-parms.r_offset_range(1))*...
    rand(parms.nareas,parms.nconds,parms.npatches);
  parms.tmp_best_th_offsets = parms.th_offset_range(1) +...
    (parms.th_offset_range(2)-parms.th_offset_range(1))*...
    rand(parms.nareas,parms.nconds,parms.npatches);
  parms.tmp_best_retmap = [];
  parms.tmp_best_rot = [0 0 0];
  parms.tmp_best_rf_sizes = parms.init_rf_sizes;
  parms.tmp_best_rf_slopes = parms.init_rf_slopes;
end;
% force some conditions to have same offsets (e.g. within quad)
for i=1:parms.nconds
  j=parms.offset_group_conds(i);
  if i==j, continue; end;
  parms.tmp_best_r_offsets(:,i,:) = parms.tmp_best_r_offsets(:,j,:);
  parms.tmp_best_th_offsets(:,i,:) = parms.tmp_best_th_offsets(:,j,:);
end;
if parms.offset_group_patches_flag % force patches to have same offsets
  for p=2:parms.npatches
    parms.tmp_best_r_offsets(:,:,p) = parms.tmp_best_r_offsets(:,:,1);
    parms.tmp_best_th_offsets(:,:,p) = parms.tmp_best_th_offsets(:,:,1);
  end;
end;
if parms.offset_group_areas_flag % force areas to have same offsets
  for a=2:parms.nareas
    parms.tmp_best_r_offsets(a,:,:) = parms.tmp_best_r_offsets(1,:,:);
    parms.tmp_best_th_offsets(a,:,:) = parms.tmp_best_th_offsets(1,:,:);
  end;
end;
for k=1:length(parms.offset_const_areas)
  % reset to best offsets for areas with constant offsets
  a = parms.offset_const_areas(k);
  parms.tmp_best_r_offsets(a,:,:) = parms.best_r_offsets(a,:,:);
  parms.tmp_best_th_offsets(a,:,:) = parms.best_th_offsets(a,:,:);
end;
parms.total_offset_mstarts = parms.total_offset_mstarts + 1;
