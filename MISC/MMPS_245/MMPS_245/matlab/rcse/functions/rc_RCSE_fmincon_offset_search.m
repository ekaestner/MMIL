function [parms,results] = rc_RCSE_fmincon_offset_search(parms,avg_data,forward)
%function [parms,results] = rc_RCSE_fmincon_offset_search(parms,avg_data,forward)
%
% Purpose: optimize dipole offset using MATLAB's fmincon
%
% Required Input:
%   parms: RCSE parms struct
%   avg_data: average data struct
%   forward: forward solution struct
%
% Created:  05/13/11 by Don Hagler
% Last Mod: 10/24/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

TolFun = 1e-4;
TolX = 1e-5;
DiffMinChange = parms.r_step;
DiffMaxChange = 10*DiffMinChange;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try % for R2009b
  fit_options = optimset('Display','iter',...
    'Algorithm','active-set' ,...
    'MaxFunEvals',Inf,...
    'MaxIter',parms.offset_niters,...
    'DiffMinChange',DiffMinChange,...
    'DiffMaxChange',DiffMaxChange,...
    'TolFun',TolFun,...
    'TolX',TolX);
catch % for R2007a
  fit_options = optimset('Display','iter',...
    'LargeScale','off',...
    'MaxFunEvals',Inf,...
    'MaxIter',parms.offset_niters,...
    'DiffMinChange',DiffMinChange,...
    'DiffMaxChange',DiffMaxChange,...
    'TolFun',TolFun,...
    'TolX',TolX);
end;

[betas,lbounds,ubounds] = init_betas(parms);

% extract prior waveform from mat file and resample to desired time
%   and number of waveforms
%   (e.g. to go from multiple contrasts to a single waveform)
if isempty(parms.fname_prior)
  parms.prior = [];
elseif ~isfield(parms,'prior') || isempty(parms.prior)
  fprintf('%s: extracting prior...\n',mfilename);
  tic;
  time = avg_data.averages(1).time;
  parms.prior = rc_extract_prior(parms.fname_prior,...
    time,parms.use_areas,parms.unique_contrasts,parms.prior_avg_flag);
  toc;
end;

[betas,min_err,exitflag,output] = ...
  fmincon(@(b) calc_err(b,parms,avg_data,forward),betas,...
  [],[],[],[],lbounds,ubounds,[],fit_options);

if exitflag~=1
  fprintf('%s: WARNING: fmincon returned with exitflag %d\n',...
    mfilename,exitflag);
end;

[err,parms,results,retforward,inverse,cond_info] = calc_err(betas,parms,avg_data,forward);
[r_offsets,th_offsets] = betas2offsets(betas,parms);
parms.tmp_min_error = err;
parms.tmp_best_retmap = retforward.retmap;
parms.tmp_best_r_offsets = r_offsets;
parms.tmp_best_th_offsets = th_offsets;
parms.tmp_best_cond_info = cond_info;

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

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err,parms,results,retforward,inverse,cond_info] = calc_err(betas,parms,avg_data,forward)
  [r_offsets,th_offsets] = betas2offsets(betas,parms);
  cond_info = rc_RCSE_set_cond_info_offsets(parms,r_offsets,th_offsets);
  retmap = rc_RCSE_retmap_from_retfit(parms,0,0,cond_info,forward);
  [parms,results,retforward,inverse]=...
    rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
  err = parms.total_error;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [betas,lbounds,ubounds] = init_betas(parms)
  % initialize betas vector with values from r_offsets and th_offsets
  r_offsets = parms.tmp_best_r_offsets;
  th_offsets = parms.tmp_best_th_offsets;
  unique_offset_conds = unique(parms.offset_group_conds)';
  offset_nconds = length(unique_offset_conds);
  if parms.offset_group_patches_flag
    offset_npatches = 1;
  else
    offset_npatches = parms.npatches;
  end;
  offset_areas = setdiff(parms.use_areas,parms.offset_const_areas);
  if parms.offset_group_areas_flag
    offset_areas = offset_areas(1);
  end;
  offset_nareas = length(offset_areas);

  nbetas = 2 * offset_nconds * offset_npatches * offset_nareas;
  betas = zeros(1,nbetas);
  lbounds = zeros(1,nbetas);
  ubounds = zeros(1,nbetas);
  j = 1;
  for a=offset_areas
    for k=1:offset_nconds
      c = unique_offset_conds(k);
      for p=1:offset_npatches
        betas(j) = r_offsets(a,c,p);
        lbounds(j) = parms.r_offset_range(1);
        ubounds(j) = parms.r_offset_range(2);
        j = j + 1;
        betas(j) = th_offsets(a,c,p);
        lbounds(j) = parms.th_offset_range(1);
        ubounds(j) = parms.th_offset_range(2);
        j = j + 1;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_offsets,th_offsets] = betas2offsets(betas,parms)
  % set r_offsets and th_offsets from betas

  r_offsets = zeros(parms.nareas,parms.nconds,parms.npatches);
  th_offsets = zeros(parms.nareas,parms.nconds,parms.npatches);

  unique_offset_conds = unique(parms.offset_group_conds);
  offset_nconds = length(unique_offset_conds);
  if parms.offset_group_patches_flag
    offset_npatches = 1;
  else
    offset_npatches = parms.npatches;
  end;
  offset_areas = setdiff(parms.use_areas,parms.offset_const_areas);
  if parms.offset_group_areas_flag
    offset_areas = offset_areas(1);
  end;
  offset_nareas = length(offset_areas);

  j = 1;
  for a=offset_areas
    for k=1:offset_nconds
      c = unique_offset_conds(k);
      for p=1:offset_npatches
        r_offsets(a,c,p) = betas(j);
        j = j + 1;
        th_offsets(a,c,p) = betas(j);
        j = j + 1;
      end;
    end;
  end;

  % force some conditions to have same offsets (e.g. within quad)
  for i=1:parms.nconds
    j=parms.offset_group_conds(i);
    if i==j, continue; end;
    r_offsets(:,i,:) = r_offsets(:,j,:);
    th_offsets(:,i,:) = th_offsets(:,j,:);
  end;
  if parms.offset_group_patches_flag % force patches to have same offsets
    for p=2:parms.npatches
      r_offsets(:,:,p) = r_offsets(:,:,1);
      th_offsets(:,:,p) = th_offsets(:,:,1);
    end;
  end;
  if parms.offset_group_areas_flag % force areas to have same offsets
    for a=2:parms.nareas
      r_offsets(a,:,:) = r_offsets(offset_areas(1),:,:);
      th_offsets(a,:,:) = th_offsets(offset_areas(1),:,:);
    end;
  end;
  for k=1:length(parms.offset_const_areas)
    % keep offsets constant (set to zero)
    a = parms.offset_const_areas(k);
    r_offsets(a,:,:) = 0;
    th_offsets(a,:,:) = 0;
  end;

return;

