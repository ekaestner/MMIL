function parms = rc_RCSE_full_offset_search(parms,avg_data,forward)
%function parms = rc_RCSE_full_offset_search(parms,avg_data,forward)
%
% Purpose: dipole optimization using specified vectors of
%   r and th offsets (exhaustive search)
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

parms.offset_errors = zeros(parms.num_r_offsets,parms.num_th_offsets);
parms.offset_err_wforms = [];
parms.offset_wforms = [];
fprintf('%s: finding best retfit offsets...\n',mfilename);
for i=1:parms.num_r_offsets
  tmp_r_offsets = ...
    parms.r_offset(i)*ones(parms.nareas,parms.nconds,parms.npatches);
  for k=1:length(parms.offset_const_areas)
    % keep offsets constant (set to zero)
    a = parms.offset_const_areas(k);
    tmp_r_offsets(a,:,:) = 0;
  end;
  tmp_r_offsets = tmp_r_offsets + parms.init_r_offsets;
  for j=1:parms.num_th_offsets
    tmp_th_offsets = ...
      parms.th_offset(j)*ones(parms.nareas,parms.nconds,parms.npatches);
    for k=1:length(parms.offset_const_areas)
      % keep offsets constant (set to zero)
      a = parms.offset_const_areas(k);
      tmp_th_offsets(a,:,:) = 0;
    end;
    tmp_th_offsets = tmp_th_offsets + parms.init_th_offsets;
    retmap = rc_RCSE_retmap_from_retfit(parms,tmp_r_offsets,tmp_th_offsets,[],forward);
    [parms,results,retforward,inverse]=...
      rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
    fprintf('%s: r_offset = %0.4f, th_offset = %0.4f, error = %f\n',...
      mfilename,parms.r_offset(i),parms.th_offset(j),parms.total_error);
    % if error is reduced, use current dipoles as benchmark
    if ~isnan(parms.total_error) &...
       ((parms.total_error<parms.min_error) | isempty(parms.best_retmap))
      parms.min_error = parms.total_error;
      parms.retmap = retmap;
      parms.best_retmap = retforward.retmap;
      parms.best_r_offset = parms.r_offset(i);
      parms.best_th_offset = parms.th_offset(j);
      fprintf('%s: ### full offset search better fit\n',mfilename);
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
    % store error and waveforms for each offset
    parms.offset_errors(i,j) = parms.total_error;
    err_wform = results.norm_var_E;
    if isempty(parms.offset_err_wforms)
      matsize = [parms.num_r_offsets,parms.num_th_offsets,length(err_wform)];
      parms.offset_err_wforms = zeros(matsize);
    end;
    parms.offset_err_wforms(i,j,:) = err_wform;
    wforms = squeeze(rc_RCSE_extract_source_wforms(parms,avg_data,results));
    if isempty(parms.offset_wforms)
      matsize = [parms.num_r_offsets,parms.num_th_offsets,size(wforms)];
      parms.offset_wforms = zeros(matsize);
    end;
    parms.offset_wforms(i,j,:) = reshape(wforms,[1 prod(size(wforms))]);
  end;
end;
fprintf('%s: ### full offset search min error = %0.4f\n',mfilename,parms.min_error);
fprintf('%s: best r_offset = %0.4f, best th_offset = %0.4f\n',...
  mfilename,parms.best_r_offset,parms.best_th_offset);

