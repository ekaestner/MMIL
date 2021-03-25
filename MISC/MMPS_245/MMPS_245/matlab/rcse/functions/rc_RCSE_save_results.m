function parms = rc_RCSE_save_results(parms,avg_data,results,retforward,inverse,forward)
%function parms = rc_RCSE_save_results(parms,avg_data,results,retforward,inverse,forward)
%
% Purpose: save results to mat and optionally, fif files
%
% Required Input:
%   parms: RCSE parms struct
%   avg_data: average data struct
%   results: results struct
%   retforward: retinotopy-constrained forward struct
%   inverse: inverse solution struct
%   forward: forward solution struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/12/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,6), return; end;

% save cond_weights
if parms.condF_thresh>0
  results.cond_weights = parms.cond_weights;
  results.condF = parms.condF;
end;

% save inverse operator
matfile=sprintf('%s/matfiles/%s_inverse.mat',...
  parms.rootoutdir,parms.prefix);
fprintf('%s: saving inverse operator to %s...\n',mfilename,matfile);
save(matfile,'inverse');

if ~isfield(parms,'best_rot')
  parms.best_rot = zeros(3,1);
end;
if isfield(parms,'best_r_offset')
  parms.best_retmap.r_offset = parms.best_r_offset;
  parms.best_retmap.r_offsets = parms.r_offset;
end;
if isfield(parms,'best_th_offset')
  parms.best_retmap.th_offset = parms.best_th_offset;
  parms.best_retmap.th_offsets = parms.th_offset;
end;
if isfield(parms,'best_r_offsets')
  parms.best_retmap.r_offsets = parms.best_r_offsets;
end;
if isfield(parms,'best_th_offsets')
  parms.best_retmap.th_offsets = parms.best_th_offsets;
end;
if isfield(parms,'offset_errors')
  parms.best_retmap.offset_errors = parms.offset_errors;
end;
if isfield(parms,'offset_err_wforms')
  parms.best_retmap.offset_err_wforms = parms.offset_err_wforms;
end;
if isfield(parms,'offset_wforms')
  parms.best_retmap.offset_wforms = parms.offset_wforms;
end;

% save retintopy forward
retforward.best_rf_sizes = parms.best_rf_sizes;
retforward.best_rf_slopes = parms.best_rf_slopes;
retforward.best_rot = parms.best_rot;
retforward.best_rscale = parms.best_rscale;
retforward.r_max = parms.r_max;
retforward.full_retmap = parms.retmap;
retforward.retmap = parms.best_retmap;
retforward.min_error = parms.min_error;
retforward.total_offset_niters = parms.total_offset_niters;
retforward.total_offset_mstarts = parms.total_offset_mstarts;
retforward.total_rot_niters = parms.total_rot_niters;
retforward.total_rscale_niters = parms.total_rscale_niters;
retforward.total_rf_niters = parms.total_rf_niters;
retforward.total_nbrhd_niters = parms.total_nbrhd_niters;
matfile=sprintf('%s/matfiles/%s_ret_forward.mat',...
  parms.rootoutdir,parms.prefix);
fprintf('%s: saving retintopy forward and retmap to %s...\n',mfilename,matfile);
save(matfile,'retforward');

% extract source time courses
fprintf('%s: extracting source time courses...\n',mfilename);
results.retmap = parms.best_retmap;
retforward.min_error = parms.min_error;
for c=1:parms.ncontrasts
  j=1;
  for i=1:results.retmap.num_areas
    if parms.loose_flag
      tmp_data = results.S(:,j:3:j+3*results.retmap.num_locs-1,c);
      j=j+3*results.retmap.num_locs;
    elseif parms.indy_locs_flag
      tmp_data = results.S(:,j:j+results.retmap.num_locs-1,c);
      j=j+results.retmap.num_locs;
    else
      tmp_data = results.S(:,j,c);
      j=j+1;
    end
    if ~isfield(results.retmap.areas,'amplitudes') | isempty(results.retmap.areas(i).amplitudes)
      results.retmap.areas(i).amplitudes = zeros([size(tmp_data),parms.ncontrasts]);
    end;
    results.retmap.areas(i).amplitudes(:,:,c) = tmp_data;
    results.retmap.areas(i).time = avg_data.averages(1).time;
  end
  s0=j-1;

  % extract source waveforms for extra dips
  for d=1:length(results.retmap.ret_dips)
    if strcmp('lh',results.retmap.ret_dips(d).hemisphere)
      j = find(results.retmap.lh_ret_dips_v==results.retmap.ret_dips(d).v);
    elseif strcmp('rh',results.retmap.ret_dips(d).hemisphere)
      j = length(results.retmap.lh_ret_dips_v) +...
          find(results.retmap.rh_ret_dips_v==results.retmap.ret_dips(d).v);
    else
      j = [];
    end;
    if isempty(j)
      error('failed to extract ret_dips');
    end;
    s = s0 + (j-1)*3*results.retmap.num_ret_dip_locs + 1;
    % quarterfields: upper right, lower right, upper left, lower left
    for q=1:results.retmap.num_ret_dip_locs
      n = results.S(:,s,c); s=s+1;
      results.retmap.ret_dips(d).n(:,q,c) = n; % norm
      t1 = results.S(:,s,c); s=s+1;
      results.retmap.ret_dips(d).t1(:,q,c) = t1; % tang1
      t2 = results.S(:,s,c); s=s+1;
      results.retmap.ret_dips(d).t2(:,q,c) = t2; % tang2

      for t=1:num_tpoints
        results.retmap.ret_dips(d).a(t,q,c) = sqrt(n(t)^2 + t1(t)^2 + t2(t)^2);
      end
    end
  end
  s0 = s0 + results.retmap.num_ret_dips*3*results.retmap.num_ret_dip_locs;
  for d=1:length(results.retmap.nonret_dips)
    if strcmp('lh',results.retmap.nonret_dips(d).hemisphere)
      j = find(results.retmap.lh_nonret_dips_v==...
          results.retmap.nonret_dips(d).v);
    elseif strcmp('rh',results.retmap.nonret_dips(d).hemisphere)
      j = length(results.retmap.lh_nonret_dips_v) +...
          find(results.retmap.rh_nonret_dips_v==...
               results.retmap.nonret_dips(d).v);
    else
      j = [];
    end;
    s = s0 + (j-1)*3 + 1;
    n = results.S(:,s,c); s=s+1;
    results.retmap.nonret_dips(d).n(:,c) = n;
    t1 = results.S(:,s,c); s=s+1;
    results.retmap.nonret_dips(d).t1(:,c) = t1;
    t2 = results.S(:,s,c); s=s+1;
    results.retmap.nonret_dips(d).t2(:,c) = t2;
    for t=1:num_tpoints
      results.retmap.nonret_dips(d).a(t,c) = sqrt(n(t)^2 + t1(t)^2 + t2(t)^2);
    end
  end
end;

% save results
matfile=sprintf('%s/matfiles/%s_results.mat',...
  parms.rootoutdir,parms.prefix);
fprintf('%s: saving results to %s...\n',mfilename,matfile);
save(matfile,'results');

if(parms.write_err_flag)
  fprintf('%s: creating error structure...\n',mfilename);
  tmp_parms = parms;
  tmp_parms.fiterr_flag = 1;
  tags = {'prefix','fiterr_flag'};
  args = mmil_parms2args(tmp_parms,tags);
  err_data = rc_synth_sensors_from_RCSE(args{:});
  matfile=sprintf('%s/matfiles/%s_err.mat',...
    parms.rootoutdir,parms.prefix);
  save(matfile,'err_data');

  % write to fif
  if(parms.write_fif_flag)
    outstem=sprintf('fifs/%s_err',parms.prefix);
    fprintf('%s: writing residual error to fif...\n',mfilename);
    ts_avg2fif(err_data,parms.template_fif,outstem,parms.evcode_flag,1);
  end
end;

if parms.write_fit_flag
  % create fitted average structures
  fprintf('%s: creating fitted average structure...\n',mfilename);
  tmp_parms = parms;
  tmp_parms.fiterr_flag = 0;
  tmp_parms.baselineflag = 0;
  tags = {'prefix','fiterr_flag','baselineflag'};
  args = mmil_parms2args(tmp_parms,tags);
  fit_data = rc_synth_sensors_from_RCSE(args{:});
  matfile=sprintf('%s/matfiles/%s_fit.mat',...
    parms.rootoutdir,parms.prefix);
  save(matfile,'fit_data');

  % write to fif
  if(parms.write_fif_flag)
    outstem=sprintf('fifs/%s_fit',parms.prefix);
    fprintf('%s: writing expected average to fif...\n',mfilename);
    ts_avg2fif(fit_data,parms.template_fif,outstem,parms.evcode_flag,1);
  end
end;

if parms.write_areafit_flag
  % create fitted average for each area
  for a=1:results.retmap.num_areas
    fprintf('%s: creating fitted average structure for area %d...\n',...
      mfilename,a);
    tmp_parms = parms;
    tmp_parms.fiterr_flag = 0;
    tmp_parms.baselineflag = 0;
    tmp_parms.areas = a;
    tags = {'prefix','fiterr_flag','baselineflag','areas'};
    args = mmil_parms2args(tmp_parms,tags);
    fit_data = rc_synth_sensors_from_RCSE(args{:});
    matfile=sprintf('%s/matfiles/%s_fit_%s.mat',...
      parms.rootoutdir,parms.prefix,results.retmap.areas(a).name);
    save(matfile,'fit_data');

    % write to fif
    if(parms.write_fif_flag)
      outstem=sprintf('fifs/%s_%s_fit',...
        parms.prefix,results.retmap.areas(a).name);
      fprintf('%s: writing expected average to fif...\n',mfilename);
      ts_avg2fif(fit_data,parms.template_fif,outstem,parms.evcode_flag,1);
    end
  end;
end;

%% todo: write to stc file (from retmap)
% for each condition, set values for each vertex in uniq_verts
% for areas, if not indy, set all vertices to same value
% for non-ret dips, set all vertices to same value
% for ret dips, set vertex for each dip to correct value
%   with qfield_flag, set quarterfield vertices to value

