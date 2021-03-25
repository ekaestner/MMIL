function [parms,results,retforward,inverse] = rc_RCSE_reweight_lsq(parms,retmap,avg_data,forward)
%function [parms,results,retforward,inverse] = rc_RCSE_reweight_lsq(parms,retmap,avg_data,forward)
%
% Purpose: robust estimation using iteratively reweighted least squares
%
% Required Input:
%   parms: RCSE parms struct
%   retmap: retinotopy mapping struct (if empty, will generate one)
%   avg_data: average data struct
%   forward: forward solution struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 03/29/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

% initialize retmap
if isempty(retmap) && ~isempty(parms.retfit_results)
  retmap = ...
    rc_RCSE_retmap_from_retfit(parms,parms.r_offset,parms.th_offset,[],forward);
end;

% calculate initial source estimates
fprintf('%s: calculating initial source estimates...\n',mfilename);
[parms,results,retforward,inverse]=...
  rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward);
parms.best_retmap = results.retmap;
wf0 = rc_RCSE_extract_source_wforms(parms,avg_data,results);

orig_retforward = retforward;
orig_avg_data = avg_data;

num_tpoints = size(results.E,1);
num_sensors = length(parms.goodchans);
num_measurements = size(retforward.F,1);

tiny_s = 1e-6 * max(std(results.Y,0,2));
if tiny_s==0
  tiny_s = 1;
end

% adjust residuals using "leverage"
if parms.reweight_leverage_flag
  X = retforward.F;
  % from Barret and Gray, 1997
  H = X*inv(X'*X)*X';
  % from matlab statrobustfit
  h = min(.9999, sum(H.*H,2));
  adjfactor = 1 ./ sqrt(1-h);
  num_stimlocs = num_measurements/num_sensors;
  adjfactor = reshape(adjfactor,[num_sensors,num_stimlocs]);
  if parms.reweight_leverage_max_flag
    adjfactor = max(adjfactor,[],1);
  else
    adjfactor = mean(adjfactor,1);
  end;
else
  adjfactor = 1;
end;

% iteratively reweighted least squares
for rw_iter=1:parms.reweight_maxiter
  fprintf('%s: reweight iter %d\n',mfilename,rw_iter);
  if parms.plotflag
    fprintf('%s: plotting...\n',mfilename);
    figure(1);
    rc_plot_sources('parms',parms,'results',results,...
      'area_names',parms.area_names,'area_colors',parms.area_colors,...
      'ylim',parms.plot_ylim,'units',parms.source_units,...
      'labelflag',0);
    figure(2);
    rc_plot_fitvar_results(parms,results,1);
    drawnow;
  end;

  % calculate average abs(error) for each condition
  if parms.fit_range_flag
    fit_range = [parms.fit_t0:parms.fit_t1];
  else
    fit_range = [1:num_tpoints];
  end;

  errmat = zeros(parms.ncontrasts,retforward.retmap.num_locs);
  j = 1;
  for c=1:parms.ncontrasts
    for i=1:retforward.retmap.num_locs
      % calculate average error for this condition
      k = j+num_sensors-1;
      err = results.E(fit_range,j:k);
      errmat(c,i) = mean(abs(err(:)));
      j = k+1;
    end;
  end
  clear err;
  errvec = mean(errmat,1);

  % calculate deviation from minimum error
  err_diff = errvec - min(errvec);

  % adjust for "leverage"
  err_diff = err_diff .* adjfactor;

  % calculate median absolute deviation
  s = max(median(err_diff) / 0.6745,tiny_s);

  % offset and normalize errors
  weights = err_diff/(s*parms.reweight_factor);
  err_norm = err_diff/s;

  % apply Tukey's bisquare weight function
  %   (greater error than median)
  weights = sqrt((abs(weights<1)) .* (1 - weights.^2).^2);

  if parms.plotflag %% todo: save these
    [err_sorted,ind_sort] = sort(err_norm);
    weights_sorted = weights(ind_sort);

    figure(3); clf;
    hist(err_norm)
    title('error histogram');

    figure(4); clf;
    plot(err_sorted,weights_sorted,'g^');
    set(gca,'YLim',[0,1]);
    xlabel('stimulus location (sorted by error)');
    title('weights vs. error');

    figure(5); clf;
    plot(err_sorted,'ro');
    title('normalized error');
    set(gca,'YLim',[-3,3]);

    figure(6); clf;
    plot(weights_sorted,'b*');
    set(gca,'YLim',[0,1]);
    xlabel('stimulus location (sorted by error)');
    title('weights');
    drawnow;
  end;

  % recalculate source estimates
  fprintf('%s: recalculating source estimates...\n',mfilename);
  parms.cond_weights = weights;
  [parms,results,retforward,inverse]=...
    rc_RCSE_fit_waveforms(parms,retmap,avg_data,forward,retforward);
  wf = rc_RCSE_extract_source_wforms(parms,avg_data,results);
  results.cond_err = errvec;
  results.cond_err_norm = err_norm;

  if rw_iter==1
    init_weights = weights;
    init_errvec = errvec;
    init_err_norm = err_norm;
  end;

  % check for convergence
  wf_diff = wf-wf0;
  wf_diff = mean(abs(wf_diff(:)))/max(abs(wf0(:)));
  if wf_diff<parms.reweight_tol
    fprintf('%s: difference in estimates is less than tolerance\n',mfilename);
    break;
  end;
  wf0 = wf;
end;

if rw_iter>=parms.reweight_maxiter
  fprintf('%s: WARNING: reweighted least squares did not converge within %d iterations\n',...
    mfilename,parms.reweight_maxiter);
end;

results.cond_weights = parms.cond_weights;
results.init_cond_weights = init_weights;
results.init_cond_err = init_errvec;
results.init_cond_err_norm = init_err_norm;

