function results = MEG_MMIL_RCSE_GroupAvg_ErrMap_Resamp(ProjID,varargin)
%function results = MEG_MMIL_RCSE_GroupAvg_ErrMap_Resamp(ProjID,[options])
%
% Purpose: compare maps of RCSE residual error across locations and subjects
%   using bootstrap resampling
%
% Required input:
%   ProjID: project ID string
%
% Optional parameters:
%   'prefix': RCSE prefix
%     {default = 'RCSE'}
%   'prefix2': second RCSE prefix to be compared to first
%     This function assumes that the two prefixes correspond to
%     mutually exclusive sets of stimulus locations (e.g. upper and lower field)
%     And currently requires that both sets have identical number of locations
%     {default = []}
%   'outdir': output directory
%     if not supplied, will place in
%       /home/<user>/MetaData/<ProjID>/Resamp_RCSE_ErrMap
%     {default = []}
%   'outstem': output file stem
%     {default = 'RCSE_ErrMap'}
%   'outfix': output file suffix
%     {default = 'err'}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Optional parameters for error maps
%   't0': start time of range used to calculate error
%     {default = -Inf}
%   't1': end time of range used to calculate error
%     {default = Inf}
%   'sign_flag': [0|1|2] whether to use abs, signed, or squared error
%     0: absolute error
%     1: signed error
%     2: squared error
%     {default = 0}
%   'std_flag': [0|1] calculate standard deviation of error across sensors
%     otherwise calculate mean
%     {default = 0} 
%   'irls_err_flag': [0|1] whether to plot normalized error calculatd
%      during iteratively reweighted least-squares (IRLS)
%      otherwise, calculate error from results.E
%     {default = 0}
%   'irls_weights_flag': [0|1] whether to plot weights from IRLS instead of err
%     {default = 0}
%   'irls_init_flag': [0|1] whether IRLS err or weights should be from initial
%     iteration (1) or last (0)
%     {default = 0}
%   'renorm_flag': [0|1|2] whether to renormalize error
%     (e.g. after combining prefix and prefix2 results)
%     0: use existing normalization
%     1: normalize to max of max_std between prefix and prefix2
%     2: normalize to median of max_std across subjects
%     {default = 1}
%
% Optional parameters for plotting:
%   'plot_flag': [0|1] whether to make plots
%     {default = 1}
%   'visible_flag': [0|1] whether plots should be visible
%     {default = 1}
%   'tif_flag': [0|1] save plots as tifs
%     {default = 1}
%   'eps_flag': [0|1] save plots as eps
%     {default = 0}
%   'err_xlim': two element vector of x-axis limits for error histograms
%     {default = [0 1]}
%   'err_ylim': two element vector of y-axis limits for error plots
%     {default = [0 1]}
%   'clim': two element vector of limits for error maps
%     {default = [0 1]}
%   'zlim': two element vector of limits for z-transformed error maps
%     {default = [-1 1]}
%   'cmap': name of colormap for error maps
%     {default = 'bone'}
%
% Other optional parameters:
%   'user': user name for account with ProjInfo
%     if not supplied, will use getenv('USER')
%     {default = []}
%   'niters': number of iterations
%     {default = 2000}
%   'randseed_flag': [0|1] generate different result every time
%     {default = 0}
%   'ci_tail': confidence interval tail size (%)
%     {default = 2.5}
%
% Output:
%   results: struct containing results of calculations
%
% Created:  03/20/12 by Don Hagler
% Last Mod: 10/14/13 by Don Hagler
%

%% todo: change to not require ProjID
%%       NOTE: requires change to MEG_MMIL_RCSE_GroupAvg to store E

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'prefix','RCSE',[],...
  'prefix2',[],[],...
  'outdir',[],[],...
  'outstem','RCSE_ErrMap',[],...
  'outfix','err',[],...
  'forceflag',false,[false true],...
... % error maps
  't0',-Inf,[],...
  't1',Inf,[],...
  'sign_flag',0,[0 1 2],...
  'std_flag',false,[false true],...
  'irls_err_flag',false,[false true],...
  'irls_weights_flag',false,[false true],...
  'irls_init_flag',false,[false true],...
  'renorm_flag',1,[0 1 2],...
... % plotting
  'plot_flag',true,[false true],...  
  'visible_flag',true,[false true],...
  'tif_flag',true,[false true],...
  'eps_flag',false,[false true],...
  'err_xlim',[0,1],[],...
  'err_ylim',[0,1],[],...
  'clim',[0,1],[],...
  'zlim',[-1,1],[],...
  'cmap','bone',[],...
... % more for plotting
  'label_flag',true,[false true],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
  'grid_unit',0.05,[],...
  'z_offset',0,[],...
  'cbar_flag',true,[false true],...
  'fg_color',0.4,[],...
  'bg_color',0,[],...
  'p1_color','k',[],...
  'p2_color','b',[],...
... % more
  'user',getenv('USER'),[],...
  'niters',2000,[1 Inf],...
  'randseed_flag',false,[false true],...
  'ci_tail',2.5,[],...
... % tags
  'errmap_tags',{'rootdir' 't0' 't1' ...
                 'sign_flag' 'std_flag' 'irls_err_flag' 'irls_weights_flag'...
                 'irls_init_flag' 'tif_flag' 'eps_flag' 'visible_flag' 'outdir' 'outstem'...
                 'cmap' 'label_flag' 'fontname' 'fontsize' 'grid_unit' 'clim' 'z_offset'...
                 'err_xlim' 'cbar_flag' 'fg_color' 'bg_color' 'plot_flag'},[],...
  'imagesc_tags',{'outdir' 'outstem' 'tif_flag' 'eps_flag' 'visible_flag' 'cbar_flag'...
                  'clim' 'cmap' 'title' 'fontname' 'fontsize' 'offset'},[],...
  'errvec_tags',{'outdir','outstem','tif_flag','eps_flag','visible_flag',...
                 'title','fontname','fontsize','wforms_err','xlabel','ylabel','time',...
                 'ylim','colors','points','condnames','legend_loc'},[],...
};
results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(ProjID,varargin,parms_filter);

% calculate RCSE per-location error for each subject
calc_subj_errmaps(parms);

matfile = sprintf('%s/%s_%s.mat',parms.matdir,parms.outstem,parms.outfix);
if ~exist(matfile,'file') || parms.forceflag
  fprintf('%s: calculating averages...\n',mfilename);

  % load and concatenate indiividual subject error
  results.data = load_data(parms);

  % calculate mean and SEM
  results = calc_avg_err(parms,results);

  % calculate difference in error between prefix and prefix2 results
  if ~isempty(parms.prefix2)
    results = calc_err_diff(parms,results);
  end;

  fprintf('%s: resampling...\n',mfilename);

  % use bootstrap resampling to calculate confidence intervals
  if parms.niters
    results = resamp_err(parms,results);
  end;

  % save results
  save(matfile,'results','-v7.3');
else
  load(matfile);
end;

if parms.plot_flag
  fprintf('%s: plotting...\n',mfilename);

  % plot maps of average error
  plot_errmaps(parms,results);

  % scatter plot of error for each location, sorted by error
  plot_err_sorted(parms,results);
  
  if ~isempty(parms.prefix2)
    % scatter plot of error for each location, unsorted by error
    plot_err_unsorted(parms,results);
  end;
end;

if ~isempty(parms.prefix2)
  % write t-test and ANOVA results to text file
  write_diff_text(parms,results);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  [parms.StudyInfo,parms.RootDirs] = ...
    MMIL_Get_StudyInfo(ProjID,'user',parms.user);
  if isempty(parms.StudyInfo)
    error('StudyInfo is missing for ProjID %s',ProjID);
  end;
  parms.nsubs = length(parms.StudyInfo);
  if isempty(parms.outdir)
    parms.outdir = sprintf('/home/%s/MetaData/%s/Resamp_RCSE_ErrMap',...
      parms.user,ProjID);
  end;
  parms.matdir = [parms.outdir '/matfiles'];
  parms.subjdir = [parms.outdir '/err_subj'];
  parms.resdir = [parms.outdir '/err_resamp'];

  if isempty(parms.prefix2)
    parms.prefix_list = {parms.prefix};
  else
    parms.prefix_list = {parms.prefix,parms.prefix2};
  end;
  parms.nprefix = length(parms.prefix_list);

  if parms.randseed_flag
    seedval = sum(100*clock);
  else
    seedval = 5489; % default
  end;
  stream = RandStream.create('mt19937ar','Seed',seedval);
  RandStream.setDefaultStream(stream);

  mmil_mkdir(parms.matdir);
  mmil_mkdir(parms.subjdir);
  mmil_mkdir(parms.resdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate RCSE per-location error for each subject
function calc_subj_errmaps(parms)
  for s=1:parms.nsubs
    VisitID = parms.StudyInfo(s).MEG_VisitID;
    parms.rootdir = sprintf('%s/%s',...
      parms.RootDirs.proc_meg,parms.StudyInfo(s).proc_meg);
    for p=1:parms.nprefix
      prefix = parms.prefix_list{p};
      matfile=sprintf('%s/%s_%s_%s.mat',...
        parms.matdir,VisitID,prefix,parms.outfix);
      if ~exist(matfile,'file') || parms.forceflag
        fprintf('%s: calculating RCSE error map for %s %s...\n',...
          mfilename,VisitID,prefix);
        tmp_parms = parms;
        tmp_parms.outdir = parms.subjdir;
        tmp_parms.outstem = ...
          sprintf('%s_%s_%s',VisitID,prefix,parms.outfix);
        args = mmil_parms2args(tmp_parms,parms.errmap_tags);
        output = rc_plot_errmap(prefix,args{:});
        save(matfile,'output');
        fprintf('%s: finished with %s %s\n',mfilename,VisitID,prefix);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load and concatenate indiividual subject error
function data = load_data(parms)
  data = init_data(parms);
  for s=1:parms.nsubs
    VisitID = parms.StudyInfo(s).MEG_VisitID;
    tmpfile=sprintf('%s/%s_%s_%s.mat',...
      parms.matdir,VisitID,parms.prefix,parms.outfix);
    tmp = load(tmpfile);
    max_std = tmp.output.max_std;
    max_var = tmp.output.max_var;
    err_norm = tmp.output.err_norm';
    err_raw = tmp.output.err_raw';
    Z = tmp.output.Z;
    W = tmp.output.W;
    M = tmp.output.M;

    % load data for second prefix and combine with first
    if ~isempty(parms.prefix2)
      tmpfile=sprintf('%s/%s_%s_%s.mat',...
        parms.matdir,VisitID,parms.prefix2,parms.outfix);
      tmp = load(tmpfile);
      Z2 = tmp.output.Z;
      W2 = tmp.output.W;
      M2 = tmp.output.M;
      err_norm2 = tmp.output.err_norm';
      err_raw2 = tmp.output.err_raw';
      max_std2 = tmp.output.max_std;
      max_var2 = tmp.output.max_var;

      % combine err maps
      Z(M2>0) = Z2(M2>0);
      W(M2>0) = W2(M2>0);
      M(M2>0) = -1*M2(M2>0);

      % combine err
      err_norm1 = err_norm;
      err_raw1 = err_raw;
      max_std1 = max_std;
      max_var1 = max_var;
      err_norm = [err_norm1;err_norm2];
      err_raw = [err_raw1;err_raw2];
      max_std = max(max_std1,max_std2);
      max_var = max(max_var1,max_var2);

      % renormalize error
      if parms.renorm_flag == 1
        if parms.sign_flag == 2
          err_norm1 = err_raw1 / max_var;
          err_norm2 = err_raw2 / max_var;
        else
          err_norm1 = err_raw1 / max_std;
          err_norm2 = err_raw2 / max_std;
        end;
      end;

      data.set1.err_norm_vals = cat(2,data.set1.err_norm_vals,err_norm1);
      data.set1.err_raw_vals = cat(2,data.set1.err_raw_vals,err_raw1);
      data.set1.max_std_vals = cat(2,data.set1.max_std_vals,max_std1);
      data.set1.max_var_vals = cat(2,data.set1.max_var_vals,max_var1);

      data.set2.err_norm_vals = cat(2,data.set2.err_norm_vals,err_norm2);
      data.set2.err_raw_vals = cat(2,data.set2.err_raw_vals,err_raw2);
      data.set2.max_std_vals = cat(2,data.set2.max_std_vals,max_std2);
      data.set2.max_var_vals = cat(2,data.set2.max_var_vals,max_var2);
    end;

    % renormalize error
    if parms.renorm_flag == 1
      if parms.sign_flag == 2
        Z = W / max_var;
        err_norm = err_raw / max_var;
      else
        Z = W / max_std;
        err_norm = err_raw / max_std;
      end;
    end;

    % concatenate err maps across subjects
    data.Z_vals = cat(3,data.Z_vals,Z);
    data.W_vals = cat(3,data.W_vals,W);
    data.M_vals = cat(3,data.M_vals,M);

    % concatenate err vectors across subjects
    data.err_norm_vals = cat(2,data.err_norm_vals,err_norm);
    data.err_raw_vals = cat(2,data.err_raw_vals,err_raw);
    data.max_std_vals = cat(2,data.max_std_vals,max_std);
    data.max_var_vals = cat(2,data.max_var_vals,max_var);
  end;
  
  data.nlocs = size(data.err_raw_vals,1);
  data.max_std = median(data.max_std_vals);
  data.max_var = median(data.max_var_vals);

  if ~isempty(parms.prefix2)
    data.set1.nlocs = size(data.set1.err_raw_vals,1);
    data.set2.nlocs = size(data.set2.err_raw_vals,1);
  end;
  %% todo: replace uses of nlocs2 with set1.nlocs and set2.nlocs
  data.nlocs2 = data.nlocs/2;

  if parms.renorm_flag == 2
    if parms.sign_flag == 2
      data.Z_vals = data.W_vals / data.max_var;
      data.err_norm_vals = data.err_raw_vals / data.max_var;
      if ~isempty(parms.prefix2)
        data.set1.err_norm_vals = data.set1.err_raw_vals / data.max_var;
        data.set2.err_norm_vals = data.set2.err_raw_vals / data.max_var;
      end;
    else
      data.Z_vals = data.W_vals / data.max_std;
      data.err_norm_vals = data.err_raw_vals / data.max_std;
      if ~isempty(parms.prefix2)
        data.set1.err_norm_vals = data.set1.err_raw_vals / data.max_std;
        data.set2.err_norm_vals = data.set2.err_raw_vals / data.max_std;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = init_data(parms)
  data = [];
  data.nsubs = parms.nsubs;
  data.Z_vals = [];
  data.W_vals = [];
  data.M_vals = [];
  data.err_norm_vals = [];
  data.err_raw_vals = [];
  data.max_std_vals = [];
  data.max_var_vals = [];
  if ~isempty(parms.prefix2)
    data.set1 = [];
    data.set1.err_norm_vals = [];
    data.set1.err_raw_vals = [];
    data.set1.max_std_vals = [];
    data.set1.max_var_vals = [];
    data.set2 = [];
    data.set2.err_norm_vals = [];
    data.set2.err_raw_vals = [];
    data.set2.max_std_vals = [];
    data.set2.max_var_vals = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate mean and SEM of error
function results = calc_avg_err(parms,results)
  avg_err = [];
  results.Z = mean(results.data.Z_vals,3);
  results.W = mean(results.data.W_vals,3);
  results.M = mean(results.data.M_vals,3);
  results.err_raw =  mean(results.data.err_raw_vals,2);
  results.err_norm =  mean(results.data.err_norm_vals,2);
  results.Z_sem  = std(results.data.Z_vals,0,3) / sqrt(results.data.nsubs);
  results.W_sem  = std(results.data.W_vals,0,3) / sqrt(results.data.nsubs);
  results.err_raw_sem = std(results.data.err_raw_vals,[],2) / sqrt(results.data.nsubs);
  results.err_norm_sem = std(results.data.err_norm_vals,[],2) / sqrt(results.data.nsubs);

  % ANOVA comparing stimulus locations
  errvals = results.data.err_norm_vals';
  [results.anova1_p,results.anova1_table,results.anova1_stats] = anova1(errvals,[],'off');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_err_diff(parms,results)
  errvals = results.data.err_norm_vals;

  % 2-way ANOVA with prefix and location as factors
  tmp_errvals = zeros(results.data.nsubs*results.data.nlocs2,2);
  j = 1;
  for i=1:results.data.nlocs2
    k = j + results.data.nsubs - 1;
    tmp_errvals(j:k,1) = errvals(i,:);
    tmp_errvals(j:k,2) = errvals(results.data.nlocs2+i,:);
    j = k + 1;
  end;
  [results.anova2_p,results.anova2_table,results.anova2_stats] = ...
    anova2(tmp_errvals,results.data.nsubs,'off');

  % difference in mean err between prefix and prefix2 conditions
  errvals1 = results.data.set1.err_norm_vals;
  errvals2 = results.data.set2.err_norm_vals;
  err1 = mean(errvals1,1);
  err2 = mean(errvals2,1);
  err_diff = err1 - err2;
  results.err_diff_mean = mean(err_diff,2);
  results.err_diff_std = std(err_diff,0,2);
  results.err_diff_sem = results.err_diff_std / sqrt(results.data.nsubs);
  [~,pval,ci,stats] = ttest(err_diff);
  results.err_diff_pval = pval;
  results.err_diff_ttest_ci = ci;
  results.err_diff_tstat = stats.tstat;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = resamp_err(parms,results)
  % resample across subjects
  ind_resamp = reshape(randi(results.data.nsubs,results.data.nsubs,parms.niters),...
    [1,results.data.nsubs*parms.niters]);

  % resample err vals across subjects
  errvals = results.data.err_norm_vals;
  errvals_resamp = reshape(errvals(:,ind_resamp),...
    [results.data.nlocs,results.data.nsubs,parms.niters]);
  errvals_resamp_mean = squeeze(mean(errvals_resamp,2));
  results.err_resamp_mean = mean(errvals_resamp_mean,2);
  results.err_resamp_lo = prctile(errvals_resamp_mean,parms.ci_tail,2);
  results.err_resamp_hi = prctile(errvals_resamp_mean,100-parms.ci_tail,2);
  results.err_resamp_ci = ...
    cat(2,results.err_resamp_lo,results.err_resamp_hi);
  results.err_resamp_sem = std(errvals_resamp_mean,1,2);

  % resample mean err (across locations) across subjects
  errvals_meanloc = mean(errvals,1);
  errvals_meanloc_resamp = reshape(errvals_meanloc(:,ind_resamp),...
    [results.data.nsubs,parms.niters]);
  errvals_meanloc_resamp_mean = mean(errvals_meanloc_resamp,1);
  errvals_meanloc_resamp_std = std(errvals_meanloc_resamp,1,1);
  results.err_meanloc_resamp_mean = mean(errvals_meanloc_resamp_mean,2);
  results.err_meanloc_resamp_lo = prctile(errvals_meanloc_resamp_mean,parms.ci_tail,2);
  results.err_meanloc_resamp_hi = prctile(errvals_meanloc_resamp_mean,100-parms.ci_tail,2);
  results.err_meanloc_resamp_ci = cat(2,...
    results.err_meanloc_resamp_lo,results.err_meanloc_resamp_hi);
  results.err_meanloc_resamp_sem = std(errvals_meanloc_resamp_mean,1,2);
  results.err_meanloc_resamp_std = mean(errvals_meanloc_resamp_std);

  % resample err maps across subjects
  Z = results.data.Z_vals;
  orig_size = size(Z);
  matsize = orig_size(1:2);
  Z = reshape(Z,[prod(matsize),orig_size(3)]);
  npoints = size(Z,1);    
  Z_resamp = Z(:,ind_resamp);
  Z_resamp = reshape(Z_resamp,[npoints,results.data.nsubs,parms.niters]);
  Z_resamp_mean = squeeze(mean(Z_resamp,2));
  Z_resamp_sem = squeeze(std(Z_resamp,1,2));
  Z_resamp_lo = prctile(Z_resamp_mean,parms.ci_tail,2);
  Z_resamp_hi = prctile(Z_resamp_mean,100-parms.ci_tail,2);
  results.Z_resamp_mean = reshape(mean(Z_resamp_mean,2),matsize);
  results.Z_resamp_lo = reshape(Z_resamp_lo,matsize);
  results.Z_resamp_hi = reshape(Z_resamp_hi,matsize);
  results.Z_resamp_ci = cat(3,results.Z_resamp_lo,results.Z_resamp_hi);
  results.Z_resamp_sem = reshape(std(Z_resamp_mean,1,2),matsize);

  if ~isempty(parms.prefix2)
    % resample difference in mean err between upper and lower field
    errvals1 = results.data.set1.err_norm_vals;
    errvals2 = results.data.set2.err_norm_vals;
    err1 = mean(errvals1,1);
    err2 = mean(errvals2,1);
    err_diff = err1 - err2;
    err_diff_resamp = reshape(err_diff(ind_resamp),[results.data.nsubs,parms.niters]);
    err_diff_resamp_mean = mean(err_diff_resamp,1);
    results.err_diff_resamp_mean = mean(err_diff_resamp_mean,2);
    results.err_diff_resamp_lo = prctile(err_diff_resamp_mean,parms.ci_tail,2);
    results.err_diff_resamp_hi = prctile(err_diff_resamp_mean,100-parms.ci_tail,2);
    results.err_diff_resamp_ci = ...
      cat(2,results.err_diff_resamp_lo,results.err_diff_resamp_hi);
    results.err_diff_resamp_sem = std(err_diff_resamp_mean,1,2);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_errmaps(parms,results)
  tmp_parms = parms;
  tmp_parms.clim = parms.zlim;
  tmp_parms.outdir = parms.resdir;
  M = results.M;

  % plot mean err map
  if parms.visible_flag, figure(1); end;
  Z = results.Z;
  Z(M~=0) = (Z(M~=0) - results.err_meanloc_resamp_mean)/results.err_meanloc_resamp_std;
  Z(M==0) = nan;
  tmp_parms.title = 'z-transformed err';
  tmp_parms.outstem = sprintf('%s_%s_errmap',parms.outstem,parms.outfix);
  args = mmil_parms2args(tmp_parms,parms.imagesc_tags);
  mmil_imagesc(Z,args{:});

  % plot mean err map, thresholded by confidence intervals
  if parms.visible_flag, figure(2); end;
  Z = results.Z;
  % set values within confidence interval to zero
  Z(Z>results.err_meanloc_resamp_lo & Z<results.err_meanloc_resamp_hi) = nan;
  % subtract mean and normalize by range of 95% confidence interval
  Z(M~=0) = (Z(M~=0) - results.err_meanloc_resamp_mean)/results.err_meanloc_resamp_std;
  Z(M==0) = nan;
  tmp_parms.title = 'thresholded, z-transformed err';
  tmp_parms.outstem = sprintf('%s_%s_errmap_thresh',parms.outstem,parms.outfix);
  args = mmil_parms2args(tmp_parms,parms.imagesc_tags);
  mmil_imagesc(Z,args{:});

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_err_sorted(parms,results)
  % plot mean err, sorted, with confidence intervals
  if parms.visible_flag, figure(3); end;
  tmp_parms = parms;
  % mean error estimate across subjects
  err_subj = results.err_resamp_mean;
  [err_subj,ind_sort] = sort(err_subj);
  % confidence intervals for err, resamped across subjects
  err_subj_lo = results.err_resamp_lo;
  err_subj_hi = results.err_resamp_hi;
  err_subj_lo = err_subj_lo(ind_sort);
  err_subj_hi = err_subj_hi(ind_sort);
  % mean error estimate across locations and subjects
  err_meanloc = results.err_meanloc_resamp_mean*ones(results.data.nlocs,1);
  % confidence intervals for err, resamped across locations and subjects
  err_meanloc_lo = results.err_meanloc_resamp_lo*ones(results.data.nlocs,1);
  err_meanloc_hi = results.err_meanloc_resamp_hi*ones(results.data.nlocs,1);
  % combined for two "conditions"
  err = cat(2,err_subj,err_meanloc);
  err_lo = cat(2,err_subj_lo,err_meanloc_lo);
  err_hi = cat(2,err_subj_hi,err_meanloc_hi);
  tmp_parms.wforms_err = cat(3,err_lo,err_hi);
  tmp_parms.xlabel = 'stimulus location (sorted by error)';
  tmp_parms.ylabel = 'normalized residual error';
  tmp_parms.title = 'normalized residual error';
  tmp_parms.colors = {parms.p1_color,'r'};
  tmp_parms.ylim = parms.err_ylim;
  tmp_parms.outdir = parms.resdir;
  tmp_parms.outstem = sprintf('%s_%s_err',parms.outstem,parms.outfix);
  tmp_parms.time = [1:results.data.nlocs];
  % points to label p1 and p2
  points = [];
  if ~isempty(parms.prefix2)
    [~,~,ind_p1] = intersect(1:results.data.nlocs2,ind_sort);
    [~,~,ind_p2] = intersect(results.data.nlocs2+1:results.data.nlocs,ind_sort);
    points(1).amplitude = err_subj(ind_p1);
    points(1).latency = tmp_parms.time(ind_p1);
    points(1).color = parms.p1_color;
    points(2).amplitude = err_subj(ind_p2);
    points(2).latency = tmp_parms.time(ind_p2);
    points(2).color = parms.p2_color;
  else
    points(1).amplitude = err_subj;
    points(1).latency = tmp_parms.time;
    points(1).color = parms.p1_color;
    points(2).amplitude = [];
    points(2).latency = [];
    points(2).color = 'k';
  end;
  tmp_parms.points = points;
  args = mmil_parms2args(tmp_parms,parms.errvec_tags);
  ts_plot_wforms(err,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_err_unsorted(parms,results)
  % plot mean err, unsorted, with confidence intervals
  if parms.visible_flag, figure(4); end;
  tmp_parms = parms;
  % mean error estimate across subjects
  err_subj = results.err_resamp_mean;
  % confidence intervals for err, resamped across subjects
  err_subj_lo = results.err_resamp_lo;
  err_subj_hi = results.err_resamp_hi;
  % mean error estimate across locations and subjects
  err_meanloc = results.err_meanloc_resamp_mean*ones(results.data.nlocs2,1);
  % confidence intervals for err, resamped across locations and subjects
  err_meanloc_lo = results.err_meanloc_resamp_lo*ones(results.data.nlocs2,1);
  err_meanloc_hi = results.err_meanloc_resamp_hi*ones(results.data.nlocs2,1);
  % combined for three "conditions"
  i_p1 = 1:results.data.nlocs2;
  i_p2 = results.data.nlocs2+1:results.data.nlocs;
  err = cat(2,err_subj(i_p1),err_subj(i_p2),err_meanloc);
  err_lo = cat(2,err_subj_lo(i_p1),err_subj_lo(i_p2),err_meanloc_lo);
  err_hi = cat(2,err_subj_hi(i_p1),err_subj_hi(i_p2),err_meanloc_hi);
  tmp_parms.wforms_err = cat(3,err_lo,err_hi);
  tmp_parms.xlabel = 'stimulus location';
  tmp_parms.ylabel = 'normalized residual error';
  tmp_parms.title = 'normalized residual error';
  tmp_parms.colors = {parms.p1_color,parms.p2_color,'r'};
  tmp_parms.ylim = parms.err_ylim;
  tmp_parms.outdir = parms.resdir;
  tmp_parms.outstem = sprintf('%s_%s_err_unsorted',parms.outstem,parms.outfix);
  tmp_parms.time = [1:results.data.nlocs2];
  % points to label upper and lower
  points = [];
  points(1).amplitude = err_subj(i_p1);
  points(1).latency = tmp_parms.time;
  points(1).color = parms.p1_color;
  points(2).amplitude = err_subj(i_p2);
  points(2).latency = tmp_parms.time;
  points(2).color = parms.p2_color;
  points(3).amplitude = [];
  points(3).latency = [];
  points(3).color = 'k';
  tmp_parms.points = points;
  tmp_parms.condnames = {parms.prefix,parms.prefix2,'overall mean'};
  tmp_parms.condnames = regexprep(tmp_parms.condnames,'_',' ');
  tmp_parms.legend_loc = 'NorthEast';
  args = mmil_parms2args(tmp_parms,parms.errvec_tags);
  ts_plot_wforms(err,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_diff_text(parms,results)
  % display whether there is significant difference between prefix and prefix2
  fname_tmp = sprintf('%s/%s_%s_diff.txt',parms.resdir,parms.outstem,parms.outfix);
  if ~exist(fname_tmp,'file') || parms.forceflag
    fid = fopen(fname_tmp,'wt');
    if fid==-1, error('failed to open %s for writing',fname_tmp); end;
    fprintf(fid,'difference between upper and lower field residual error:\n');
    fprintf(fid,'\n');
    fprintf(fid,'mean: %0.3f\n',results.err_diff_resamp_mean);
    fprintf(fid,'bootstrap confidence interval: [%0.3f,%0.3f]\n',...
      results.err_diff_resamp_ci(1),results.err_diff_resamp_ci(2));
    fprintf(fid,'\n');
    fprintf(fid,'t-stat: %0.3f\n',results.err_diff_tstat);
    fprintf(fid,'p-value: %0.3e\n',results.err_diff_pval);
    fprintf(fid,'t-test confidence interval: [%0.3f,%0.3f]\n',...
      results.err_diff_ttest_ci(1),results.err_diff_ttest_ci(2));
    fprintf(fid,'\n');
    fprintf(fid,'anova1 F: %0.3f\n',results.anova1_table{2,5});
    fprintf(fid,'anova1 p-value: %0.3e\n',results.anova1_table{2,6});
    fprintf(fid,'\n');
    fprintf(fid,'anova2 F: %0.3f\n',results.anova2_table{2,5});
    fprintf(fid,'anova2 p-value: %0.3e\n',results.anova2_table{2,6});
    fprintf(fid,'\n');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

