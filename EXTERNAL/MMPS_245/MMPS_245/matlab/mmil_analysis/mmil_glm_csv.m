function results = mmil_glm_csv(fname_data,varargin)
%function results = mmil_glm_csv(fname_data,[options])
%
% Purpose: calculate general linear model fit to data
%   contained in csv spreadsheet file and create scatter plots
%
% Usage:
%  results = mmil_glm(fname,'key1', value1,...);
%
% Required Input:
%  fname_data: name of input data file
%     may be either csv spreadsheet (containing data and regressors) or
%       mgh/mgz file (containing concatenated surface data)
%
% Optional Input:
%  'outstem': output file stem
%    if full path is given, will override outdir
%    {default = 'glm'}
%  'outdir': output directory
%    {default = pwd}
%  'fname_reg': name of csv spreadsheet file containing covariate regressors
%    number of rows (minus header) must match number of rows in csv fname_data
%    if fname_reg is supplied, they will be joined based on first column
%       or merge field, by which they must be sorted
%    {default = []}
%  'merge_field': column header used to merge fname_data and fname_reg
%     if fname_data is csv file and fname_reg is also supplied
%     if merge_field not supplied, will use first column header in common
%    {default = []}
%  'data_labels': cell array of column headers for dependent variables
%    only relevant if fname_data is a csv file
%    if empty, will use all numeric columns except those in reg_labels
%    {default = []}
%  'data_label_patterns': cell array of partial column headers
%    only relevant if fname_data is a csv file
%    if supplied, will use regexp to match to column headers
%    {default = []}
%  'reg_labels': cell array of column headers for covariate regressor values
%    {default = []}
%  'demean_flag': [0|1] for numeric regressors, subtract mean
%    and normalize by standard deviation
%    {default = 0}
%  'demean_labels': cell array of regressors to de-mean
%    must be numeric and included in reg_labels
%    if empty, and demean_flag = 1, de-mean all numeric regressors
%    {default = []}
%  'quad_labels': cell array of regressors to be modeled
%     with linear and quadratic terms
%     must be numeric and included in reg_labels
%    {default = []}
%  'cubic_labels': cell array of regressors to be modeled
%     with linear, quadratic, and cubic terms
%     must be numeric and included in reg_labels
%    {default = []}
%  'interaction_labels': cell matrix of regressor names (size = [n,2])
%     to be used as interaction terms
%     must be numeric and included in reg_labels
%    {default = []}
%  'outlier_labels': cell array of regressor or data labels
%     for which to exclude outliers
%       outliers defined as absolute difference from mean greater than
%       outlier_thresh * stdev
%    {default = []}
%  'outlier_thresh': threshold applied for outlier regressors or data points
%     multiple of standard deviation
%    {default = 3}
%  'categ_contrast_flag': [0|1] include GLM tests for direct contrasts
%     between values of categorical variables
%     {default = 0}
%  'statlist': cell array of statistics for which to create output files
%     valid: 'mean','stdv','ssq','r','tstat','Fstat',
%            'pval','log10_pval','ssq_data','ssq_fit','ssq_err'
%            'var_data','var_fit','var_err','exp_var'
%     {default = {'mean','tstat'}}
%  'plotflag': [0|1] create scatter plots
%    {default = 1}
%  'plot_reg_labels': cell array of reg_labels for which to create plots
%    if empty, create plots for all numeric reg_labels
%    {default = []}
%  'rmbaseflag': [0|1] remove baseline values
%    {default = 0}
%  'residflag': [0|1] plot residualized values, removing all other regressors
%    {default = 1}
%  'title_stem': stem of plot title
%    {default = 'GLM results'}
%  'ylabel': y-axis label
%    {default = 'data'}
%  'yrange': y-axis range limits (vector of two numbers)
%    if empty, auto-scale plots
%    {default = []}
%  'legend_flag': [0|1] display legend
%    {default = 1}
%  'legloc': location of legend
%    {default = 'EastOutside'}
%  'pval_flag': [0|1] display minimum p-value for each plot
%    {default = 1}
%  'text_xpos': x-position of p-value text label relative to left side
%    {default = 0.6}
%  'text_ypos': y-position of p-value text label relatiev to bottom
%    {default = 0.05}
%  'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   results: structure containing GLM betas, contrasts, etc.
%
% Created:  12/17/12 by Don Hagler
% Last Mod: 05/26/15 by Don Hagler
%

%% todo: allow weights input -- extra column?

%% todo: put in header once all issues are resolved (see mmil_glm_calc)
%  'reweight_flag' [0|1]: perform iteratively reweighted least squares (IRLS)
%    {default = 0}
%  'reweight_fact': reweighting factor (multiple of median absolute deviation)
%    {default = 4.7}
%  'reweight_maxiter': maximum number of reweighting iterations
%    {default = 100}
%  'reweight_tol': tolerance value for change in weights (stop iterations)
%    {default = 10^-7}
%  'reweight_leverage_flag': [0|1] for IRLS, adjust residuals using leverage
%    {default = 1}

%% todo: scatter plots with each group as different color and regression line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(fname_data,varargin);

args = mmil_parms2args(parms,parms.glm_tags);
results = mmil_glm(fname_data,args{:});

if parms.plotflag
  regnames = intersect(results.orig_regnames,parms.plot_reg_labels);
  for r=1:length(regnames)
    regname = regnames{r};
    ind = find(strcmp(regname,results.orig_regnames));
    x = results.orig_regressors(:,ind);
    plot_data_fit(parms,results,regname,x);
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_data,options)
  parms = mmil_args2parms(options,{,...
    'outstem','glm',[],...
    'outdir',pwd,[],...
    'fname_reg',[],[],...
    'merge_field','SubjID',[],...
    'data_labels',[],[],...
    'data_label_patterns',[],[],...
    'reg_labels',[],[],...
    'quad_labels',[],[],...
    'cubic_labels',[],[],...
    'interaction_labels',[],[],...
    'outlier_labels',[],[],...
    'outlier_thresh',3,[],...
    'demean_flag',true,[false true],...
    'demean_labels',[],[],...
    'return_data_flag',true,[false true],...
    'categ_contrast_flag',false,[false true],...
    'statlist',{'mean','tstat','pval','Fstat','r'},...
               {'mean','stdv','ssq','r','tstat','Fstat',...
                'pval','log10_pval',...
                'ssq_data','ssq_fit','ssq_err',...
                'var_data','var_fit','var_err','exp_var'},...
  ... % IRLS
    'reweight_flag',false,[false true],...
    'reweight_fact',4.7,[0.1,100],...
    'reweight_maxiter',100,[1,1000],...
    'reweight_tol',1e-7,[1e-100,1e-2],...
    'reweight_leverage_flag',true,[false true],...
  ... % plotting
    'plotflag',true,[false true],...
    'plot_reg_labels',[],[],...
    'rmbaseflag',false,[false true],...
    'residflag',true,[false true],...
    'title_stem','GLM results',[],...
    'ylabel','data',[],...
    'yrange',[],[],...
    'legend_flag',true,[false true],...
    'legloc','EastOutside',[],...
    'pval_flag',true,[false true],...
    'text_xpos',0.6,[0 1],...
    'text_ypos',0.05,[0 1],...
    'forceflag',false,[false true],...
  ...
    'glm_tags',{'outstem','outdir','fname_reg','merge_field',...
      'data_labels','data_label_patterns','reg_labels','demean_flag',...
      'demean_labels','quad_labels','cubic_labels','interaction_labels',...
      'outlier_labels','outlier_thresh','categ_contrast_flag',...
      'statlist','weights','reweight_flag','reweight_fact',...
      'reweight_maxiter','reweight_tol','reweight_leverage_flag',...      
      'forceflag','outext','return_data_flag'},[],...
  });

  % check input data file
  if ~exist(fname_data,'file')
    error('file %s not found',fname_data);
  end;
  [tmp,tmp,parms.inext] = fileparts(fname_data);
  if ~strcmp(parms.inext,'.csv')
    error('invalid input data type: %s',parms.inext);
  end;

  if isempty(parms.plot_reg_labels)
    parms.plot_reg_labels = parms.reg_labels;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_data_fit(parms,results,regname,x)
  % set labels and title string stem
  regname_label = regexprep(regname,'_',' ');
  title_stem = [parms.title_stem ' vs ' regname_label];
  legstr = regexprep(results.data_labels,'_',' ');
  % get parameters, data, and fit from results
  betas = results.betas;
  data = results.data;
  X = results.X;
  fit = X*betas;
  % index for linear regressor
  ind_lin = find(strcmp(results.regnames,regname));
  % p-values for linear regressor
  pvals = results.contrasts(1+ind_lin).pval;
  % indices for linear and nonlinear regressors
  ind_all = ind_lin;
  if ismember(regname,parms.quad_labels)
    ind_tmp = find(strcmp(results.regnames,[regname '^2']));
    ind_all = [ind_all,ind_tmp];
  end;
  if ismember(regname,parms.cubic_labels)
    ind_tmp = find(strcmp(results.regnames,[regname '^3']));
    ind_all = [ind_all,ind_tmp];
  end;
  if isempty(results.groupregs)
    % adjust for extra baseline column
    ind_all = 1+ind_all;
  end;
  % calculate baseline values for each variable
  baseline = results.contrasts(1).vec*betas;
  % calculate residualized data for plotting
  tmp_data = data;
  if parms.residflag
    % calculate fit excluding this set of regressors
    tmp_betas = betas;
    tmp_betas(ind_all,:) = 0;
    tmp_fit_excl = X*tmp_betas;
    tmp_data = data - tmp_fit_excl;
    if parms.rmbaseflag
      % remove baseline differences between data variables
      title_str = [title_stem ' with baseline and other effects removed'];
    else
      % keep baseline differences between data variables
      title_str = [title_stem ' with other effects removed'];
      tmp_data = bsxfun(@plus,tmp_data,baseline);      
    end;
  else
    if parms.rmbaseflag
      % remove baseline differences between variables
      title_str = [title_stem ' with baseline removed'];
      tmp_data = bsxfun(@minus,tmp_data,baseline);      
    else
      % keep baseline differences between variables
      title_str = title_stem;
    end;
  end;

  % calculate fit for plotting
  tmp_betas = betas(ind_all,:);
  tmp_X = X(:,ind_all);
  tmp_fit = tmp_X*tmp_betas;
  if ~parms.rmbaseflag
    tmp_fit = bsxfun(@plus,tmp_fit,baseline);
  end;

  xrange = [min(x),max(x)];

  % plot data and fit
  outstem = [parms.outstem '_vs_' regname];
  if parms.residflag
    outstem = [outstem '_resid'];
  end;
  if parms.rmbaseflag
    outstem = [outstem '_rmbase'];
  end;
  plot_scatter(x,tmp_data,tmp_fit,...
    title_str,outstem,...
    xrange,parms.yrange,...
    regname_label,parms.ylabel,legstr,pvals,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_scatter(x,y,yfit,title_str,outstem,xrange,yrange,xlab,ylab,legstr,pvals,parms);
  fname_out = [parms.outdir '/' regexprep(outstem,' ','_') '.tif'];
  if ~exist(fname_out,'file') || parms.forceflag
    figure; clf;
    set(gcf,'visible','off');
    hold on;
    plot(x,y,'*');
    [x_sort,ind_sort] = sort(x);
    yfit_sort = yfit(ind_sort,:);    
    plot(x_sort,yfit_sort,'-','LineWidth',2);
    xlabel(xlab);
    ylabel(ylab);
    set(gca,'xlim',xrange);
    if ~isempty(yrange)
      set(gca,'ylim',yrange);
    else
      yrange = get(gca,'ylim');
    end;
    title(regexprep(title_str,'_',' '));
    if parms.legend_flag
      legend(legstr,'Location',parms.legloc);
    end;
    if parms.pval_flag
      min_pval = min(pvals);
      if min_pval < 0.01
        col = [1 0 0];
      elseif min_pval < 0.05
        col = [0 1 0];
      else
        col = [0 0 1];
      end;
      tx = xrange(1) + parms.text_xpos * range(xrange);
      ty = yrange(1) + parms.text_ypos * range(yrange);
      text(tx,ty,sprintf('min p = %0.5f',min_pval),'Color',col);
    end;
    print(gcf,'-dtiff',fname_out,'-r 300');
    close(gcf);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

