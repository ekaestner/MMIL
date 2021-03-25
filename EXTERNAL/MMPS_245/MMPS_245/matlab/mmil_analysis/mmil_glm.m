function results = mmil_glm(fname_data,varargin)
%function results = mmil_glm(fname_data,[options])
%
% Purpose: calculate general linear model fit to data
%   contained in csv spreadsheet file or concatenated mgh file
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
%     or number of frames in mgh/mgz fname_data
%    if fname_data is csv file and fname_reg is supplied, they will be
%     joined based on first column, by which they must be sorted
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
%  'write_mgh_flag': [0|1] save results as mgh/mgz files
%     only applies if fname_data is mgh/mgz file
%     {default = 1}
%  'mgh_contrast_labels': cell array of contrast names to save as mgh/mgz files
%     if empty, save all contrasts as mgh/mgz
%     only applies if fname_data is mgh/mgz file
%     {default = []}
%  'statlist': cell array of statistics for which to create output files
%     valid: 'mean','stdv','ssq','r','tstat','Fstat',
%            'pval','log10_pval','ssq_data','ssq_fit','ssq_err'
%            'var_data','var_fit','var_err','exp_var'
%     {default = {'mean','tstat'}}
%  'verbose': [0|1] display warning messages
%     {default = 1}
%  'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   results: structure containing GLM betas, contrasts, etc.
%
% Created:  12/06/12 by Don Hagler
% Prev Mod: 05/04/16 by Don Hagler
% Last Mod: 05/17/17 by Don Hagler
%

%% todo: allow weights input -- extra column?
%%    resize along with data

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: r values are likely incorrect for multiple regression

%% todo: allow fname_mask to reduce number of tests for surf and vol data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin, 1), return; end;

% check input parameters
parms = check_input(fname_data,varargin);

if ~exist(parms.fname_results,'file') || parms.forceflag
  % check input csv files
  parms = check_csv(parms);

  % prepare data matrix
  parms = prep_data(parms);

  % prepare regressor matrix
  parms = prep_reg(parms);

  % exclude rows containing NaNs
  parms = exclude_nans(parms);

  % exclude rows containing outliers
  parms = exclude_outliers(parms);

  % prepare design matrix, contrast struct
  parms = prep_glm(parms);

  % calculate GLM fit, contrasts
  results = calc_glm(parms);

  % store indices of included data rows
  results.ind_included = parms.ind_included;
  
  % store size of non-csv data file
  if ~parms.csv_flag
    results.volsz = parms.volsz;
  end;

  % save results to mat file
  save(parms.fname_results,'results');
else
  load(parms.fname_results);
  if ~parms.csv_flag
    parms.volsz = results.volsz;
  end;
end;

% write text file with info
write_info(parms,results);

% write stat files
if parms.csv_flag
  write_csv_stats(parms,results);
elseif parms.write_mgh_flag
  write_mgh_stats(parms,results);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_data,options)
  parms = mmil_args2parms(options,{...
    'fname_data',fname_data,[],...
  ...
    'outstem','glm',[],...
    'outdir',pwd,[],...
    'fname_reg',[],[],...
    'merge_field',[],[],...
    'data_labels',[],[],...
    'data_label_patterns',[],[],...
    'reg_labels',[],[],...
    'demean_flag',false,[false true],...
    'demean_labels',[],[],...
    'quad_labels',[],[],...
    'cubic_labels',[],[],...
    'interaction_labels',[],[],...
    'outlier_labels',[],[],...
    'outlier_thresh',3,[],...
    'categ_contrast_flag',false,[false true],...
    'write_mgh_flag',true,[false true],...
    'mgh_contrast_labels',[],[],...
    'statlist',{'mean','tstat'},...
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
  ...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'outext','.mgz',{'.mgh','.mgz'},...
    'return_data_flag',false,[false true],...
  ...
    'calc_tags',{'weights','reweight_flag','reweight_fact',...
                'reweight_maxiter','reweight_tol','reweight_leverage_flag',...
                'singularity_flag','contrast_vectors','contrast_names',...
                'nconds'},[],...
  });

  % check input data file
  if ~exist(parms.fname_data,'file')
    error('file %s not found',parms.fname_data);
  end;
  [tmp,tmp,parms.inext] = fileparts(parms.fname_data);
  switch parms.inext
    case '.csv'
      % csv type
      parms.csv_flag = 1;
      parms.hemistr = [];
    case {'.mgh','.mgz'}
      % mgh type
      parms.csv_flag = 0;
      parms.hemistr = get_hemi(parms.fname_data);
    otherwise
      error('invalid input data type: %s',parms.inext);
  end;

  % check input reg file
  if ~isempty(parms.fname_reg) && ~exist(parms.fname_reg,'file')
    error('file %s not found',parms.fname_reg);
  end;

  % check input cell arrays
  parms = check_arrays(parms);

  % check outstem
  if ~mmil_isrelative(parms.outstem)
    parms.outdir = fileparts(parms.outstem);
  else
    parms.outstem = [parms.outdir '/' parms.outstem];
  end;
  mmil_mkdir(parms.outdir);

  % set output file name
  parms.fname_results = sprintf('%s_results%s.mat',...
    parms.outstem,parms.hemistr);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_csv(parms)
  % initialize variables
  parms.ind_num_reg = [];
  parms.ind_cat_reg = [];
  parms.ind_demean_reg = [];
  parms.ind_quad_reg = [];
  parms.ind_cubic_reg = [];
  parms.ind_outlier_reg = [];

  % read csv file(s)
  parms = read_csv(parms);

  % classify column labels and compare against input labels
  if ~isempty(parms.col_labels)
    % determine which columns are numeric and which are categorical
    parms = classify_columns(parms);

    % get column indices for regressors
    parms = check_reg_labels(parms);

    % get column indices for data
    if parms.csv_flag
      parms = check_data_labels(parms);
      % check for overlap between data_labels and reg_labels
      if ~isempty(parms.data_labels) && ~isempty(parms.reg_labels)
        tmp_labels = intersect(parms.data_labels,parms.reg_labels); 
        if ~isempty(tmp_labels)
          error('overlap between data_labels and reg_labels: %s',...
            sprintf('%s ',tmp_labels{:}));
        end;
      end;
    end;
  end;
return;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hemistr = get_hemi(fname_data)
  hemistr = [];
  [tpath,tstem,text] = fileparts(fname_data);
  n = regexp(tstem,'.+-(?<hemi>[lr]h)','names');
  if ~isempty(n)
    hemistr = ['-' n.hemi];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_arrays(parms)
  % check that input cell arrays are cell arrays
  if ~isempty(parms.data_labels)
    if ~iscell(parms.data_labels)
      parms.data_labels = {parms.data_labels};
    end;
    parms.data_labels = mmil_rowvec(parms.data_labels);
  end;
  if ~isempty(parms.data_label_patterns) && ~iscell(parms.data_label_patterns)
    parms.data_label_patterns = {parms.data_label_patterns};
  end;
  if ~isempty(parms.reg_labels) && ~iscell(parms.reg_labels)
    parms.reg_labels = {parms.reg_labels};
  end;
  if ~isempty(parms.demean_labels) && ~iscell(parms.demean_labels)
    parms.demean_labels = {parms.demean_labels};
  end;
  if ~isempty(parms.quad_labels) && ~iscell(parms.quad_labels)
    parms.quad_labels = {parms.quad_labels};
  end;
  if ~isempty(parms.cubic_labels) && ~iscell(parms.cubic_labels)
    parms.cubic_labels = {parms.cubic_labels};
  end;
  if ~isempty(parms.interaction_labels)
    if ~iscell(parms.interaction_labels) || size(parms.interaction_labels,2)~=2
      error('interaction_labels must be a cell matrix with size = [n,2]');
    end;
    parms.num_interactions = size(parms.interaction_labels,1);
  end;
  if ~isempty(parms.outlier_labels) && ~iscell(parms.outlier_labels)
    parms.outlier_labels = {parms.outlier_labels};
  end;
  if ~isempty(parms.mgh_contrast_labels) && ~iscell(parms.mgh_contrast_labels)
    parms.mgh_contrast_labels = {parms.mgh_contrast_labels};
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = read_csv(parms)
  parms.all_vals = [];
  % read fname_data
  if parms.csv_flag
    parms.all_vals = mmil_readtext(parms.fname_data);
    parms.ind_included = [1:size(parms.all_vals,1)];
  end;
  % read fname_reg
  if ~isempty(parms.fname_reg)
    reg_vals = mmil_readtext(parms.fname_reg);
    if parms.csv_flag
      [parms.all_vals,ind_orig] = mmil_merge_cells(parms.all_vals,reg_vals,...
        parms.merge_field);
      parms.ind_included = parms.ind_included(ind_orig);
    else
      parms.all_vals = reg_vals;
    end;
  end;
  % separate column labels from values
  if ~isempty(parms.all_vals) 
    parms.col_labels = parms.all_vals(1,:);
    parms.all_vals = parms.all_vals(2:end,:);
    if parms.csv_flag
      parms.ind_included = parms.ind_included(2:end)-1;
    end;
  else
    parms.col_labels = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = classify_columns(parms)
  [parms.all_vals,parms.ind_num,parms.ind_cat] = ...
    mmil_classify_columns(parms.all_vals,parms.col_labels,parms.verbose);
  parms.num_labels = parms.col_labels(parms.ind_num);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_reg_labels(parms)
  % set indices for different regressors
  if ~isempty(parms.reg_labels)
    % get indices for all regressors excluding mixed columns
    [tmp,parms.ind_reg] = intersect(parms.col_labels,parms.reg_labels);
    % check for invalid (missing) reg column labels
    if length(parms.ind_reg) < length(parms.reg_labels)
      inval_labels = setdiff(parms.reg_labels,parms.col_labels);
      if ~isempty(inval_labels)
        fprintf('%s: WARNING: ignoring invalid reg_labels: %s\n',...
          mfilename,sprintf('"%s" ',inval_labels{:}));    
      end;
      parms.reg_labels = parms.col_labels(parms.ind_reg);
    end;
    % get indices for numeric regressors
    parms.ind_num_reg = intersect(parms.ind_reg,parms.ind_num);
    parms.num_reg_labels = parms.col_labels(parms.ind_num_reg);
    % get indices for categorical regressors
    parms.ind_cat_reg = intersect(parms.ind_reg,parms.ind_cat);
    parms.cat_reg_labels = parms.col_labels(parms.ind_cat_reg);
    % check labels for regressors to be de-meaned
    if parms.demean_flag
      if isempty(parms.demean_labels)
        parms.demean_labels = parms.num_reg_labels;
      end;
      % remove any members that are not in num_reg_labels
      parms = exclude_inval_labels(parms,'demean');
    end;
    % check labels for regressors to be fitted with quadratic function
    if ~isempty(parms.quad_labels)
      % remove any members that are also be in cubic_labels
      if ~isempty(parms.cubic_labels)
        parms.quad_labels = setdiff(parms.quad_labels,parms.cubic_labels);
      end;
      % remove any members that are not in num_reg_labels
      parms = exclude_inval_labels(parms,'quad');
    end;
    % check labels for regressors to be fitted with cubic function
    if ~isempty(parms.cubic_labels)
      % remove any members that are not in num_reg_labels
      parms = exclude_inval_labels(parms,'cubic');
    end;
    % check labels for regressors to be filtered for outliers
    if ~isempty(parms.outlier_labels)
      % remove any members that are not in num_reg_labels
      parms = exclude_inval_labels(parms,'outlier',1);
    end;    
    % check labels for regressors to be used for interaction terms
    if ~isempty(parms.interaction_labels)
      % error if any are not in num_reg_labels
      for i=1:numel(parms.interaction_labels)
        if ~ismember(parms.interaction_labels{i},parms.num_reg_labels)
          error('interaction label %s not a numeric regressor',...
            parms.interaction_labels{i});
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = exclude_inval_labels(parms,label_type,any_flag)
  if ~exist('any_flag','var') || isempty(any_flag), any_flag = 0; end;
  % note: if any_flag = 1, allow label to be any numeric label
  %   otherwise, must be a regressor
  tag = [label_type '_labels'];
  % define set of valid labels
  if any_flag
    valid_labels = parms.num_labels;
  else
    valid_labels = parms.num_reg_labels;
  end;
  % remove any members that are not in valid_labels
  inval_labels = setdiff(parms.(tag),valid_labels);
  if ~isempty(inval_labels)
    fprintf('%s: WARNING: ignoring invalid %s: %s\n',...
      mfilename,tag,sprintf('"%s" ',inval_labels{:}));    
    parms.(tag) = ...
      intersect(parms.(tag),valid_labels);        
  end;
  tag_ind = ['ind_' label_type '_reg'];
  [tmp,parms.(tag_ind)] = ...
    intersect(parms.col_labels,parms.(tag));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_data_labels(parms)
  if ~isempty(parms.data_label_patterns)
    % add data_labels based on matches to data_label_patterns
    parms = match_label_patterns(parms);  
  end;
  if isempty(parms.data_labels)
    % set data_labels if empty
    parms.ind_data = setdiff(parms.ind_num,parms.ind_num_reg);
    parms.data_labels = parms.col_labels(parms.ind_data);
  else
    % get indices of data columns
    [tmp,ind_cols,ind_data] = intersect(parms.col_labels,parms.data_labels);
    % undo alphabetical reordering done by intersect
    [ind_data_sorted,ind_sort] = sort(ind_data);
    parms.ind_data = ind_cols(ind_sort);
    % check for invalid (missing) data column labels
    if length(parms.ind_data) < length(parms.data_labels)
      inval_labels = setdiff(parms.data_labels,parms.col_labels);
      fprintf('%s: WARNING: ignoring invalid data_labels: %s\n',...
        mfilename,sprintf('"%s" ',inval_labels{:}));    
      parms.data_labels = parms.col_labels(parms.ind_data);
    end;
    % check for nonnumeric data column labels
    ind_nonnum = find(~ismember(parms.ind_data,parms.ind_num));
    if ~isempty(ind_nonnum)
      fprintf('%s: WARNING: excluding non-numeric data columns: %s\n',...
        mfilename,sprintf('"%s" ',parms.col_labels{parms.ind_data(ind_nonnum)}));
      [tmp,ind_data] = intersect(parms.ind_data,parms.ind_num);
      parms.ind_data = parms.ind_data(sort(ind_data));
      parms.data_labels = parms.col_labels(parms.ind_data);
    end;
  end;
  if isempty(parms.ind_data)
    error('no valid data columns specified');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = match_label_patterns(parms)
  for i=1:length(parms.data_label_patterns)
    pat = parms.data_label_patterns{i};
    res = regexp(parms.col_labels,pat);
    ind = find(~cellfun(@isempty,res));
    if isempty(ind)
      fprintf('%s: WARNING: no match for data label pattern "%s"\n',...
        mfilename,pat);
    else
      parms.data_labels = cat(2,parms.data_labels,parms.col_labels(ind));
    end;
  end;
  [tmp,ind_uniq] = unique(parms.data_labels,'first');
  ind_uniq = sort(ind_uniq);
  parms.data_labels = parms.data_labels(ind_uniq);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_data(parms)
  if parms.csv_flag
    parms.data = cell2mat(parms.all_vals(:,parms.ind_data));
  else
    parms.data = fs_load_mgh(parms.fname_data);
    parms.volsz = size(parms.data);
    if length(parms.volsz)<4
      error('multiple frames required');
    end;
    parms.ind_included = [1:parms.volsz(4)];
    %% todo: select voxels/vertices from mask
    parms.data = reshape(parms.data,[prod(parms.volsz(1:3)),parms.volsz(4)])';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_reg(parms)
  parms.regressors = [];
  parms.regnames = [];
  parms.orig_regressors = [];
  parms.orig_regnames = [];
  parms.groupregs = [];
  nvals = size(parms.all_vals,1);

  % create columns for numerical variables
  if ~isempty(parms.ind_num_reg)
    for j=1:length(parms.ind_num_reg)
      ind = parms.ind_num_reg(j);
      vals = cell2mat(parms.all_vals(:,ind));
      regname = parms.col_labels{ind};
      parms.orig_regressors = cat(2,parms.orig_regressors,vals);
      if isempty(parms.orig_regnames), parms.orig_regnames = {}; end;
      parms.orig_regnames = cat(2,parms.orig_regnames,regname);
      if parms.demean_flag && ismember(ind,parms.ind_demean_reg) &&...
          ~all(isnan(vals))
        % subtract mean and normalize by standard deviation
        mean_val = mean(vals(~isnan(vals)));
        std_val = std(vals(~isnan(vals)));
        vals = (vals - mean_val)/std_val;
      end;
      % expand quadratic and cubic terms
      if ismember(ind,parms.ind_cubic_reg)
        np = 3;
      elseif ismember(ind,parms.ind_quad_reg)
        np = 2;
      else
        np = 1;
      end;
      for p=1:np
        if p>1
          tmp_vals = vals.^p;
          tmp_regname = sprintf('%s^%d',regname,p);
        else
          tmp_vals = vals;
          tmp_regname = regname;
        end;
        parms.regressors = cat(2,parms.regressors,tmp_vals);
        if isempty(parms.regnames), parms.regnames = {}; end;
        parms.regnames = cat(2,parms.regnames,tmp_regname);
      end;
    end;
  end;

  % create columns for categorical variables
  if ~isempty(parms.ind_cat_reg)
    for j=1:length(parms.ind_cat_reg)
      ind = parms.ind_cat_reg(j);
      vals = parms.all_vals(:,ind);
      regname = parms.col_labels{ind};
      % expand categorical regressors
      uniq_vals = unique(vals(~strcmpi(vals,'nan')));
      uniq_vals = regexprep(uniq_vals,'\s+','_');
      num_uniq = length(uniq_vals);
      tmp_groupregs = [];
      if num_uniq>1
        for k=1:num_uniq
          tmp_val = uniq_vals{k};
          tmp_regname = [regname ':' tmp_val];
          tmp_reg = zeros(nvals,1);
          tmp_reg(strcmp(vals,tmp_val)) = 1;
          tmp_reg(strcmpi(vals,'nan')) = NaN;
          parms.regressors = cat(2,parms.regressors,tmp_reg);
          parms.regnames = cat(2,parms.regnames,tmp_regname);
          tmp_groupregs = cat(2,tmp_groupregs,length(parms.regnames));
        end;
        parms.groupregs{end+1} = tmp_groupregs;
      end;
    end;
  end;

  % create interaction terms
  if ~isempty(parms.interaction_labels)
    for i=1:parms.num_interactions
      tmp_label1 = parms.interaction_labels{i,1};
      tmp_label2 = parms.interaction_labels{i,2};
      tmp_regname = sprintf('%s_X_%s',tmp_label1,tmp_label2);
      ind1 = find(strcmp(tmp_label1,parms.regnames));
      ind2 = find(strcmp(tmp_label2,parms.regnames));
      tmp_vals1 = parms.regressors(:,ind1);
      tmp_vals2 = parms.regressors(:,ind2);
      tmp_vals = tmp_vals1 .* tmp_vals2;
      parms.regressors = cat(2,parms.regressors,tmp_vals);
      parms.regnames = cat(2,parms.regnames,tmp_regname);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = exclude_nans(parms)
  % exclude subjects with NaNs
  nrows = size(parms.regressors,1);
  tmp = sum(cat(2,parms.data,parms.regressors),2);
  ind_valid = find(~isnan(tmp));
  nvalid = length(ind_valid);
  if nvalid==0
    error('no valid rows without NaNs');
  elseif nvalid<nrows
    fprintf('%s: WARNING: excluding %d rows containing NaNs...\n',...
      mfilename,nrows - nvalid);
    parms.data = parms.data(ind_valid,:);
    parms.all_vals = parms.all_vals(ind_valid,:);
    if ~isempty(parms.regressors)
      parms.regressors = parms.regressors(ind_valid,:);
    end;
    if ~isempty(parms.orig_regressors)
      parms.orig_regressors = parms.orig_regressors(ind_valid,:);
    end;
    parms.ind_included = parms.ind_included(ind_valid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = exclude_outliers(parms)
  if isempty(parms.ind_outlier_reg), return; end;
  % exclude subjects with outliers
  nrows = size(parms.regressors,1);
  valid_flags = ones(nrows,1);
  for i=1:length(parms.ind_outlier_reg)
    ind_col = parms.ind_outlier_reg(i);
    vals = cell2mat(parms.all_vals(:,ind_col));
    mean_val = mean(vals(~isnan(vals)));
    std_val = std(vals(~isnan(vals)));
    vals(isnan(vals)) = Inf;
    ind_invalid = find((abs(vals - mean_val)/std_val) > parms.outlier_thresh);
    valid_flags(ind_invalid) = 0;
  end;
  ind_valid = find(valid_flags);
  nvalid = length(ind_valid);
  if nvalid==0
    error('no valid rows without outliers');
  elseif nvalid<nrows
    fprintf('%s: WARNING: excluding %d rows containing outliers...\n',...
      mfilename,nrows - nvalid);
    parms.data = parms.data(ind_valid,:);
    parms.all_vals = parms.all_vals(ind_valid,:);
    if ~isempty(parms.regressors)
      parms.regressors = parms.regressors(ind_valid,:);
    end;
    if ~isempty(parms.orig_regressors)
      parms.orig_regressors = parms.orig_regressors(ind_valid,:);
    end;
    parms.ind_included = parms.ind_included(ind_valid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_glm(parms)
  nvals = size(parms.all_vals,1);
  if isempty(parms.groupregs)
    baseflag = 1;
  else
    baseflag = 0;
  end;
  parms.X = mmil_glm_design(parms.regressors,'nsubs',nvals,...
    'baseflag',baseflag);
  contrasts = mmil_glm_contrasts(parms.regressors,...
    'group_contrast_flag',parms.categ_contrast_flag,...
    'regnames',parms.regnames,'groupregs',parms.groupregs,...
    'baseflag',baseflag);
  parms.contrast_vectors = contrasts.contrast_vectors;
  parms.contrast_names = contrasts.contrast_names;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_glm(parms)
  args = mmil_parms2args(parms,parms.calc_tags);
  results = mmil_glm_calc(parms.X,parms.data,args{:});
  % save additional information
  if ~isempty(parms.data_labels)
    results.data_labels = parms.data_labels;
  end;
  if parms.return_data_flag
    results.data = parms.data;
  end;
  results.orig_regressors = parms.orig_regressors;
  results.orig_regnames = parms.orig_regnames;
  results.regressors = parms.regressors;
  results.regnames = parms.regnames;
  results.groupregs = parms.groupregs;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_info(parms,results)
  % write info to text file
  fname = sprintf('%s_info%s.txt',parms.outstem,parms.hemistr);
  if ~exist(fname,'file') || parms.forceflag
    fid = fopen(fname,'wt');
    if fid==-1
      error('failed to open %s for writing');
    end;
    if iscell(results.regnames)
      nregs = length(results.regnames);
    else
     nregs = 1;
    end;
    fprintf(fid,'N = %d\n',results.nsubs);
    fprintf(fid,'# of regressors = %d\n',nregs);
    fprintf(fid,'# of parameters = %d\n',results.nparams);
    fprintf(fid,'dof = %d\n',results.dof);
    fprintf(fid,'# of independent tests = %d\n',results.nvals);
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_csv_stats(parms,results)
  % save stats to csv files
  for s=1:length(parms.statlist)
    statname = parms.statlist{s};
    fname = [parms.outstem '_' statname '.csv'];
    if ~exist(fname,'file') || parms.forceflag
      row_labels = results.data_labels;
      firstcol_label = 'dependent_variable';
      switch statname
        case {'ssq_data','ssq_fit','ssq_err',...
              'var_data','var_fit','var_err','exp_var'}
          col_labels = {statname};
          nr = length(row_labels);
          nc = length(col_labels);
          data = zeros(nr,nc);
          data(:,1) = results.(statname);
        otherwise
          col_labels = {results.contrasts.name};
          nr = length(row_labels);
          nc = length(col_labels);
          data = zeros(nr,nc);
          for c=1:nc
            data(:,c) = results.contrasts(c).(statname);
          end;
      end;
      mmil_write_csv(fname,data,...
        'row_labels',row_labels,...
        'col_labels',col_labels,...
        'firstcol_label',firstcol_label);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_mgh_stats(parms,results)
  % set string for hemisphere
  % save stats to mgh/mgz files
  for c=1:length(results.contrasts)
    contname = results.contrasts(c).name;
    % only save selected
    if ~isempty(parms.mgh_contrast_labels) &&...
       ~ismember(contname,parms.mgh_contrast_labels)
      continue;
    end;
    for s=1:length(parms.statlist)
      statname = parms.statlist{s};
      fname = sprintf('%s_%s_%s%s%s',...
        parms.outstem,contname,statname,parms.hemistr,parms.outext);
      if ~exist(fname,'file') || parms.forceflag
        vals = mmil_colvec(results.contrasts(c).(statname));
        %% todo: if fname_mask supplied, restore to full volume/surface
        vals = reshape(vals,parms.volsz(1:3));
        fs_save_mgh(vals,fname);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

