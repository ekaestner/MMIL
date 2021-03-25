function results = mmil_twin_analysis(fstem_data,fname_info,varargin)
%function results = mmil_twin_analysis(fstem_data,fname_info,[options])
%
% Purpose: use MX to create cortical surface maps of genetic relationships
%
% Usage:
%  results = mmil_twin_analysis(fname,'key1', value1,...);
%
% Required Input:
%  fstem_data: file stem of input data files (left and right hemispheres)
%    full name will be <fstem_data>-<hemi><inext>
%    alternative: input as matrix with [nvox,nsubj] instead of string
%  fname_info: name of csv spreadsheet file containing subject info
%    including twin pairs and covariates
%    number of rows (minus header) must match number of frames in data files
%
% Optional MX parameters:
%   'mxbindir': directory containing MX binary
%     may be empty if MX is in path or if mxbin is supplied as full path
%     {default = []}
%   'mxscriptdir': directory containing MX scripts
%     ignored if mxscript is supplied as full path
%     if empty, will use $MMPS_DIR/mx
%     {default = []}
%   'mxbin': name of MX executable binary
%     {default = 'mxt167'}
%   'mxscript': name of MX script
%     {default = 'fast_univariate_ACE.mx'}
%   'mxtmpdir': directory in which to create temporary output files
%     {default = '/tmp/mxtmp'}
%
% Optional input data parameters:
%  'indir': input directory
%    ignored if fstem_data or fname_info contain absolute paths
%    {default = pwd}
%  'inext': input file extension
%    {default = '.mgz'}
%  'fname_reg': name of csv spreadsheet file containing additional regressors
%    number of rows (minus header) must match number of rows in fname_info
%    {default = []}
%  'merge_field': column header used to merge fname_info and fname_reg
%     if merge_field not supplied, will use first column header in common
%    {default = []}
%  'fnames_label': cell array of label files (one for each hemisphere)
%    if empty, will use label_name
%    {default = []}
%  'label_name': name of label file used to limit calculations
%    use 'none' to include all vertices in calculations
%    {default = 'cortex'}
%  'subjname': name of freesurfer subject with label files
%    {default = 'fsaverage'}
%  'subjdir': full path of directory containing subjname
%    {default = $FREESURFER_HOME/subjects}
%  'init_data_flag': [0|1] load data from mgz and save to mat file, then quit
%    {default = 0}
%
% Optional preprocessing parameters:
%  'scalefact': scaling factor applied to input data
%    {default = 1}
%  'global_flag': [0|1|2] remove global effects from data
%    0: do not remove global effects
%    1: subtract mean value of all vertices (e.g. for thickness)
%    2: divide by sum of all vertices (e.g. for area)
%    {default = 1}
%
% Optional regressor parameters:
%  'seed_flag': [0|1] calculate bivariate correlations between seed and all vertices
%    if true, mxscript for bivariate analysis must be specified
%    {default = 0}
%  'seed_vert': seed point vertex number
%    {default = 1}
%  'seed_hemi': seed point hemisphere ('lh' or 'rh')
%    {default = 'lh'}
%  'bivar_label': column header for bivariate regressor values
%    if supplied, mxscript for bivariate analysis must be specified
%    ignored if seed_flag = 1
%    {default = []}
%  'pair_label': column header for pair IDs
%    each pair of twin subjects must share the same pair ID
%    {default = 'pair'}
%  'twin_label': column header for twin IDs
%    if supplied, each twin should have one of two different values
%    e.g. 'A' or 'B', 1 or 2
%    {default = []}
%  'zyg_label': column header for zygosity values
%    with 1 for MZ and 2 for DZ
%    {default = 'zyg'}
%  'resid_labels': cell array of column headers for covariate regressor values
%    to use for pre-residualization of data
%    {default = []}
%  'quad_labels': cell array of regressors to be modeled
%     with linear and quadratic terms
%     must be numeric and included in resid_labels
%    {default = []}
%  'cubic_labels': cell array of regressors to be modeled
%     with linear, quadratic, and cubic terms
%     must be numeric and included in resid_labels
%    {default = []}
%  'outlier_labels': cell array of regressor or data labels
%     for which to exclude outliers
%       outliers defined as absolute difference from mean greater than
%       outlier_thresh * stdev
%    {default = []}
%  'outlier_thresh': threshold applied for outlier regressors or data points
%     multiple of standard deviation
%    {default = 3}
%
% Optional output data parameters:
%  'outdir': output directory
%    {default = pwd}
%  'outstem': output file stem
%    if full path is given, will override outdir
%    {default = 'twin'}
%  'outext': output file extension
%    {default = '.mgz'}
%  'pheno_flag': [0|1] calculate phenotypic correlations
%    {default = 1}
%  'pheno_xcorr_flag':  [0|1] calculate phenotypic cross-correlations between vertices
%    {default = 0}
%  'mx_flag': [0|1] calculate genetic correlations
%    {default = 1}
%  'write_mgh_flag': [0|1] save results as mgh/mgz files
%    {default = 1}
%  'verbose': [0|1] display frequent status messages
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Output:
%   results: structure containing results
%
% Created:  02/21/14 by Don Hagler
% Prev Mod: 09/26/17 by Don Hagler
% Last Mod: 10/17/17 by Don Hagler
%

%% todo: render images?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin, 1), return; end;

% check input parameters
parms = check_input(fstem_data,fname_info,varargin);

if ~exist(parms.fname_results,'file') || parms.forceflag
  % load info csv files
  parms = load_info(parms);

  % load data matrix
  if ~isempty(parms.fstem_data)
    parms = load_data(parms);
  end;
  if parms.init_data_flag, return; end;
  
  % check input seed index
  parms = check_seed(parms);

  % exclude rows containing NaNs
  parms = exclude_nans(parms);

  % exclude rows containing outliers
  parms = exclude_outliers(parms);

  % exclude subjects without a twin pair or zygosity info
  parms = exclude_nontwins(parms);

  % demean numeric residualization regressors and bivariate regressor
  parms = demean_regs(parms);

  % scale data
  parms = scale_data(parms);

  % remove global effects from data
  parms = remove_global_data(parms);

  % residualize data for resid_labels
  parms = residualize_data(parms);

  % z-transform data
  parms = zscore_data(parms);

  % initialize results struct
  results = init_results(parms);

  % calculate phenotypic correlations
  if parms.pheno_flag
    results = calc_pheno_corr(parms,results);
  end;

  % calculate phenotypic cross-correlation (between vertices)
  if parms.pheno_xcorr_flag
    results = calc_pheno_xcorr(parms,results);
  end;
  
  % run MX calculations
  if parms.mx_flag
    results = calc_MX(parms,results);
  end;

  % save results to mat file
  if parms.verbose
    fprintf('%s: saving results to %s...\n',mfilename,parms.fname_results);
  end;
  save(parms.fname_results,'results','-v7.3');
else
  load(parms.fname_results);
end;

% save output
if parms.write_mgh_flag
  write_mgh_results(parms,results);
end;

% remove mxtmpdir
if parms.cleanup_flag
  cleanup(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fstem_data,fname_info,options)
  parms = mmil_args2parms(options,{...
    'fstem_data',fstem_data,[],...
    'fname_info',fname_info,[],...
  ... % MX
    'mxbindir',[],[],...
    'mxscriptdir',[],[],...
    'mxbin','mxt167',[],...
    'mxscript','fast_univariate_ACE.mx',[],...
    'mxtmpdir','/tmp/mxtmp',[],...
  ... % input
    'indir',pwd,[],...
    'inext','.mgz',{'.mgh','.mgz'},...
    'fname_reg',[],[],...
    'merge_field',[],[],...
    'fnames_label',[],[],...
    'label_name','cortex',[],...
    'subjdir',[],[],...
    'subjname','fsaverage',[],...
    'init_data_flag',false,[false true],...
  ... % preprocessing
    'scalefact',1,[1e-10,1e10],...
    'global_flag',1,[0:2],...
  ... % regressors
    'seed_flag',false,[false true],...
    'seed_vert',1,[1,1e10],...
    'seed_hemi','lh',{'lh','rh'},...
    'bivar_label',[],[],...
    'pair_label','pair',[],...
    'twin_label',[],[],...
    'zyg_label','zyg',[],...
    'resid_labels',[],[],...
    'quad_labels',[],[],...
    'cubic_labels',[],[],...
    'outlier_labels',[],[],...
    'outlier_thresh',3,[],...
  ... % output
    'outdir',pwd,[],...
    'outstem','glm',[],...
    'outext','.mgz',{'.mgh','.mgz'},...
    'pheno_flag',true,[false true],...
    'pheno_xcorr_flag',false,[false true],...
    'mx_flag',true,[false true],...
    'write_mgh_flag',true,[false true],...
    'cleanup_flag',true,[false true],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % undocumented
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'maxnvals',20000,[1,Inf],...
    'dump_mxout_flag',false,[false true],...
    'bivar_pheno_outlist',{'rMZ','rDZ','rP','rP_log10pval'},[],...
    'univar_pheno_outlist',{'rMZ','rDZ'},[],...
    'bivar_mx_outlist',{'A1','A2','E','F'},[],...
    'univar_mx_outlist',{'A','C','E','F'},[],...
  ...
    'Ns',[42 162 642 2562 10242 40962 163842],[],...
    'ico',7,[1:7],...
  });

  % check mxbin
  if mmil_isrelative(parms.mxbin) && ~isempty(parms.mxbindir)
    parms.mxbin = [parms.mxbindir '/' parms.mxbin];
  end;
  [s,r] = unix(sprintf('which %s',parms.mxbin));
  if ~isempty(regexp(r,'Command not found'))
    error('MX bin file not found:\n%s',r);
  end;

  % check mxscript
  if mmil_isrelative(parms.mxscript)
    if isempty(parms.mxscriptdir)
      parms.mxscriptdir = [getenv('MMPS_DIR') '/mx'];
    end;
    parms.mxscript = [parms.mxscriptdir '/' parms.mxscript];
  end;
  if ~exist(parms.mxscript,'file')
    error('MX script file %s not found',parms.mxscript);
  end;

  % check input data
  if isempty(parms.fstem_data)
    error('fstem_data is empty');
  elseif ~isstr(parms.fstem_data)
    parms.data = parms.fstem_data'; % input must be [nvox,nsubj]
    parms.fstem_data = [];
    parms.nframes = size(parms.data,1); % nsubj
    parms.nverts = [];
    parms.v_labels = [];
    parms.nvals = size(parms.data,2); % nvox
    parms.ind_included = [1:parms.nframes];
    parms.write_mgh_flag = 0;
  else
    % check input data files
    if mmil_isrelative(parms.fstem_data)
      parms.fstem_data = sprintf('%s/%s',parms.indir,parms.fstem_data);
    end;
    parms.nhemi = length(parms.hemilist);
    parms.fnames_data = [];
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      parms.fnames_data{h} = sprintf('%s-%s%s',...
        parms.fstem_data,hemi,parms.inext);
      if ~exist(parms.fnames_data{h},'file')
        error('file %s not found',parms.fnames_data{h});
      end;
    end;
  end;

  % check input info file
  if mmil_isrelative(parms.fname_info)
    parms.fname_info = sprintf('%s/%s',parms.indir,parms.fname_info);
  end;
  if ~exist(parms.fname_info,'file')
    error('file %s not found',parms.fname_info);
  end;

  % check input reg file
  if ~isempty(parms.fname_reg)
    if mmil_isrelative(parms.fname_reg)
      parms.fname_reg = sprintf('%s/%s',parms.indir,parms.fname_reg);
    end;
    if ~exist(parms.fname_reg,'file')
      error('file %s not found',parms.fname_reg);
    end;
  end;

  % check input label files (e.g. cortex)
  if isempty(parms.fnames_label) && ~isempty(parms.fstem_data)
    if ~strcmpi(parms.label_name,'none')
      if isempty(parms.subjdir)
        fshomedir = getenv('FREESURFER_HOME');
        if isempty(fshomedir)
          error('FREESURFER_HOME environment variable undefined');
        end;
        parms.subjdir = [fshomedir '/subjects'];
      end;
      if ~exist(parms.subjdir,'dir')
        error('FreeSurfer subject dir %s not found',parms.subjdir);
      end;
      fspath = [parms.subjdir '/' parms.subjname];
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        parms.fnames_label{h} = sprintf('%s/label/%s.%s.label',...
          fspath,hemi,parms.label_name);
      end;
    else
      parms.fnames_label = [];
    end;
  end;
  if ~isempty(parms.fnames_label)
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      if ~exist(parms.fnames_label{h},'file')
        error('file %s not found',parms.fnames_label{h});
      end;
    end;
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

  % set output data file name
  parms.outstem_data = parms.outstem;
  parms.fname_data = sprintf('%s_data.mat',parms.outstem_data);

  % set output results file name
  if parms.seed_flag
    parms.outstem = sprintf('%s_sv_%s%d',...
      parms.outstem,parms.seed_hemi,parms.seed_vert);
    if ~isempty(parms.bivar_label)
      fprintf('%s: WARNING: seed_flag = %d, so ignoring bivar_label %s\n',...
        mfilename,parms.seed_flag,parms.bivar_label);
      parms.bivar_label = [];
    end;
  end;
  parms.fname_results = sprintf('%s_results.mat',parms.outstem);

  % set names of output types depending on bivariate vs univariate analysis
  if parms.seed_flag || ~isempty(parms.bivar_label)
    parms.mx_outlist = parms.bivar_mx_outlist;
    parms.pheno_outlist = parms.bivar_pheno_outlist;
  else
    parms.mx_outlist = parms.univar_mx_outlist;
    parms.pheno_outlist = parms.univar_pheno_outlist;
  end;
  if parms.pheno_flag && parms.mx_flag
    parms.outlist = cat(2,parms.pheno_outlist,parms.mx_outlist);
  elseif parms.pheno_flag
    parms.outlist = parms.pheno_outlist;
  elseif parms.mx_flag
    parms.outlist = parms.mx_outlist;
  else
    parms.outlist = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_info(parms)
  if parms.verbose
    fprintf('%s: loading info csv file...\n',mfilename);
  end;

  % merge fname_info and fname_reg
  if ~isempty(parms.fname_reg)
    fname_info = sprintf('%s_info.csv',parms.outstem_data);
    mmil_merge_csv(parms.fname_info,parms.fname_reg,fname_info,...
      parms.merge_field,1,parms.forceflag);
    parms.fname_info = fname_info;
  end;

  % initialize variables
  parms.ind_num_reg = [];
  parms.ind_cat_reg = [];
  parms.ind_quad_reg = [];
  parms.ind_cubic_reg = [];
  parms.ind_outlier_reg = [];
  parms.ind_bivar = [];

  % read csv file(s)
  parms = read_csv(parms);

  % classify column labels and compare against input labels
  if ~isempty(parms.col_labels)
    % determine which columns are numeric and which are categorical
    parms = classify_columns(parms);

    % get column indices for twin pair info
    parms = check_twin_labels(parms);

    % get column index for bivariate regressor
    parms = check_bivar_label(parms);

    % get column indices for redisualization regressors
    parms = check_resid_labels(parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_arrays(parms)
  % check that input cell arrays are cell arrays
  if ~isempty(parms.resid_labels) && ~iscell(parms.resid_labels)
    parms.resid_labels = {parms.resid_labels};
  end;
  if ~isempty(parms.outlier_labels) && ~iscell(parms.outlier_labels)
    parms.outlier_labels = {parms.outlier_labels};
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = read_csv(parms)
  parms.all_vals = [];
  % read fname_info
  if ~isempty(parms.fname_info)
    parms.all_vals = mmil_readtext(parms.fname_info);
    parms.col_labels = parms.all_vals(1,:);
    parms.all_vals = parms.all_vals(2:end,:);
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

function parms = check_twin_labels(parms)
  % get index for pairID column
  parms.ind_pair = find(strcmp(parms.pair_label,parms.col_labels));
  if isempty(parms.ind_pair)
    error('missing column with label "%s"\n',parms.pair_label);
  end;
  % get index for twinID column
  if ~isempty(parms.twin_label)
    parms.ind_twin = find(strcmp(parms.twin_label,parms.col_labels));
    if isempty(parms.ind_twin)
      error('missing column with label "%s"\n',parms.twin_label);
    end;
    % check that there are only two values
    twinIDs = parms.all_vals(:,parms.ind_twin);
    if ismember(parms.ind_twin,parms.ind_cat)
      uniq_vals = unique(twinIDs(~strcmpi(twinIDs,'nan')));
    else
      twinIDs = cell2mat(twinIDs);
      uniq_vals = unique(twinIDs(~isnan(twinIDs)));
    end;
    if length(uniq_vals)~=2
      error('number of unique twinIDs is %d (must be 2)',...
        length(uniq_vals));
    end;
  else
    parms.ind_twin = [];
  end;
  % get index for zygosity column
  parms.ind_zyg = find(strcmp(parms.zyg_label,parms.col_labels));
  if isempty(parms.ind_zyg)
    error('missing column with label "%s"\n',parms.zyg_label);
  end;
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_bivar_label(parms)
  if ~isempty(parms.bivar_label)
    parms.ind_bivar = find(strcmp(parms.bivar_label,parms.col_labels));
    if isempty(parms.ind_bivar) || ~ismember(parms.ind_bivar,parms.ind_num)
      error('bivar_label %s is not a valid numeric column label',...
        parms.bivar_label);
    end;
    % check for overlap between bivar_label and resid_labels
    if ismember(parms.bivar_label,parms.resid_labels)
      error('bivar_label %s also included in resid_labels',...
        parms.bivar_label);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_resid_labels(parms)
  % set indices for different regressors
  if ~isempty(parms.resid_labels)
    % get indices for all regressors excluding mixed columns
    [tmp,parms.ind_reg] = intersect(parms.col_labels,parms.resid_labels);
    parms.ind_reg = mmil_rowvec(parms.ind_reg);
    % check for invalid (missing) reg column labels
    if length(parms.ind_reg) < length(parms.resid_labels)
      inval_labels = setdiff(parms.resid_labels,parms.col_labels);
      if ~isempty(inval_labels)
        fprintf('%s: WARNING: ignoring invalid resid_labels: %s\n',...
          mfilename,sprintf('"%s" ',inval_labels{:}));    
      end;
      parms.resid_labels = parms.col_labels(parms.ind_reg);
    end;
    % get indices for numeric regressors
    parms.ind_num_reg = intersect(parms.ind_reg,parms.ind_num);
    parms.num_resid_labels = parms.col_labels(parms.ind_num_reg);
    % get indices for categorical regressors
    parms.ind_cat_reg = intersect(parms.ind_reg,parms.ind_cat);
    parms.cat_resid_labels = parms.col_labels(parms.ind_cat_reg);
    % demean all num_resid_labels and bivar_label
    parms.demean_labels = union(parms.num_resid_labels,parms.bivar_label);
    % remove any members that are not in num_reg_labels
    parms = exclude_inval_labels(parms,'demean',1);
    % check labels for regressors to be fitted with quadratic function
    if ~isempty(parms.quad_labels)
      % remove any members that are also be in cubic_labels
      if ~isempty(parms.cubic_labels)
        parms.quad_labels = setdiff(parms.quad_labels,parms.cubic_labels);
      end;
      % remove any members that are not in num_resid_labels
      parms = exclude_inval_labels(parms,'quad');
    end;
    % check labels for regressors to be fitted with cubic function
    if ~isempty(parms.cubic_labels)
      % remove any members that are not in num_resid_labels
      parms = exclude_inval_labels(parms,'cubic');
    end;    
    % check labels for regressors to be filtered for outliers
    if ~isempty(parms.outlier_labels)
      % remove any members that are not in num_resid_labels
      parms = exclude_inval_labels(parms,'outlier');
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
    valid_labels = parms.num_resid_labels;
  end;
  % remove any members that are not in valid_labels
  inval_labels = setdiff(parms.(tag),valid_labels);
  if ~isempty(inval_labels)
    if parms.verbose
      fprintf('%s: WARNING: ignoring invalid %s: %s\n',...
        mfilename,tag,sprintf('"%s" ',inval_labels{:}));
    end;
    parms.(tag) = ...
      intersect(parms.(tag),valid_labels);        
  end;
  tag_ind = ['ind_' label_type '_reg'];
  [tmp,parms.(tag_ind)] = ...
    intersect(parms.col_labels,parms.(tag));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_data(parms)
  data = []; nframes = []; nverts = []; v_labels = []; nvals = [];
  % load data from mgz files or existing mat file
  if ~exist(parms.fname_data,'file') || parms.forceflag
    for h=1:parms.nhemi
      % load label file
      if ~isempty(parms.fnames_label)
        if parms.verbose
          fprintf('%s: loading label file %s...\n',...
            mfilename,parms.fnames_label{h});
        end;
        v_labels{h} = sort(fs_read_label(parms.fnames_label{h}));
      end;
      % load data file
      if parms.verbose
        fprintf('%s: loading data from %s...\n',mfilename,parms.fnames_data{h});
      end;
      tmp_data = fs_load_mgh(parms.fnames_data{h});
      volsz = size(tmp_data);
      if length(volsz)<4, error('multiple frames required'); end;
      tmp_data = reshape(tmp_data,[prod(volsz(1:3)),volsz(4)])';
      if isempty(nframes), nframes = volsz(4); end;
      if volsz(4) ~= nframes, error('mismatch in number of frames'); end;
      nverts(h) = size(tmp_data,2);
      % optionally reduce ico level by sub-selection
      if parms.ico ~= 7
        nverts_out = parms.Ns(parms.ico);
        if isempty(v_labels)
          v_labels{h} = [1:nverts_out];
        else
          v_labels{h} = intersect(v_labels{h},[1:nverts_out]);
        end;
      end;
      % select vertices in cortex label
      if ~isempty(v_labels)
        tmp_data = tmp_data(:,v_labels{h});
      end;
      data = cat(2,data,tmp_data);
    end;
    nvals = size(data,2);
    if parms.verbose
      fprintf('%s: saving data to %s...\n',mfilename,parms.fname_data);
    end;
    save(parms.fname_data,...
      'data','nframes','nverts','v_labels','nvals','-v7.3');
  else
    if parms.init_data_flag, return; end;
    if parms.verbose
      fprintf('%s: loading data from %s...\n',mfilename,parms.fname_data);
    end;
    load(parms.fname_data);
  end;
  parms.data = data;
  parms.nframes = nframes;
  parms.nverts = nverts;
  parms.v_labels = v_labels;
  parms.nvals = nvals;
  clear data;
  parms.ind_included = [1:nframes];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_seed(parms)
  if parms.seed_flag
    h = find(strcmp(parms.seed_hemi,parms.hemilist));
    parms.v_seed = find(parms.v_labels{h} == parms.seed_vert);
    if isempty(parms.v_seed)
      error('seed_vert %d not found in %s vertex list',...
        parms.seed_vert,parms.seed_hemi);
    end;
    if h==2, parms.v_seed = length(parms.v_labels{1}) + parms.v_seed; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = exclude_nans(parms)
  if parms.verbose
    fprintf('%s: checking for rows with NaNs...\n',mfilename);
  end;
  ind_required = cat(2,parms.ind_reg,parms.ind_bivar,...
                       parms.ind_pair,parms.ind_twin,parms.ind_zyg);
  nrows = size(parms.all_vals,1);
  valid_flags = ones(nrows,1);
  for i=1:length(ind_required)
    ind_col = ind_required(i);
    if ismember(ind_col,parms.ind_num)
      vals = cell2mat(parms.all_vals(:,ind_col));
      ind_invalid = find(isnan(vals));
      valid_flags(ind_invalid) = 0;
    elseif ismember(ind_col,parms.ind_cat)
      vals = parms.all_vals(:,ind_col);
      ind_invalid = find(strcmpi(vals,'nan'));
      valid_flags(ind_invalid) = 0;
    end;
  end;
  % check data
  ind_invalid = find(isnan(sum(parms.data,2)));  
  valid_flags(ind_invalid) = 0;
  % set valid  
  ind_valid = find(valid_flags);
  % select valid subjects
  parms = select_rows(parms,ind_valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = select_rows(parms,ind_valid)
  nrows = size(parms.all_vals,1);
  nvalid = length(ind_valid);
  if nvalid==0
    error('no valid rows');
  elseif nvalid<nrows
    if parms.verbose
      fprintf('%s: WARNING: excluding %d rows...\n',mfilename,nrows-nvalid);
    end;
    parms.data = parms.data(ind_valid,:);
    parms.all_vals = parms.all_vals(ind_valid,:);
    parms.ind_included = parms.ind_included(ind_valid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = exclude_outliers(parms)
  if isempty(parms.ind_outlier_reg), return; end;
  if parms.verbose
    fprintf('%s: checking for rows with outliers...\n',mfilename);
  end;
  % exclude subjects with outliers
  nrows = size(parms.all_vals,1);
  valid_flags = ones(nrows,1);
  for i=1:length(parms.ind_outlier_reg)
    ind_col = parms.ind_outlier_reg(i);
    vals = cell2mat(parms.all_vals(:,ind_col));
    mean_val = mean(vals);
    std_val = std(vals);
    ind_invalid = find((abs(vals - mean_val)/std_val) > parms.outlier_thresh);
    valid_flags(ind_invalid) = 0;
  end;
  ind_valid = find(valid_flags);
  parms = select_rows(parms,ind_valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = exclude_nontwins(parms)
  if parms.verbose
    fprintf('%s: checking twin information...\n',mfilename);
  end;
  nrows = size(parms.all_vals,1);
  % get pairIDS
  pairIDs = parms.all_vals(:,parms.ind_pair);
  % if pairIDs are numeric, convert to matrix
  if ismember(parms.ind_pair,parms.ind_num)
    pairIDs = cell2mat(pairIDs);
  end;
  % get twinIDs
  if ~isempty(parms.ind_twin)
    twinIDs = parms.all_vals(:,parms.ind_twin);
    % if twinIDs are numeric, convert to matrix
    if ismember(parms.ind_twin,parms.ind_num)
      twinIDs = cell2mat(twinIDs);
    end;
  else
    twinIDs = [];
  end;
  % get zygosity values (note: MZ: zyg = 1; DZ: zyg = 2)
  zygvals = parms.all_vals(:,parms.ind_zyg);
  if ismember(parms.ind_zyg,parms.ind_num)
    zygvals = cell2mat(zygvals);
  end;
  % identify twin pairs
  [uniq_pairIDs,ind_orig,ind_uniq] = unique(pairIDs);
  % exclude subjects without a twin
  ind_twins = [];
  for i=1:length(uniq_pairIDs)
    pairID = uniq_pairIDs(i);
    ind = find(ismember(pairIDs,pairID));
    if length(ind)<2, continue; end;
    if ~isempty(twinIDs) && length(unique(twinIDs(ind)))<2, continue; end;
    if length(unique(zygvals(ind)))==2
      if parms.verbose
        if ~ischar(pairID)
          pairIDstr = num2str(pairID);
        else
          pairIDstr = pairID;
        end;
        fprintf('%s: WARNING: twin pair %s has mismatched zygosity!\n',...
          mfilename,pairIDstr);
      end;
      continue;
    end;
    ind = find(ismember(pairIDs,pairID));
    ind_twins = cat(1,ind_twins,ind);
  end;
  parms = select_rows(parms,ind_twins);
  pairIDs = pairIDs(ind_twins);
  if ~isempty(twinIDs)
    twinIDs = twinIDs(ind_twins);
  end;
  zygvals = zygvals(ind_twins);
  % create vectors of indices for MZa, MZb, DZa, and DZb
  parms.ind_MZa = [];
  parms.ind_MZb = [];
  parms.ind_DZa = [];
  parms.ind_DZb = [];
  uniq_pairIDs = unique(pairIDs);
  for i=1:length(uniq_pairIDs)
    pairID = uniq_pairIDs(i);
    ind = find(ismember(pairIDs,pairID));
    if ~isempty(twinIDs)
      tmp_twinIDs = twinIDs(ind);
    else
      tmp_twinIDs = [1,2];
    end;
    [tmp,ind_sort] = sort(tmp_twinIDs);
    ind = ind(ind_sort);
    tmp_zygval = zygvals(ind(1));
    switch tmp_zygval
      case 1
        parms.ind_MZa = cat(2,parms.ind_MZa,ind(1));
        parms.ind_MZb = cat(2,parms.ind_MZb,ind(2));
      case 2
        parms.ind_DZa = cat(2,parms.ind_DZa,ind(1));
        parms.ind_DZb = cat(2,parms.ind_DZb,ind(2));
      otherwise
        error('invalid zygosity value (%d)',tmp_zygval);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = demean_regs(parms)
  for j=1:length(parms.ind_demean_reg)
    ind = parms.ind_demean_reg(j);
    vals = cell2mat(parms.all_vals(:,ind));
    % subtract mean and normalize by standard deviation
    mean_val = mean(vals);
    std_val = std(vals);
    vals = (vals - mean_val)/std_val;
    parms.all_vals(:,ind) = num2cell(vals);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = scale_data(parms)
  if parms.scalefact==1, return; end;
  if parms.verbose
    fprintf('%s: scaling data...\n',mfilename);
  end;
  parms.data = parms.scalefact * parms.data;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = remove_global_data(parms)
  switch parms.global_flag
    case 1
      % subtract the mean value across vertices
      if parms.verbose
        fprintf('%s: subtracting global mean from data...\n',mfilename);
      end;
      mean_vals = mean(parms.data,2);
      parms.data = bsxfun(@minus,parms.data,mean_vals);
    case 2
      % divide by the sum across vertices
      if parms.verbose
        fprintf('%s: dividing by global sum of data...\n',mfilename);
      end;
      sum_vals = sum(parms.data,2);
      parms.data = bsxfun(@rdivide,parms.data,sum_vals)*mean(sum_vals);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = residualize_data(parms)
  if isempty(parms.ind_num_reg) && isempty(parms.ind_cat_reg)
    return;
  end;
  if parms.verbose
    fprintf('%s: residualizing data...\n',mfilename);
  end;
  % prepare regressor matrix
  parms = prep_reg(parms);
  % construct forward matrix
  nsubs = size(parms.all_vals,1);
  baseflag = isempty(parms.groupregs);
  X = mmil_glm_design(parms.regressors,...
    'nsubs',nsubs,'baseflag',baseflag);
  % calculate inverse matrix
  Xi = pinv(X);
  % calculate beta-weights
  beta = Xi*parms.data;
  % subtract fit from data
  parms.data = parms.data - X*beta;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = zscore_data(parms)
  % z transform data
  parms.data = zscore(parms.data,0,1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_reg(parms)
  parms.regressors = [];
  parms.regnames = {};
  parms.orig_regressors = [];
  parms.orig_regnames = {};
  parms.groupregs = [];
  nvals = size(parms.all_vals,1);
  for j=1:length(parms.ind_num_reg)
    ind = parms.ind_num_reg(j);
    vals = cell2mat(parms.all_vals(:,ind));
    regname = parms.col_labels{ind};
    parms.orig_regressors = cat(2,parms.orig_regressors,vals);
    parms.orig_regnames = cat(2,parms.orig_regnames,regname);
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
      parms.regnames = cat(2,parms.regnames,tmp_regname);
    end;
  end;
  for j=1:length(parms.ind_cat_reg)
    ind = parms.ind_cat_reg(j);
    vals = parms.all_vals(:,ind);
    regname = parms.col_labels{ind};
    % expand categorical regressors
    uniq_vals = unique(vals(~strcmpi(vals,'nan')));
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
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(parms)
  results.nframes = parms.nframes;
  results.nverts = parms.nverts;
  results.v_labels = parms.v_labels;
  results.nvals = parms.nvals;
  results.ind_included = parms.ind_included;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_pheno_corr(parms,results)
  if parms.verbose
    fprintf('%s: calculating MZ and DZ correlations...\n',mfilename);
  end;
  % calculate correlations between MZ twins
  results.rMZ = mmil_corr_mat(parms.data(parms.ind_MZa,:),...
                              parms.data(parms.ind_MZb,:))';
  % calculate correlations between DZ twins
  results.rDZ = mmil_corr_mat(parms.data(parms.ind_DZa,:),...
                              parms.data(parms.ind_DZb,:))';
  % calculate correlation between data and bivar_label values
  if parms.seed_flag || ~isempty(parms.bivar_label)
    if parms.seed_flag
      if parms.verbose
        fprintf('%s: calculating phenotypic correlations with seed %s %d...\n',...
          mfilename,parms.seed_hemi,parms.seed_vert);
      end;
      bivar_vals = parms.data(:,parms.v_seed);
    else
      if parms.verbose
        fprintf('%s: calculating phenotypic correlations with %s...\n',...
          mfilename,parms.bivar_label);
      end;
      bivar_vals = cell2mat(parms.all_vals(:,parms.ind_bivar));
    end;
    [rho,pval] = mmil_corr_mat(bivar_vals,parms.data);
    results.rP = rho';
    results.rP_log10pval = (-sign(rho).*log10(pval))';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_pheno_xcorr(parms,results)
  if parms.verbose
    fprintf('%s: calculating phenotypic cross-correlations...\n',mfilename);
  end;
  results.pheno_xcorr = corr(parms.data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_MX(parms,results)
  % initialize results vectors
  results = init_mx_results(parms,results);
  % create temporary output directory
  mmil_mkdir(parms.mxtmpdir);
  % calculate number of batches
  maxnvals = min(parms.maxnvals,parms.nvals);
  nbatches = ceil(parms.nvals/maxnvals);
  % loop over batches of vertices
  for i=1:nbatches
    if parms.verbose
      fprintf('%s: preparing MX batch %d of %d...\n',mfilename,i,nbatches);
      tic;
    end;
    j = 1 + (i-1)*maxnvals;
    k = min(j + maxnvals - 1,parms.nvals);
    % initialize output files (remove if necessary)  
    [fname_cor,fname_res] = init_mx_output(parms);
    % calculate MZ and DZ covariance for each vertex
    write_mx_input(parms,j,k);
    if parms.dump_mxout_flag
      cmd = sprintf('cd %s; %s < %s > /dev/null',...
        parms.mxtmpdir,parms.mxbin,parms.mxscript);
    else
      cmd = sprintf('cd %s; %s < %s',...
        parms.mxtmpdir,parms.mxbin,parms.mxscript);
    end;
    if parms.verbose
      toc;
    end;
    % run MX calculations
    if parms.verbose
      fprintf('%s: running MX...\n',mfilename);
      tic;
    end;
    [s,r] = unix(cmd);
    if s || ~isempty(strfind(r,'Error'))
      error('MX cmd %s failed:\n%s',cmd,r);
    end;
    if parms.verbose
      toc;
    end;
    if parms.verbose
      fprintf('%s: collecting MX output...\n',mfilename);
    end;
    % read output files
    batch_results = read_mxout(parms,fname_cor,fname_res);
    % repackage results
    results = copy_mx_results(parms,batch_results,results,j,k);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_cor,fname_res] = init_mx_output(parms)
  [tmp,fstem_mxscript] = fileparts(parms.mxscript);
  fname_cor = sprintf('%s/gencor.txt',parms.mxtmpdir);
  fname_res = sprintf('%s/results.txt',parms.mxtmpdir);
  fnames = {fname_cor,fname_res};
  for i=1:length(fnames)
    if exist(fnames{i},'file')
      delete(fnames{i});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_mx_input(parms,j,k)
  % write MZn file
  fname_MZn = sprintf('%s/mzn.txt',parms.mxtmpdir);
  fid = fopen(fname_MZn,'w');
  fprintf(fid,'(*)\n');
  fprintf(fid,'%d\n',length(parms.ind_MZa));
  fclose(fid);
  % write DZn file
  fname_DZn = sprintf('%s/dzn.txt',parms.mxtmpdir);
  fid = fopen(fname_DZn,'w');
  fprintf(fid,'(*)\n');
  fprintf(fid,'%d\n',length(parms.ind_DZa));
  fclose(fid);
  % get covariate values
  if parms.seed_flag
    bivar_vals = parms.data(:,parms.v_seed);
  elseif ~isempty(parms.ind_bivar)
    bivar_vals = cell2mat(parms.all_vals(:,parms.ind_bivar));
  else
    bivar_vals = [];
  end;
  if ~isempty(bivar_vals)
    bivar_MZa = bivar_vals(parms.ind_MZa);
    bivar_MZb = bivar_vals(parms.ind_MZb);
    bivar_DZa = bivar_vals(parms.ind_DZa);
    bivar_DZb = bivar_vals(parms.ind_DZb);
  end;
  % loop over vertices
  kk = j + parms.maxnvals - 1;
  for v=j:kk
    if v<=k
      % get MZ data
      data_MZa = parms.data(parms.ind_MZa,v);
      data_MZb = parms.data(parms.ind_MZb,v);
      % get DZ data
      data_DZa = parms.data(parms.ind_DZa,v);
      data_DZb = parms.data(parms.ind_DZb,v);
      % calculate MZ covariance
      if ~isempty(bivar_vals)
        X = cat(2,bivar_MZa,data_MZa,bivar_MZb,data_MZb);
      else
        X = cat(2,data_MZa,data_MZb);
      end;
      MZcov = nancov(X);
      % calculate DZ covariance
      if ~isempty(bivar_vals)
        X = cat(2,bivar_DZa,data_DZa,bivar_DZb,data_DZb);
      else
        X = cat(2,data_DZa,data_DZb);
      end;
      DZcov = nancov(X);
    else
      MZcov = [1 0; 0 1];
      DZcov = [1 0; 0 1];
    end;
    % write to MZcov file
    fname_MZcov = sprintf('%s/mzcov%d.txt',parms.mxtmpdir,1+v-j);
    fid = fopen(fname_MZcov,'w');
    fprintf(fid,'(*)\n');
    fprintf(fid,'%f\n',MZcov(:));
    fclose(fid);
    % write to DZcov file
    fname_DZcov = sprintf('%s/dzcov%d.txt',parms.mxtmpdir,1+v-j);
    fid = fopen(fname_DZcov,'w');
    fprintf(fid,'(*)\n');
    fprintf(fid,'%f\n',DZcov(:));
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = read_mxout(parms,fname_cor,fname_res)
  results = [];
  % read genetic correlation output file
  if parms.seed_flag || ~isempty(parms.ind_bivar)
    fid = fopen(fname_cor,'r');
    if fid==-1
      error('failed to load file %s',fname_cor);
    end;
    % discard first line
    tmp = fgetl(fid);
    % get values
    mxout = fscanf(fid,'%f %f %f %f',[4 inf]);
    fclose(fid);
    % save results 
    results.A1 = mxout(1,:);
    results.A2 = mxout(2,:);
    % get values
    mxtmp = mmil_readtext(fname_res,'\t');
    % discard first line
    mxtmp = mxtmp(2:end);
    nvals = length(mxtmp);
    % extract values
    mxout = zeros(10,nvals);
    for i=1:nvals
      tmp = regexp(mxtmp{i},'\s+','split');
      tmp = regexprep(tmp,'*+','0');
      tmp = cellfun(@str2num,tmp);
      mxout(:,i) = tmp;
    end;    
    % save results
    Z11 = mxout(7,:);
    Z21 = mxout(8,:);
    Z22 = mxout(9,:);
    % unique environmental
    results.E = (Z11.*Z21)./sqrt((Z11.*Z11).*((Z21.*Z21)+(Z22.*Z22)));
    results.F = mxout(10,:);
  else
    % read results output file
    fid = fopen(fname_res,'r');
    if fid==-1
      error('failed to load file %s',fname_res);
    end;
    % discard first line
    tmp = fgetl(fid);
    mxout = fscanf(fid,'%f %f %f %f',[4 inf]);
    fclose(fid);
    results.A = mxout(1,:);
    results.C = mxout(2,:);
    results.E = mxout(3,:);
    results.F = mxout(4,:);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_mx_results(parms,results)
  for i=1:length(parms.mx_outlist)
    results.(parms.mx_outlist{i}) = nan(parms.nvals,1);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = copy_mx_results(parms,batch_results,results,j,k)
  ind = [j:k];
  nvals = length(ind);
  for i=1:length(parms.mx_outlist)
    results.(parms.mx_outlist{i})(ind) = ...
      batch_results.(parms.mx_outlist{i})(1:nvals);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_mgh_results(parms,results)
  for i=1:length(parms.outlist)
    output = parms.outlist{i};
    fnames = [];
    write_flag = 0;
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      fnames{h} = sprintf('%s_%s-%s%s',parms.outstem,output,hemi,parms.outext);
      if ~exist(fnames{h},'file') || parms.forceflag
        write_flag = 1;
      end;
    end;
    if ~write_flag, continue; end;
    if parms.verbose
      fprintf('%s: writing %s results to mgh...\n',mfilename,parms.outlist{i});
    end;
    vals = results.(parms.outlist{i});
    j = 1;
    for h=1:parms.nhemi
      nverts = results.nverts(h);
      % account for label used to select vertices (e.g. cortex)
      if ~isempty(results.v_labels)
        v_label = results.v_labels{h};
      else
        v_label = [1:nverts];
      end;
      tmp_vals = zeros(nverts,1);
      k = j + length(v_label) - 1;
      % select values for each hemisphere
      tmp_vals(v_label) = vals(j:k);
      % save values to mgh/mgz file
      fs_save_mgh(tmp_vals,fnames{h});
      j = k + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup(parms)
  if ~exist(parms.mxtmpdir,'dir'), return; end;
  cmd = sprintf('rm -r %s',parms.mxtmpdir);
  [s,r] = unix(cmd);
  if s, error('failed to remove mxtmpdir %s:\n%s',parms.mxtmpdir,r); end;
return;

