function errcode = BOLD_MMIL_GLM_ROI_Analysis(ContainerPath,varargin)
%function errcode = BOLD_MMIL_GLM_ROI_Analysis(ContainerPath,[options])
%
% Purpose: Calculate ROI averages from results of GLM analysis on fMRI data
%
% Required Input:
%  ContainerPath: Full path of MRIPROCESSED directory containing event-related fMRI scans
%
% Optional Parameters:
%  'outdir': output directory
%    if empty, will place results in GLM analysis subdirectory in ContainerPath
%    {default = []}
%  'outstem': output file stem
%    if outdir is empty, will append to input GLM analysis file stem
%    {default = []}
%  'label_dir': full path of directory containing one or more label files
%    {default = ContainerPath}
%  'label_stem': file stem for batch of labels (e.g. lh.{stem}_{name}.label)
%    {default = []}
%  'snums': vector of scan numbers used to run GLM analysis
%    If empty, will assume all valid 'for' or 'rev' scans
%     (depending on SessInfo.revflag from BOLD_MMIL_Get_ScanInfo)
%    {default = []}
%  'snums_valid': vector of BOLD scan numbers used in processing
%    If empty, will assume all scans
%    {default = []}
%  'concat_flag' - [0|1] analysis was done after concatenating across multiple runs
%    if 0, will extract values for each scan individually
%    {default: 1}
%  'infix': BOLD file name infix (e.g. '', 'corr', 'corr_resBOLD')
%    {default = []}
%  'norm_flag': [0|1] normalize response function to last condition
%    {default = 1}
%  'iresp_flag': load GLM analysis results from iresp files (impulse response)
%    0: use summary statistics (i.e. area under curve)
%    1: use average of selected iresp time points
%    2: use entire iresp in estimation of response function
%    {default = 0}
%  'iresp_t0': first time point of iresp to include in average
%    {default = 1}
%  'iresp_t1': last time point of iresp to include in average
%    {default = Inf}
%  'iresp_baseline_flag': [0|1] subtract first time point from iresp timecourse
%    {default = 0}
%  'TR': repetition time (only used to calculate time vector for iresp)
%    {default = 1}
%  'weights_flag': [0|1] use weights loaded from mgh files as spatial prior
%    if retfit_flag = 1, weights will be calculated from retfit
%    {default = 0}
%  'weights_dir': full path of directory containing one or more weights files
%    must have corresponding weights file for each label in label_dir
%    If empty, will be set to label_dir
%    {default = []}
%  'combine_ROIs_flag': [0|1] combine results for ROIs with same name
%    {default = 0}
%  'ROI_names': cell array of ROI names to be combined
%    {default = {'v1','v2','v3'}}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional Parameters for using ROIs from retfit
%  'retfit_flag': [0|1] create visual area labels from retfit
%    If 1, label_dir and label_stem will be ignored
%    {default = 0}
%  'retfit_dir': full path of directory containing retfit results
%    {default = [ContainerPath '/retfit']}
%  'retfit_stem': file stem for retfit results
%    {default = 'retfit'}
%  'retfit_thresh': threshold applied to weights (relative to max for one area)
%    {default = 0}
%  'FSPath': FreeSurfer Container Path (required if retfit_flag = 1)
%    {default = []}
%  'cond_info' - struct array containing condition information
%    see rc_read_cond_info
%    If cond_info and fname_conds are empty,
%      will create labels for entire visual areas
%    {default: []}
%  'fname_conds' - full path of csv file containing condition information
%    Ignored if cond_info is supplied
%    {default: []}
%  'r_vec' - vector of eccentricities (degrees visual angle)
%    Ignored if cond_info or fname_conds are supplied
%    {default = 7}
%  'th_vec' vector of polar angles (degrees)
%    Ignored if cond_info or fname_conds are supplied
%    {default = [45,135,225,315]}
%  'ecc_width' - eccentricity width of stimuli (deg. vis. ang.)
%    Ignored if cond_info has ecc_width column
%    {default = 10}
%  'theta_width' - polar angle width of stimuli (degrees)
%    Ignored if cond_info has theta_width column
%    will be ignored if fname_conds has theta_width column
%    {default = 90}
%  'r_max': maximum radius (degrees visual angle) used for eccentricity mapping
%    determines phase for a given eccentricity
%    {default = 12.5}
%  'area_names': names of visual areas
%    {default = {'v1','v2','v3'}}
%  'rf_sizes': vector of receptive field sizes (degrees visual angle)
%    -- one for each visual area
%    {default = [0.66,1.03,1.88]}
%  'rf_slopes': vector of slopes of linear trend of receptive field sizes
%    w.r.t. ecc for each visual area
%    Intercept is assumed to be half of r_max
%    {default = [0.06,0.10,0.15]}
%  'surround_flag': [0|1] whether to model inhibitory surround with difference
%     of Gaussians; Requires vf2ctx_flag = 1
%    {default = 0}
%  'surround_rf_fact': size of surround receptive field relative to center
%    {default = 1.5}
%  'surround_amp_fact': amplitude of surround response relative to center
%    {default = 0.6}
%
% Output:
%   errcode: returns 1 if unsuccessful
%   NOTE: this function creates .mat and .csv files containing roi results
%
% Created:  03/22/11 by Don Hagler
% Last Mod: 12/06/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

errcode = 0;
if ~mmil_check_nargs(nargin,2), return; end;
parms = check_input(ContainerPath,varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main

[ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,...
  'snums',parms.snums_valid,'fnamestem',parms.fnamestem);
if errcode || isempty(ScanInfo), return; end;

% construct cell array of sets of scan numbers
if SessInfo.nscans>1 && parms.concat_flag
  snum_sets = {parms.snums};  
else
  nscans = length(parms.snums);
  snum_sets = cell(nscans,1);
  for i=1:nscans
    snum_sets{i} = parms.snums(i);
  end;
end;

% load label and weights files (create if retfit_flag = 1)
if parms.retfit_flag
  % create label and weights files from retfit
  [ROIs,errcode] = create_retfit_ROIs(ContainerPath,parms);
else
  % find existing label and weights files
  [ROIs,errcode] = load_ROIs(parms);
end;
if errcode, return; end;

% load over sets of scan numbers
for s=1:length(snum_sets)
  [GLM_info,errcode] = load_GLM_info(ContainerPath,snum_sets{s},parms);
  if errcode, return; end;
  
  if isempty(parms.outdir)
    outdir = GLM_info.GLM_dir;
    outstem = GLM_info.GLM_stem;
    if ~isempty(parms.outstem)
      outstem = [outstem '_' parms.outstem];
    end;
  else
    outdir = parms.outdir;
    outstem = parms.outstem;
  end;

  if ~isempty(outstem)
    outstem = [outstem '_'];
  end;

  fname_mat = sprintf('%s/%s%sroi_data.mat',outdir,outstem,parms.label_stem);
  fname_csv = sprintf('%s/%s%sroi_data.csv',outdir,outstem,parms.label_stem);
  if ~exist(fname_mat,'file') || ~exist(fname_csv,'file') || parms.forceflag
    % load data
    [data,errcode] = load_data(ROIs,GLM_info,parms);
    if errcode, return; end;

    % get values for each ROI
    [ROI_data,errcode] = extract_roi_data(data,ROIs,parms);
    if errcode, return; end;

    % combine ROIs
    if parms.combine_ROIs_flag
      [ROI_data,errcode] = combine_ROI_data(ROI_data,parms);
      if errcode, return; end;
    end;

    % calculate responses for each condition
    results = calc_response_function(ROI_data,parms);

    % combine all data into results struct
    results.ROI_data = ROI_data;
    results.data = data;
    results.ROIs = ROIs;
    results.GLM_info = GLM_info;
    results.parms = parms;

    % save results to mat file
    save(fname_mat,'results');

    % save results to csv file
    save_results_csv(fname_csv,results);
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ContainerPath,options)
  parms_filter = { ...
    'outdir',[],[],...
    'outstem',[],[],...
    'label_dir',ContainerPath,[],...
    'label_stem',[],[],...
    'snums',[],[],...
    'snums_valid',[],[],...
    'concat_flag',true,[false true],...
    'infix',[],[],...
    'norm_flag',true,[false true],...
    'iresp_flag',0,[0 1 2],...
    'iresp_t0',1,[1,Inf],...
    'iresp_t1',Inf,[1,Inf],...
    'iresp_baseline_flag',false,[false true],...
    'TR',1,[0.1,100],...
    'weights_flag',false,[false true],...
    'weights_dir',[],[],...
    'combine_ROIs_flag',false,[false true],...
    'ROI_names',{'v1','v2','v3'},[],...
    'forceflag',false,[false true],...
  ... % for ROIs from retfit
    'retfit_flag',false,[false true],...
    'retfit_dir',[ContainerPath '/retfit'],[],...
    'retfit_stem','retfit',[],...
    'retfit_thresh',0,[],...
    'FSPath',[],[],...
    'cond_info',[],[],...
    'fname_conds',[],[],...
    'r_vec',7,[],...
    'th_vec',[45,135,225,315],[],...
    'ecc_width',10,[0,100],...
    'theta_width',90,[0,360],...
    'r_max',12.5,[0,Inf],...
    'rf_sizes',[0.66,1.03,1.88],[],...
    'rf_slopes',[0.06,0.10,0.15],[0,10],...
    'vf2ctx_flag',true,[false true],...
    'stimres',100,[50,1000],...
    'w_thresh',0.01,[0,100],...
    'surround_flag',false,[false true],...
    'surround_rf_fact',1.5,[],...
    'surround_amp_fact',0.6,[],...
    'restrict_hemi_flag',true,[false true],...
    'restrict_uplow_flag',true,[false true],...
  ... % undocumented
    'fnamestem','BOLD',[],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
  ...
    'retfit_label_tags',{'retfit_stem' 'outdir' 'outstem' 'weights_flag'...
                         'cond_info' 'fname_conds' 'r_vec' 'th_vec'...
                         'ecc_width' 'theta_width' 'r_max'...
                         'restrict_hemi_flag' 'restrict_uplow_flag'...
                         'area_names' 'rf_sizes' 'rf_slopes'...
                         'vf2ctx_flag' 'stimres' 'w_thresh' 'retfit_thresh'...
                         'surround_flag' 'surround_rf_fact' ...
                         'surround_amp_fact' 'subjdir' 'forceflag'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if parms.retfit_flag
    if ~exist(parms.retfit_dir,'dir')
      fprintf('%s: ERROR: retfit_flag = 1 but retfit_dir %s not found\n',...
        mfilename,parms.retfit_dir);
      errcode = 1;
      return;
    end;
    if isempty(parms.FSPath)
      fprintf('%s: ERROR: retfit_flag = 1 but FSPath is empty\n',mfilename);
      errcode = 1;
      return;
    end;
    if ~exist(parms.FSPath,'dir')
      fprintf('%s: ERROR: retfit_flag = 1 but FSPath %s not found\n',...
        mfilename,parms.FSPath);
      errcode = 1;
      return;
    end;
    parms.label_stem = parms.retfit_stem;
    parms.area_names = parms.ROI_names;
  end;
  
  if ~isempty(parms.label_stem)
    parms.label_stem = [parms.label_stem '_'];
  end;

  tmp_label_stem = regexprep(parms.label_stem,'_','\_');
  parms.label_regpat = ['(?<hemi>[lr]h).' tmp_label_stem  '(?<stem>[\w\+-\.]+).label'];

  if isempty(parms.weights_dir), parms.weights_dir = parms.label_dir; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GLM_info,errcode] = load_GLM_info(ContainerPath,snums,parms)
  GLM_info = [];
  errcode = 0;
  
  % set GLM analysis file stem
  [GLM_info.full_GLM_stem,errcode] = BOLD_MMIL_Set_GLM_Stem(ContainerPath,...
    'snums',snums,'infix',parms.infix);
  if errcode
    fprintf('%s: ERROR: failed to set GLM stem\n',mfilename);
    return;
  end;
  [GLM_info.GLM_dir,GLM_info.GLM_stem] = fileparts(GLM_info.full_GLM_stem);

  % load stats_info mat file
  clear stats_info;
  fname_info = [GLM_info.full_GLM_stem '.mat'];
  if ~exist(fname_info,'file')
    fprintf('%s: ERROR: file %s not found\n',mfilename,fname_info);
    errcode = 1;
    return;
  end;
  load(fname_info);
  if ~exist('stats_info','var')
    fprintf('%s: ERROR: file %s did not contain ''stats_info'' as expected\n',...
      mfilename,fname_info);
    errcode = 1;
    return;
  end;
  GLM_info.stats_info = stats_info;

  % find coefficient indices and number of conditions
  GLM_info.ind_coef = find(strcmp('coef',{GLM_info.stats_info.type}));
  GLM_info.ind_stat = find(strcmp('fstat',{stats_info.type}) &...
    ~strcmp('Full',{stats_info.name}));
  GLM_info.nconds = length(GLM_info.ind_coef);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ROIs,errcode] = create_retfit_ROIs(ContainerPath,parms)
  ROIs = [];
  errcode = 0;

  % set parameters
  [parms.subjdir,subj,text] = fileparts(parms.FSPath);
  subj = [subj text];
  if isempty(parms.outdir)
    parms.outdir = ContainerPath;
  end;
  parms.outdir = [parms.outdir '/retfit_ROIs'];
  parms.outstem = parms.retfit_stem;
  if parms.surround_flag
    parms.outstem = sprintf('%s_srf%0.1f_saf%0.1f',...
      parms.outstem,parms.surround_rf_fact,parms.surround_amp_fact);
  end;
  if parms.retfit_thresh>0
    parms.outstem = sprintf('%s_thresh%0.2f',parms.outstem,parms.retfit_thresh);
  end;
  parms.retfit_stem = [parms.retfit_dir '/matfiles/' parms.retfit_stem];
  args = mmil_parms2args(parms,parms.retfit_label_tags);

  try
    [label_files,weights_files] = rc_create_retfit_labels(subj,args{:});
  catch
    fprintf('%s: ERROR: failed to create retfit labels:\n%s\n',...
      mfilename,lasterr);
    errcode = 1;
    return;
  end;
  
  % load labels
  for r=1:length(label_files)
    [tmp,fstem,fext] = fileparts(label_files{r});
    fname = [fstem fext];
    [ROIs(r).name,ROIs(r).hemi] = ...
      get_ROI_name(fname,parms.label_regpat);
    ROIs(r).fname_label = label_files{r};
    ROIs(r).verts = fs_read_label(ROIs(r).fname_label);
  end;

  % load weights
  if parms.weights_flag
    for r=1:length(ROIs)
      ROIs(r).fname_weights = weights_files{r};
      ROIs(r).weights = fs_load_mgh(ROIs(r).fname_weights);
      ROIs(r).weights = ROIs(r).weights(ROIs(r).verts);
    end;
  end;

  if isempty(ROIs)
    fprintf('%s: ERROR: ROIs structure is empty\n',mfilename);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ROIs,errcode] = load_ROIs(parms)
  ROIs = [];
  errcode = 0;

  r = 1;
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};

    % find all label files with label_stem
    flist = dir(sprintf('%s/%s.%s*.label',...
      parms.label_dir,hemi,parms.label_stem));
    if isempty(flist)
      fprintf('%s: ERROR: %s.%s label files not found in %s',...
        hemi,parms.label_stem,parms.label_dir);
      errcode = 1;
      return;
    end;
    fnames = {flist.name};

    % exclude labels that have been resampled to ico
    n = regexp(fnames,'ico\d\.label$');
    if ~isempty(n)
      ind_keep = find(cellfun(@isempty,n));
      fnames = fnames(ind_keep);
    end;
    if isempty(fnames)
      error('no non-ico labels found in %s',parms.label_dir);
    end;

    % load labels
    for f=1:length(fnames)
      [ROIs(r).name,ROIs(r).hemi] = ...
        get_ROI_name(fnames{f},parms.label_regpat);
      ROIs(r).fname_label = [parms.label_dir '/' fnames{f}];
      ROIs(r).verts = fs_read_label(ROIs(r).fname_label);
      r = r + 1;
    end;
  end;
  
  % load weights
  if parms.weights_flag
    for r=1:length(ROIs)
      % find corresponding weights file for each label
      ROIs(r).fname_weights = sprintf('%s/%s-%s.mgh',...
        parms.weights_dir,ROIs(r).name,ROIs(r).hemi);
      if ~exist(ROIs(r).fname_weights,'file')
        fprintf('%s: ERROR: weights file %s not found in %s\n',...
          mfilename,ROIs(r).fname_weights,parms.weights_dir);
        errcode = 1;
        return;
      end;
      ROIs(r).weights = fs_load_mgh(ROIs(r).fname_weights);
      ROIs(r).weights = ROIs(r).weights(ROIs(r).verts);
    end;
  end;

  if isempty(ROIs)
    fprintf('%s: ERROR: ROIs structure is empty\n',mfilename);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [name,hemi] = get_ROI_name(fname,regpat)
  name = []; hemi = [];
  n = regexp(fname,regpat,'names');
  if isempty(n)
    fprintf('%s: WARNING: ROI file %s does not match expected pattern (%s)',...
      mfilename,fname,regpat);
    name = fname;
  else
    name = n.stem;
    hemi = n.hemi;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,errcode] = load_data(ROIs,GLM_info,parms)
  data = [];
  errcode = 0;
  
  data.nconds = GLM_info.nconds;
  data.hemilist = parms.hemilist;
  
  % loop over hemispheres
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};

    % load statistics (F-stats and area under curve coefficients)
    fname_data = [GLM_info.full_GLM_stem '-' hemi '.mgh'];
    if ~exist(fname_data,'file')
      fprintf('%s: ERROR: file %s not found',mfilename,fname_data);
      errcode = 1;
      return;
    end;
    data.hemi_data(h).fname = fname_data;
    data.hemi_data(h).hemi = hemi;
    vals = squeeze(fs_load_mgh(fname_data));
    data.hemi_data(h).coefs = vals(:,GLM_info.ind_coef);
    data.hemi_data(h).stats = vals(:,GLM_info.ind_stat);
    data.hemi_data(h).nverts = size(data.hemi_data(h).coefs,1);

    if parms.iresp_flag
      for c=1:GLM_info.nconds
        fname_data = sprintf('%s_iresp_cond%d-%s.mgh',...
          GLM_info.full_GLM_stem,c,hemi);
        if ~exist(fname_data,'file')
          fprintf('%s: ERROR: file %s not found',mfilename,fname_data);
          errcode = 1;
          return;
        end;
        data.hemi_data(h).iresp_data(c).fname = fname_data;
        data.hemi_data(h).iresp_data(c).vals = squeeze(fs_load_mgh(fname_data));
        [data.hemi_data(h).iresp_data(c).nverts,...
         data.hemi_data(h).iresp_data(c).ntpoints] = ...
            size(data.hemi_data(h).iresp_data(c).vals);
      end;
    end;
  end;

  if isempty(data)
    fprintf('%s: ERROR: data structure is empty\n',mfilename);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ROI_data,errcode] = extract_roi_data(data,ROIs,parms)
  ROI_data = ROIs;
  errcode = 0;
  
  % extract values for each ROI
  for r=1:length(ROI_data)
    ROI_data(r).nconds = data.nconds;
    v = ROI_data(r).verts;
    h = find(strcmp(ROI_data(r).hemi,parms.hemilist));
    ROI_data(r).nverts = length(v);
    ROI_data(r).coefs = data.hemi_data(h).coefs(v,:);
    ROI_data(r).stats = data.hemi_data(h).stats(v,:);
    ROI_data(r).mean_coefs = mean(ROI_data(r).coefs,1);
    ROI_data(r).mean_stats = mean(ROI_data(r).stats,1);

    if parms.norm_flag % normalize to last condition
      ROI_data(r).mean_coefs = ...
        ROI_data(r).mean_coefs/ROI_data(r).mean_coefs(data.nconds);
    end;
    
    if parms.iresp_flag
      for c=1:data.nconds
        ntpoints = data.hemi_data(h).iresp_data(c).ntpoints;
        t0 = max(parms.iresp_t0,1);
        t1 = min(parms.iresp_t1,ntpoints);
        roi_vals = ...
          data.hemi_data(h).iresp_data(c).vals(v,:);
        if parms.iresp_baseline_flag
          base_vals = roi_vals(:,1);
          roi_vals = roi_vals - base_vals*ones(1,ntpoints);
        end;
        ROI_data(r).iresp_data(c).vals = roi_vals;
        ROI_data(r).iresp_data(c).ntpoints = ntpoints;
        ROI_data(r).iresp_data(c).time = ...
          parms.TR*[0:data.hemi_data(h).iresp_data(c).ntpoints-1];
        ROI_data(r).iresp_data(c).mean_vals = mean(roi_vals(:,t0:t1),2);
        ROI_data(r).iresp_data(c).mean_timecourse = mean(roi_vals,1);
      end;    
    end;
  end;
  if isempty(ROI_data)
    fprintf('%s: ERROR: ROI_data is empty\n',mfilename);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [new_ROI_data,errcode] = combine_ROI_data(ROI_data,parms)
  % combine across ROIs with same name
  new_ROI_data = [];
  errcode = 0;
  
  ROI_names = {ROI_data.name};
  nconds = ROI_data(1).nconds;
  for a=1:length(parms.ROI_names)
    new_ROI_data(a).name = parms.ROI_names{a};
    n = regexp(ROI_names,new_ROI_data(a).name);
    ind = find(~cellfun(@isempty,n));
    if isempty(ind)
      fprintf('%s: WARNING: no ROIs matching %s\n',...
        mfilename,new_ROI_data.name);
      continue;
    end;
    new_ROI_data(a).nconds = nconds;
    new_ROI_data(a).nverts = 0;
    new_ROI_data(a).coefs = [];
    new_ROI_data(a).stats = [];
    new_ROI_data(a).hemis = [];
    for r=ind
      h = find(strcmp(ROI_data(r).hemi,parms.hemilist));
      hemis = h*ones(ROI_data(r).nverts,1);
      new_ROI_data(a).hemis = cat(1,new_ROI_data(a).hemis,hemis);
      new_ROI_data(a).coefs = cat(1,new_ROI_data(a).coefs,ROI_data(r).coefs);
      new_ROI_data(a).stats = cat(1,new_ROI_data(a).stats,ROI_data(r).stats);
    end;
    new_ROI_data(a).nverts = size(new_ROI_data(a).coefs,1);
    new_ROI_data(a).mean_coefs = mean(new_ROI_data(a).coefs,1);
    new_ROI_data(a).mean_stats = mean(new_ROI_data(a).stats,1);

    if parms.norm_flag % normalize to last condition
      new_ROI_data(a).mean_coefs = ...
        new_ROI_data(a).mean_coefs/new_ROI_data(a).mean_coefs(nconds);
    end;
    
    if parms.iresp_flag
      for c=1:nconds
        ntpoints = ROI_data(r).iresp_data(c).ntpoints;
        t0 = max(parms.iresp_t0,1);
        t1 = min(parms.iresp_t1,ntpoints);
        time = ROI_data(r).iresp_data(c).time;
        roi_vals = [];
        for r=ind
          roi_vals = ...
            cat(1,roi_vals,ROI_data(r).iresp_data(c).vals);
        end;
        new_ROI_data(a).iresp_data(c).vals = roi_vals;
        new_ROI_data(a).iresp_data(c).ntpoints = ntpoints;
        new_ROI_data(a).iresp_data(c).time = time;
        new_ROI_data(a).iresp_data(c).mean_vals = mean(roi_vals(:,t0:t1),2);
        new_ROI_data(a).iresp_data(c).mean_timecourse = mean(roi_vals,1);
      end;
    end;

    if parms.weights_flag
      roi_weights = [];
      for r=ind
        roi_weights = ...
          cat(1,roi_weights,ROI_data(r).weights);
      end;
      new_ROI_data(a).weights = roi_weights;
    end;
  end;
  if isempty(new_ROI_data)
    fprintf('%s: ERROR: ROI_data is empty\n',mfilename);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate response function
function results = calc_response_function(ROI_data,parms)
  results = [];
  nroi = length(ROI_data);
  nconds = ROI_data(1).nconds;

  % create model of hemodynamic response function from average across all areas
  if parms.iresp_flag==2
    ntpoints = ROI_data(1).iresp_data(1).ntpoints;
    mean_hrf = zeros(1,ntpoints);
    j = 1;
    for r=1:nroi
      for c=1:nconds
        tmp_hrf = ROI_data(r).iresp_data(c).mean_timecourse;
        mean_hrf = mean_hrf + tmp_hrf;
        j = j + 1;
      end;
    end;
    mean_hrf = mean_hrf / j;  
    mean_hrf = (mean_hrf - min(mean_hrf)) / (max(mean_hrf) - min(mean_hrf));
    mean_hrf(end) = 0;
  else
    ntpoints = 1;
    mean_hrf = [];
  end;

  resp = zeros(nconds,nroi);
  for r=1:nroi
    nverts = ROI_data(r).nverts;
    if parms.iresp_flag==2
      % construct forward matrix with hemodynamic response function
      nmeas = nverts*nconds*ntpoints;
      X = zeros(nmeas,nconds);
      j = 1;
      for c=1:nconds
        for t=1:ntpoints
          k = j + nverts - 1;
          if parms.weights_flag
            % temporal prior from hemodynamic response
            %  and spatial prior from weights
            X(j:k,c) = mean_hrf(t)*ROI_data(r).weights;
          else
            % temporal prior from hemodynamic response
            X(j:k,c) = mean_hrf(t);
          end;
          j = j + nverts;
        end;
      end;
      % combine data from multiple conditions and time points
      Y = zeros(nmeas,1);
      j = 1;
      for c=1:nconds
        for t=1:ntpoints
          k = j + nverts - 1;
          Y(j:k) = ROI_data(r).iresp_data(c).vals(:,t);
          j = j + nverts;
        end;
      end;
    else
      % construct forward matrix for single time point
      nmeas = nverts*nconds;
      X = zeros(nmeas,nconds);
      for c=1:nconds
        j = 1 + (c-1)*nverts;
        k = j + nverts - 1;
        if parms.weights_flag
          % spatial prior from weights
          X(j:k,c) = ROI_data(r).weights;
        else
          % simple average
          X(j:k,c) = 1;
        end;
      end;  
      Y = zeros(nmeas,1);
      for c=1:nconds
        j = 1 + (c-1)*nverts;
        k = j + nverts - 1;
        if parms.iresp_flag
          Y(j:k) = ROI_data(r).iresp_data(c).mean_vals;
        else
          Y(j:k) = ROI_data(r).coefs(:,c);
        end;
      end;  
    end;
    beta = pinv(X)*Y;
    if parms.norm_flag
      resp(:,r) = beta/beta(nconds); % normalize to last condition
    else
      resp(:,r) = beta;
    end;
  end;
  
  results.nroi = nroi;
  results.nconds = nconds;
  results.ntpoints = ntpoints;
  results.mean_hrf = mean_hrf;
  results.resp = resp;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_results_csv(fname_csv,results)
  fid = fopen(fname_csv,'wt');
  if fid<0, error('failed to open %s for writing',fname_csv); end;
  fprintf(fid,'"ROI name"');
  nstats = length(results.GLM_info.stats_info);
  for i=results.GLM_info.ind_coef
    fprintf(fid,',"%s-resp"',...
      results.GLM_info.stats_info(i).name);
  end;
  for i=results.GLM_info.ind_coef
    fprintf(fid,',"%s-%s"',...
      results.GLM_info.stats_info(i).name,results.GLM_info.stats_info(i).type);
  end;
  for i=results.GLM_info.ind_stat
    fprintf(fid,',"%s-%s"',...
      results.GLM_info.stats_info(i).name,results.GLM_info.stats_info(i).type);
  end;
  fprintf(fid,',"nverts"\n');
  for r=1:results.nroi
    fprintf(fid,'"%s"',results.ROI_data(r).name);
    for c=1:results.nconds
      fprintf(fid,',%0.4f',results.resp(c,r));
    end;
    for c=1:results.nconds
      fprintf(fid,',%0.4f',results.ROI_data(r).mean_coefs(c));
    end;
    for c=1:results.nconds
      fprintf(fid,',%0.4f',results.ROI_data(r).mean_stats(c));
    end;
    fprintf(fid,',%d\n',results.ROI_data(r).nverts);
  end;
  fclose(fid);
return;


