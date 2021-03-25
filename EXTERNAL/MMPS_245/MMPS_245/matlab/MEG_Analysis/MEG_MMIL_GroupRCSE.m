function MEG_MMIL_GroupRCSE(ProjID,varargin)
%function MEG_MMIL_GroupRCSE(ProjID,varargin)
%
% Purpose: Retinotopy Constrained Source Estimation for multiple subjects
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%    Must have 'MEG_VisitID' to indicate name of the orig_meg data directory
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include this field: proc_meg
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'hiC_flag': [0|1] use source results for highest condition only
%     if 0, all subjects must have same number of condition levels
%    {default = 0}
%  'reweight_flag': [0|1] whether to use iteratively re-weighted least squares
%    {default = 0}
%  'reweight_init_flag': [0|1] whether to use single subject IRLS weights
%     for initial RCSE calculations
%    {default = 1}
%  'reweight_factor': tunable reweighting factor or IRLS
%    {default = 2}
%  'inverse_type': [0|1|2] type of inverse calculations
%    0: unregularized pseudo-inverse (fast, very little memory)
%    1: regularized psuedo-inverse with identity matrix for noise covariance
%    2: regularized psuedo-inverse with noise covariance matrix concatenated
%       across subjects (huge amount of memory (i.e. > 48 GB))
%    {default = 1}
%  'SNR': estimated signal to noise ratio (for regularization)
%    {default = 10}
%  'resamp_flag': [0|1] whether to use bootstrap resampling
%      to calculate standard error and 95% confidence intervals
%    {default = 0}
%  'resamp_niters': number of bootstrap resampling iterations
%    {default = 100}
%  'randseed_flag': [0|1] for resampling, generate different result every time
%    {default = 0}
%  'fill_err_flag': [0|1] when plotting group waveforms, fill confidence
%     interval as shaded region; otherwise use error bars
%    {default = 0}
%  'outdir': output directory
%    full path or relative to /home/{user}/MetaData/{ProjID}
%    if ProjID is empty, will be relative to current working directory
%    {default = 'GroupRCSE'}
%  'outstem': output file stem (relative to outdir)
%    {default = 'RCSE'}
%  'qcflag': [0|1] exclude subjects with StudyInfo.QC = 0
%    {default = 1}
%  'plotflag': [0|1] plot results
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional Input from ProjInfo, command line, or StudyInfo
%  'RCSE_prefix': RCSE prefix
%    {default = 'RCSE'}
%  'RCSE_infix': extra string attached to RCSE output prefix
%    {default = []}
%  'RCSE_infix_flag': whether infix was automatically generated
%    ignored if RCSE_infix is not empty
%    {default = 0}
%  'RCSE_fstem_conds': stem of csv file containing condition information
%    {default = 'cond_info'}
%
% Created:  03/07/11 by Don Hagler
% Last Mod: 03/30/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
...
  'RCSE_prefix','RCSE',[],...
  'RCSE_infix',[],[],...
  'RCSE_infix_flag',false,[false true],...
  'RCSE_fstem_conds','cond_info',[],...
...
  'hiC_flag',false,[false true],...
  'reweight_flag',false,[false true],...
  'reweight_init_flag',true,[false true],...
  'reweight_factor',2,[],...
  'inverse_type',1,[0,1,2],...
  'SNR',10,[],...
  'resamp_flag',false,[false true],...
  'resamp_niters',100,[1,100000],...
  'randseed_flag',false,[false true],...
  'fill_err_flag',false,[false true],...
  'ci_tail',2.5,[0.1,50],...
  'outdir','GroupRCSE',[],...
  'outstem','RCSE',[],...
  'qcflag',true,[false true],...
  'plotflag',true,[false true],...
  'forceflag',false,[false true],...  
... for waveform analysis
  'smooth_sigma',0,[0,1000],...
  'auc_range',[0,350],[],...
  'auc_nbins',1,[1,100],...
  'auc_baseline_flag',true,[false true],...
  'auc_baseline_range',[-100,0],[],...
  'peak_range',[50,150],[],...
  'peak_pol',-1,[-1,1,1],...
  'peak_mindiff',0.5,[0,Inf],...
  'normflag',1,[0 1 2],...
  'powerflag',2,[0 1 2],...
  'firstflag',false,[false true],...
  'sfreq',1000,[],...
  't0',-100,[],...
  't1',350,[],...
  'time',[],[],...
... % for plotting analysis results
  'condition_values',[],[],...
  'condition_label','Stimulus Condition',[],...
  'eps_flag',false,[false true],...
  'visible_flag',false,[false true],...
  'linewidth',1,[],...
  'min_linewidth',1,[],...
  'max_linewidth',2.5,[],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
  'ylim_auc',[0,1.3],[],...
  'ylim_peak',[0,1.3],[],...
  'ylim_latency',[50,150],[],...
  'ylim_wform',[-25,10],[],...
  'xlim',[0,1.05],[],...
  'units_wform','nA M',[],...
  'logx_flag',false,[false true],...
  'stderrflag',true,[false true],...
  'area_names',{'V1','V2','V3'},[],...
  'area_colors',{'b','g','r'},[],...
  'legend_loc','SouthEast',[],...
  'label_flag',true,[false true],...
...
  'reweight_maxiter',100,[],...
  'reweight_tol',1e-7,[],...
  'reweight_leverage_flag',false,[false true],... % uses a lot of memory
  'reweight_leverage_max_flag',true,[false true],...
...
  'plot_text_dx',0.01,[],...
  'plot_text_dy',0,[],...
...
  'required_containers',{'proc_meg'},[],...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[RootDirs,StudyInfo,parms] = check_input(ProjID,parms_filter,varargin);

fname_results = [parms.outstem '_results.mat'];
results = [];
if ~exist(fname_results,'file') || parms.forceflag
  data = [];
  data.time = parms.time;
  data.ntpoints = parms.ntpoints;

  % load parms, data, forward for all subjects
  fname = [parms.outstem '_data.mat'];
  if ~exist(fname,'file') || parms.forceflag
    % check RCSE parameters for each subject
    [subj_parms,StudyInfo] = check_subj_RCSE(RootDirs,StudyInfo,parms);
    if isempty(subj_parms), return; end;

    parms.nareas = subj_parms(1).nareas;
    data.nareas = parms.nareas;
    data.areas = subj_parms(1).use_areas;
    parms.ncontrasts_orig = subj_parms(1).ncontrasts;
    data.contrasts = subj_parms(1).unique_contrasts;
    if parms.hiC_flag
      parms.ncontrasts = 1;
      data.contrasts = data.contrasts(end);
    else
      parms.ncontrasts = parms.ncontrasts_orig;
    end;
    data.ncontrasts_orig = parms.ncontrasts_orig;
    data.ncontrasts = parms.ncontrasts;
    data.nwforms = data.nareas * data.ncontrasts;
    data.nsubjects = length(subj_parms);
    data.MEG_VisitIDs = {StudyInfo.MEG_VisitID};

    % load and concatenate data for all subjects
    [data.Y,data.condvec,data.subjvec] = ...
      concat_data(RootDirs,StudyInfo,parms,subj_parms);

    % load and concatenate forward matrices for all subjects
    data.F = concat_forward(RootDirs,StudyInfo,parms);

    save(fname,'data','subj_parms','StudyInfo','-v7.3');
  else
    load(fname);
    parms.nareas = data.nareas;
    parms.ncontrasts = data.ncontrasts;
  end;

  fprintf('%s: calculating source estimates...\n',mfilename);
  if parms.reweight_flag
    results = calc_waveforms_IRLS(data,parms,subj_parms);
    % plot IRLS error
    if parms.plotflag
      plot_IRLS_weights(results,parms);
    end;
  else
    results = calc_waveforms(data,parms,subj_parms);
  end;

  if parms.resamp_flag
    S = zeros(results.ntpoints,results.nareas,...
      results.ncontrasts,parms.resamp_niters);
    for i=1:parms.resamp_niters
      subj_inds = randi(parms.nsubs,parms.nsubs,1);
      [tmp_data,tmp_subj_parms] = bootstr_resamp_data(subj_inds,data,subj_parms);
      fprintf('%s: calculating source estimates iteration %d...\n',mfilename,i);
      if parms.reweight_flag
        tmp_results = calc_waveforms_IRLS(tmp_data,parms,tmp_subj_parms);
      else
        tmp_results = calc_waveforms(tmp_data,parms,tmp_subj_parms);
      end;
      % compile results
      S(:,:,:,i) = tmp_results.S;
    end;
    % calculate confidence intervals
    results.S_resamp = S;
    results.S_sem = std(results.S_resamp,1,4);
    results.S_lo = prctile(results.S_resamp,parms.ci_tail,4);
    results.S_hi = prctile(results.S_resamp,100-parms.ci_tail,4);
    results.S_lo_sem = prctile(results.S_resamp,100*normcdf(-1),4);
    results.S_hi_sem = prctile(results.S_resamp,100*normcdf(1),4);
  end;

  % plot waveforms
  fprintf('%s: plotting waveforms...\n',mfilename);
  plot_waveforms(results,parms);

  % analyze waveforms
  fprintf('%s: analyzing waveforms...\n',mfilename);
  results = analyze_waveforms(results,parms);

  fprintf('%s: saving results to %s...\n',mfilename,fname_results);
  save(fname_results,'results');
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RootDirs,StudyInfo,parms] = check_input(ProjID,parms_filter,options)
  parms = mmil_args2parms(options,parms_filter);
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if ~isempty(ProjInfo)
    % For arg names present in both varargin and ProjInfo
    % the varargin values will appear in merged_args
    ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
    merged_args = mmil_merge_args(options,ProjInfo_args);
    % check that parameters fit allowed range, use defaults if not supplied
    parms = mmil_args2parms(merged_args,parms_filter);
  end;

  parms.nsubs  = length(StudyInfo);
  if ~parms.nsubs
    fprintf('%s: ERROR: no valid subjects in StudyInfo\n');
    return;
  end;

  if isfield(StudyInfo,'GroupAvg')
    for i=1:parms.nsubs
      if isempty(StudyInfo(i).GroupAvg)
        StudyInfo(i).GroupAvg = 0;
      end;
    end;
    ind_GroupAvg = find([StudyInfo.GroupAvg]);
    StudyInfo = StudyInfo(ind_GroupAvg);
    parms.nsubs  = length(StudyInfo);
    if ~parms.nsubs
      fprintf('%s: ERROR: no valid subjects with StudyInfo.GroupAvg=1\n');
      return;
    end;
  end;

  if ~isfield(StudyInfo,'RCSE_fstem_conds')
    for i=1:parms.nsubs
      StudyInfo(i).RCSE_fstem_conds = parms.RCSE_fstem_conds;
    end;
  end;

  if ~isfield(StudyInfo,'RCSE_infix')
    for i=1:parms.nsubs
      StudyInfo(i).RCSE_infix = parms.RCSE_infix;
    end;
  end;

  % set RCSE_prefix in case emtpy (as would be done by rc_RCSE)
  for s=1:parms.nsubs
    infix = StudyInfo(s).RCSE_infix;
    if parms.RCSE_infix_flag && isempty(infix)
      infix = rc_RCSE_set_infix(parms);
    end;
    if ~isfield(StudyInfo,'RCSE_prefix') || isempty(StudyInfo(s).RCSE_prefix)
      StudyInfo(s).RCSE_prefix = parms.RCSE_prefix;
    end;
    if ~isempty(infix)
      StudyInfo(s).RCSE_prefix = [StudyInfo(s).RCSE_prefix '_' infix];
    end;
  end;

  if isempty(parms.time)
    samp_dur = 1000/parms.sfreq;
    parms.time = [parms.t0:samp_dur:parms.t1]/1000;
  end;
  parms.ntpoints = length(parms.time);

  % change outstem to be full path
  if mmil_isrelative(parms.outdir)
    if isempty(ProjID)
      parms.outdir = [pwd '/' parms.outdir];
    else
      parms.outdir = [RootDirs.home '/MetaData/' ProjID '/' parms.outdir];
    end;
  end;
  mmil_mkdir(parms.outdir);
  parms.outstem = [parms.outdir '/' parms.outstem];

  if parms.randseed_flag
    seedval = sum(100*clock);
  else
    seedval = 5489; % default
  end;
  stream = RandStream.create('mt19937ar','Seed',seedval);
  RandStream.setDefaultStream(stream);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [subj_parms,StudyInfo] = check_subj_RCSE(RootDirs,StudyInfo,parms)
  fprintf('%s: checking RCSE parameters and results for each subject...\n',...
    mfilename);
  infix_list = {'parms','avg_data','forward_prep','ret_forward'};
  subj_parms = [];
  keep_flags = ones(1,parms.nsubs);
  for s=1:parms.nsubs
    % check input files exist
    indir = [RootDirs.proc_meg '/' StudyInfo(s).proc_meg];
    for i=1:length(infix_list)
      fname = sprintf('%s/matfiles/%s_%s.mat',...
        indir,StudyInfo(s).RCSE_prefix,infix_list{i});
      if ~exist(fname,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
        keep_flags(s) = 0;
        break;
      end;
    end;
    if ~keep_flags(s), continue; end;
    % load RCSE parameters
    fname = sprintf('%s/matfiles/%s_parms.mat',indir,StudyInfo(s).RCSE_prefix);
    tmp = load(fname);
    tmp_parms = tmp.parms;
    % get scale matrix to be applied to data
    fname = sprintf('%s/matfiles/%s_forward_prep.mat',indir,StudyInfo(s).RCSE_prefix);
    tmp = load(fname);
    tmp_parms.scale_matrix = single(tmp.forward.scale_matrix);
    % get RCSE results for IRLS weights
    if parms.reweight_init_flag && parms.reweight_flag && tmp_parms.reweight_flag
      fname = sprintf('%s/matfiles/%s_results.mat',...
        indir,StudyInfo(s).RCSE_prefix);
      tmp = load(fname);
      tmp_parms.cond_weights = tmp.results.cond_weights;
    end;
    subj_parms = [subj_parms tmp_parms];
  end;
  StudyInfo = StudyInfo(keep_flags>0);

  % check that all subjects have right number of visual areas
  nareas_vec = [subj_parms.nareas];
  if any(nareas_vec~=nareas_vec(1))
    fprintf('%s: ERROR: number of visual areas do not match\n',...
      mfilename);
    return;
    subj_parms = [];
  end;

  % check that all subjects have same number of contrast levels
  ncontrasts_vec = [subj_parms.ncontrasts];
  if min(ncontrasts_vec)~=max(ncontrasts_vec)
    fprintf('%s: ERROR: number of contrast levels varies between subjects (min = %d, max = %d)\n',...
      mfilename,min(ncontrasts_vec),max(ncontrasts_vec));
    return;
    subj_parms = [];
  end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = concat_forward(RootDirs,StudyInfo,parms)
  fprintf('%s: concatenating forward solutions for all subjects...\n',mfilename);
  F = [];
  for s=1:parms.nsubs
    indir = [RootDirs.proc_meg '/' StudyInfo(s).proc_meg];
    fname = sprintf('%s/matfiles/%s_ret_forward.mat',...
      indir,StudyInfo(s).RCSE_prefix);
    if ~exist(fname,'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
      continue;
    end;
    load(fname);
    F = cat(1,F,single(retforward.F));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = calc_inverse(F,parms,subj_parms)
  fprintf('%s: calculating inverse matrix...\n',mfilename);
  switch parms.inverse_type
    case 0
      W = rc_calc_simple_inverse(F);
    case 1
      W = rc_calc_ret_inverse_operator_noprior(F,parms.SNR);
    case 2
      W = [];
      num_measurements = size(F,1);
      num_sources = size(F,2);
      R = speye(num_sources,num_sources);
      C = sparse(num_measurements,num_measurements);
      for s=1:length(subj_parms)
        tmpC = subj_parms(s).noisecovar;
        num_sensors = length(subj_parms(s).goodchans);
        num_conds = length(subj_parms(s).unique_location_conds);
        for i=1:num_conds
          j = 1 + (i-1)*num_sensors;
          k = j + num_sensors - 1;
          C(j:k,j:k) = tmpC;
        end;
      end;
      W = rc_calc_ret_inverse_operator(F,parms.SNR,R,C);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y,condvec,subjvec] = concat_data(RootDirs,StudyInfo,parms,subj_parms)
  Y = [];
  condvec = [];
  subjvec = [];
  if parms.hiC_flag
    contrast_levels = parms.ncontrasts_orig;
  else
    contrast_levels = [1:parms.ncontrasts];
  end;
  for c=contrast_levels
    if parms.ncontrasts_orig>1
      fprintf('%s: loading data for each subject for contrast level %d...\n',...
        mfilename,c);
    else
      fprintf('%s: loading data for each subject...\n',mfilename);
    end;
    contY = [];
    offset = 0;
    for s=1:parms.nsubs
      [subjY,tmp_condvec] = load_subj_data(RootDirs,StudyInfo(s),parms,subj_parms(s),c);
      contY = cat(2,contY,subjY);
      if c==contrast_levels(1)
        tmp_condvec = tmp_condvec + offset;
        condvec = cat(1,condvec,tmp_condvec);
        offset = max(condvec);
        tmp_subjvec = ones(length(tmp_condvec),1)*s;
        subjvec = cat(1,subjvec,tmp_subjvec);
      end;
    end;
    Y = cat(3,Y,contY);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,subj_parms] = bootstr_resamp_data(subj_inds,data,subj_parms)
  % determine which sensors belong to which subjects
  sens_inds = zeros(size(data.condvec));
  j = 1;
  for i=1:length(subj_inds)
    s = subj_inds(i);
    tmp = find(data.subjvec==s);
    k = j + length(tmp) - 1;
    sens_inds(j:k) = tmp;
    j = k + 1;    
  end;
  sens_inds = sens_inds(sens_inds>0);

  % reorder subject info according to subj_inds
  data.MEG_VisitIDs = data.MEG_VisitIDs(subj_inds);
  subj_parms = subj_parms(subj_inds);

  % reorder data according to sens_inds
  data.Y = data.Y(:,sens_inds,:);
  data.condvec = data.Y(sens_inds);
  data.subjvec = data.Y(sens_inds);
  data.F = data.F(sens_inds,:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y,condvec] = load_subj_data(RootDirs,StudyInfo,parms,subj_parms,c)
  Y = [];
  indir = [RootDirs.proc_meg '/' StudyInfo.proc_meg];
  fname = sprintf('%s/matfiles/%s_avg_data.mat',...
    indir,StudyInfo.RCSE_prefix);
  load(fname);

  % get conditions for this contrast level
  ind_contrast = find(subj_parms.contrasts == subj_parms.unique_contrasts(c));
  cond_numbers = [subj_parms.cond_info.cond_number];
  cond_order = cond_numbers(ind_contrast);
  nconds = length(cond_order);
  
  % extract sensor waveforms from avg_data
  Y = rc_load_ret_avgt(avg_data,cond_order,...
    subj_parms.goodchans,subj_parms.scale_matrix);
  Y = single(Y);

  % vector of stimulus location numbers (for IRLS)
  nsensors = length(subj_parms.goodchans);
  condvec = reshape(ones(nsensors,1)*[1:nconds],[nconds*nsensors,1]);

  % resample to common time samples with spline
  subj_time = avg_data.averages(1).time;
  if parms.ntpoints~=length(subj_time) || any(parms.time~=subj_time)
    nwforms = size(Y,2);
    Yres = zeros(parms.ntpoints,nwforms);
    for w=1:nwforms
      Yres(:,w) = spline(subj_time,Y(:,w),parms.time);
    end;
    Y = Yres;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_waveforms(data,parms,subj_parms)
  % initialize output
  results = init_results(data);

  % calculate inverse
  results.W = calc_inverse(data.F,parms,subj_parms);

  for c=1:results.ncontrasts
    % select data for this contrast level
    tmpY = data.Y(:,:,c);

    % calculate source estimates
    tmpS = (results.W*tmpY')';
    % insert into complete matrix
    j = 1 + (c-1)*parms.nareas;
    k = j + parms.nareas - 1;
    results.S(:,:,c) = tmpS;
    
    % calculate fit
    tmpYfit = (data.F*tmpS')';
    results.Yfit(:,:,c) = tmpYfit;

    % calculate error
    tmpE = tmpY - tmpYfit;
    results.E(:,:,c) = tmpE;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_waveforms_wtd(data,results,parms,subj_parms)
  % initialize output
  results = init_results(data,results);

  % apply IRLS weights to F and Y
  [F_wtd,Y_wtd] = apply_IRLS_weights(data,results);

  % calculate inverse
  results.W = calc_inverse(F_wtd,parms,subj_parms); % weighted forward solution

  for c=1:results.ncontrasts
    % select data for this contrast level
    tmpY = data.Y(:,:,c); % unweighted data
    tmpY_wtd = Y_wtd(:,:,c); % weighted data

    % calculate source estimates
    tmpS = (results.W*tmpY_wtd')'; % weighted data
    % insert into complete matrix
    j = 1 + (c-1)*parms.nareas;
    k = j + parms.nareas - 1;
    results.S(:,:,c) = tmpS;
    
    % calculate fit
    tmpYfit = (data.F*tmpS')'; % unweighted forward solution
    results.Yfit(:,:,c) = tmpYfit;

    % calculate error
    tmpE = tmpY - tmpYfit; % unweighted data
    results.E(:,:,c) = tmpE;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(data,results)
  if ~exist('results','var'), results = []; end;
  results.time = data.time;
  results.ntpoints = data.ntpoints;
  results.areas = data.areas;
  results.nareas = data.nareas;
  results.nwforms = data.nwforms;
    results.contrasts = data.contrasts;
  results.ncontrasts = data.ncontrasts;
  results.condvec = data.condvec;
  results.subjvec = data.subjvec;
  results.nsubjects = data.nsubjects;
  results.MEG_VisitIDs = data.MEG_VisitIDs;
  actual_nareas = length(results.areas);
  if results.nareas>actual_nareas
    results.nwforms = results.nwforms*actual_nareas/results.nareas;
    results.nareas = actual_nareas;
  end;
  results.S = zeros(results.ntpoints,results.nareas,results.ncontrasts);
  results.Yfit = zeros(size(data.Y));
  results.E = zeros(size(data.Y));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F_wtd,Y_wtd] = apply_IRLS_weights(data,results)
  weights = mmil_getfield(results,'weights',[]);
  F_wtd = data.F;
  Y_wtd = data.Y;
  if ~isempty(weights)
    % apply weights to F and Y, using condvec
    num_conds = length(weights);
    for i=1:num_conds
      ind = find(data.condvec==i);
      F_wtd(ind,:) = weights(i)*F_wtd(ind,:);
      Y_wtd(:,ind,:) = weights(i)*Y_wtd(:,ind,:);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_waveforms_IRLS(data,parms,subj_parms)
  % calculate initial source estimates
  fprintf('%s: calculating initial source estimates...\n',mfilename);

  if parms.reweight_init_flag && any([subj_parms.reweight_flag])
    if all([subj_parms.reweight_flag])
      weights = [subj_parms.cond_weights];
    else
      weights = [];
      for s=1:parms.nsubs
        if ~isempty(subj_parms(s).cond_weights)
          weights = [weights subj_parms(s).cond_weights];
        else
          weights = [weights ones(1,length(subj_parms(s).unique_location_conds))];
        end;
      end;
    end;
    results = [];
    results.weights = weights;
    results = calc_waveforms_wtd(data,results,parms,subj_parms);
  else
    results = calc_waveforms(data,parms,subj_parms);
  end;
  
  S_init = results.S;

  % number of stimulus locations x number of subjects
  num_conds = max(data.condvec);

  % choose small value relative to data amplitude
  tiny_s = 1e-6 * max(max(std(data.Y,0,2),[],3));
  % std across conditions, max across contrasts, max across time
  if tiny_s==0, tiny_s = 1; end

  % adjust residuals using "leverage"
  if parms.reweight_leverage_flag
    results.adjfactor = calc_IRLS_leverage(data.F,...
      data.condvec,parms.reweight_leverage_max_flag);
  else
    results.adjfactor = 1;
  end;

  % iteratively reweighted least squares
  for rw_iter=1:parms.reweight_maxiter
    fprintf('%s: reweight iter %d\n',mfilename,rw_iter);
    % calculate average abs(error) for each condition
    results.errvec = [];
    for i=1:num_conds
      ind = find(data.condvec==i);
      % calculate average error for this location/subject
      err = results.E(:,ind,:); % 3rd dim is contrast
      results.errvec = [results.errvec mean(abs(err(:)))];
    end
    clear err;

    % calculate weights from errvec
    [results.weights,results.err_norm] = ...
      calc_IRLS_weights(results.errvec,parms.reweight_factor,...
        results.adjfactor,tiny_s);

    % recalculate source estimates
    fprintf('%s: recalculating source estimates...\n',mfilename);
    results = calc_waveforms_wtd(data,results,parms,subj_parms);

    if rw_iter==1
      results.init_weights = results.weights;
      results.init_errvec = results.errvec;
      results.init_err_norm = results.err_norm;
    end;

    % check for convergence
    S_diff = results.S-S_init;
    S_diff = mean(abs(S_diff(:)))/max(abs(S_init(:)));
    if S_diff<parms.reweight_tol
      fprintf('%s: difference in estimates is less than tolerance\n',mfilename);
      break;
    end;
    S_init = results.S;
  end;

  if rw_iter>=parms.reweight_maxiter
    fprintf('%s: WARNING: reweighted least squares did not converge within %d iterations\n',...
      mfilename,parms.reweight_maxiter);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adjfactor = calc_IRLS_leverage(X,condvec,max_flag);
  fprintf('%s: calculating IRLS leverage...\n',mfilename);
  X = single(X);
  % from Barret and Gray, 1997
  H = X*inv(X'*X)*X'; % this uses a lot of memory (e.g. 15GB)
  % from matlab statrobustfit
  h = min(.9999, sum(H.*H,2)); % this doubles the memory usage (e.g. 30GB)
  tmp_adjfactor = 1 ./ sqrt(1-h);
  clear H;

  % loop over num_conds, get mean or max for each
  num_conds = max(condvec);
  adjfactor = ones(1,num_conds);
  for i=1:num_conds
    if max_flag
      adjfactor(i) = max(tmp_adjfactor(condvec==i));
    else
      adjfactor(i) = mean(tmp_adjfactor(condvec==i));
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [weights,err_norm] = calc_IRLS_weights(errvec,reweight_factor,adjfactor,tiny_s)
  % calculate deviation from minimum error
  err_diff = errvec - min(errvec);

  % adjust for "leverage"
  err_diff = err_diff .* adjfactor;

  % calculate median absolute deviation
  s = max(median(err_diff) / 0.6745,tiny_s);

  % offset and normalize errors
  weights = err_diff/(s*reweight_factor);
  err_norm = err_diff/s;

  % apply Tukey's bisquare weight function
  %   (greater error than median)
  weights = sqrt((abs(weights<1)) .* (1 - weights.^2).^2);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_IRLS_weights(results,parms)
  if isempty(results.err_norm)
    fprintf('%s: WARNING: IRLS err_norm is empty\n',mfilename);
    return;
  end;

  [err_sorted,ind_sort] = sort(results.err_norm);
  weights_sorted = results.weights(ind_sort);

  figure(2); clf;
  hist(results.err_norm)
  title('error histogram');
  xlabel('average error per condition');
  ylabel('number of conditions');
  save_plot(parms,[parms.outstem '_err_hist']);

  figure(3); clf;
  plot(err_sorted,weights_sorted,'r^');
  set(gca,'YLim',[0,1]);
  xlabel('normalized error');
  ylabel('IRLS weights');
  title('weights vs. error');
  save_plot(parms,[parms.outstem '_weights_vs_err']);

  figure(4); clf;
  plot(err_sorted,'go');
  title('normalized error');
  xlabel('stimulus location (sorted by error)');
  ylabel('average error per condition');
  set(gca,'YLim',[0,3]);
  save_plot(parms,[parms.outstem '_norm_err']);

  figure(5); clf;
  plot(weights_sorted,'b*');
  set(gca,'YLim',[0,1]);
  xlabel('stimulus location (sorted by error)');
  ylabel('IRLS weights');
  title('weights');
  save_plot(parms,[parms.outstem '_weights']);

  % calculate average err and weight for each subject
  subj_err_norm = zeros(1,results.nsubjects);
  subj_weights  = zeros(1,results.nsubjects);
  for s=1:results.nsubjects
    ind_subj = find(results.subjvec==s);
    ind_subj_conds = unique(results.condvec(ind_subj));
    subj_weights(s) = mean(results.weights(ind_subj_conds));
    subj_err_norm(s) = mean(results.err_norm(ind_subj_conds));
  end;
  [subj_err_sorted,ind_sort] = sort(subj_err_norm);
  subj_weights_sorted = subj_weights(ind_sort);

  figure(6); clf;
  plot(subj_err_sorted,subj_weights_sorted,'r^');
  set(gca,'YLim',[-0.05,1.05]);

  set(gca,'XLim',[0,max(subj_err_sorted)+0.3]);
  xlabel('average error per subject');
  ylabel('IRLS weights');
  title('subject weights vs. error');
  % mark plot with MEG_VisitIDs
  hold on;
  for s=1:results.nsubjects
    VisitID = regexprep(results.MEG_VisitIDs{s},'_',' ');
    x = subj_err_norm(s) + parms.plot_text_dx;
    y = subj_weights(s) + parms.plot_text_dy;
    text(x,y,VisitID);
  end;
  save_plot(parms,[parms.outstem '_subj_weights_vs_err']);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_waveforms(results,parms)
  plot_parms.outstem = [parms.outstem '_group_wforms'];
  plot_parms.title = 'Multiple Subject Constrained Source Waveforms';
  plot_parms.fontsize = parms.fontsize;
  plot_parms.fontname = parms.fontname;
  plot_parms.ylabel = sprintf('Source Amplitude (%s)',parms.units_wform);
  plot_parms.label_flag = parms.label_flag;
  plot_parms.ylim = parms.ylim_wform;
  plot_parms.legend_loc = parms.legend_loc;
  plot_parms.eps_flag = parms.eps_flag;
  plot_parms.visible_flag = parms.visible_flag;
  plot_parms.fill_err_flag = parms.fill_err_flag;
  plot_parms.time = 1000*parms.time;

  % vary linewidth for varying contrast levels
  if results.ncontrasts > 1
    linewidths = parms.min_linewidth + ...
                        ([1:results.ncontrasts]-1)*(parms.max_linewidth - parms.min_linewidth)/...
                        (results.ncontrasts-1);
  else
    linewidths = parms.linewidth;
  end;

  plot_parms.colors = cell(results.nwforms,1);
  plot_parms.linewidth = zeros(results.nwforms,1);
  plot_parms.condnames = cell(results.nwforms,1);
  k = 1;
  for c=1:results.ncontrasts
    for a=1:results.nareas
      plot_parms.colors{k} = parms.area_colors{a};
      plot_parms.linewidth(k) = linewidths(c);
      plot_parms.condnames{k} = sprintf('%s cont %0.2f',...
        parms.area_names{a},results.contrasts(c));
      k = k + 1;
    end;
  end;

  wforms = reshape(results.S,[results.ntpoints,results.nwforms]);
  if parms.resamp_flag
    plot_parms.wforms_err = cat(3,...
      reshape(results.S_lo,[results.ntpoints,results.nwforms]),...      
      reshape(results.S_hi,[results.ntpoints,results.nwforms]));
  else
    plot_parms.wforms_err = [];
  end;

  plot_args = mmil_parms2args(plot_parms);
  figure(1);
  ts_plot_wforms(wforms,plot_args{:});  

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = analyze_waveforms(results,parms)
  % find peaks, calculate auc, make plots
  tmp_parms = parms;
  tmp_parms.condition_values = unique(results.contrasts);
  tmp_parms.time = 1000*parms.time;
  wforms = results.S;
  args = MMIL_Args(tmp_parms,'rc_analyze_wforms');
  results.analysis = rc_analyze_wforms(wforms,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_plot(parms,outstem)
  if ~parms.visible_flag, set(gcf,'Visible','off'); end;
  print(gcf,'-dtiff',[outstem '.tif']);
  if parms.eps_flag, mmil_printeps(gcf,[outstem '.eps']); end;
  if ~parms.visible_flag, close(gcf); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
