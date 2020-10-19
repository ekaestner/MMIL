function MMIL_SurfGroupGLM_BOLD(RootDirs,varargin)
%function MMIL_SurfGroupGLM_BOLD(RootDirs,[options])
%
% Required Input:
%  RootDirs:
%    a struct which must contain the following fields:
%         proc_bold, fsurf
%    and may contain the following fields:
%         orig, raw, proc, prob_dti, proc_bold
%         fsurf, fsico, fsclean, long
%    these specify the locations of data
%
% Optional Parameters:
%  'infix': if empty, will look for files like 'BOLD1.mgz'
%     otherwise, input file will be sprintf('BOLD%d_%s.mgz',snum,infix)
%     example infix = 'corr_resBOLD'
%     {default = []}
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%      must contain these fields:
%        SubjID
%        StudyDate
%      may contain these fields
%        proc_bold
%        fsurf
%      if proc_bold and fsurf are unspecified, will look for Containers
%        with SubjID and StudyDate
%       (will choose first one if more than one)
%      if empty, use all subjects found in RootDirs.proc_bold and RootDirs.fsurf
%     {default = []}
%  'ProjID': project name inserted into output path
%     {default = []}
%  'outdir': output directory
%    {default = [RootDirs.home '/MetaData/' ProjID '/BOLD_SurfGroupGLM']}
%  'outstem': output file stem
%    {default = 'BOLD_SurfGroupGLM'}
%  'fname_regressors': full path name of csv file containing subject IDs
%     and regressor variables in additional columns
%     First column should have subject IDs, additional columns should contain
%       regressor values (e.g. group, test scores, etc.)
%     First row should be column headers, number of additional rows
%       is number of subjects included in analysis
%  'groupcols': which columns in csv file correspond to group info
%     These columns should have values of 0 or 1 (in group or not)
%     e.g. if columns 2 and 3 are group info, then groupcols=[2 3]
%    {default = []}
%  'aparc_mask_flag': [0|1] whether to mask non-cortical regions using aparc
%    {default = 0}
%  'multirunflag': [0|1] whether to combine across conditions
%    {default = 1}
%  'coef_frame' - frame number for GLM coefficient in 3dDeconv results
%     {default = 2}
%  'forceflag' - [0|1] whether to run calculations even if output log file exist
%    {default = 0}
%
% Created:  06/22/09 by Don Hagler
% Last Mod: 12/31/12 by Don Hagler
%

%% todo: accept ProjID as required input, RootDirs as optional
%% todo: accept QC_BOLD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms(varargin, { ...
  'infix',[],[],...
  'StudyInfo',[],[],...
  'ProjID',[],[],...
  'fname_regressors',[],[],...
  'groupcols',[],[],...
  'outdir',[],[],...
  'outstem','BOLD_SurfGroupGLM',[],...
  'aparc_mask_flag',false,[false true],...
  'coef_frame',2,[1 Inf],...
  'forceflag',false,[false true],...
...
  'smoothsteps',0,[0 Inf],...
  'sphsmoothsteps',10,[0 Inf],...
  'multirunflag',true,[false true],... % whether to combine across conditions (runs)
...
  'statlist',{'mean','stdv','tstat','log10_pval'},{'mean','stdv','tstat','log10_pval'},...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'required_containers',{'proc_bold','fsurf'},[],...
  'qcflag',true,[false true],...
  'modality','MRI',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter StudyInfo or generate StudyInfo struct if none supplied
args = MMIL_Args(parms,'MMIL_Check_StudyInfo');
[StudyInfo,RootDirs] = MMIL_Check_StudyInfo(parms.StudyInfo,RootDirs,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

% set outputdir
if isempty(parms.outdir)
  parms.outdir = [RootDirs.home '/MetaData/' parms.ProjID '/BOLD_SurfGroupGLM'];
end;
if ~exist(parms.outdir,'dir')
  [s,m]=mkdir(parms.outdir);
  if ~s, error('failed to create output dir %s:\n%s',parms.outdir,m); end;
end;

% read regressor file if supplied
if ~isempty(parms.fname_regressors)
  if ~exist(parms.fname_regressors,'file')
    error('file %s not found',parms.fname_regressors);
  end;
  % read fname_regressors
  try
    input_info = mmil_readtext(parms.fname_regressors);
    SubjIDs = {input_info{2:end,1}};
  catch
    error('failed to get info from regressors file %s -- check format',...
      parms.fname_regressors);
  end;
  regnames = input_info(1,2:end);
  try
    regressors = cell2mat(input_info(2:end,2:end));
  catch
    error('list of regressors in file %s must all be integers',...
      parms.fname_regressors);
  end

  reg_subjIDs =  {input_info{2:end,1}};
else
  regressors = [];
end;
orig_regressors = regressors;
if isempty(regressors), regnames=[]; end;
orig_regnames = regnames;


% find average subject for parcellation-based masking of non-cortical regions
if parms.aparc_mask_flag
  subjdir = getenv('FREESURFER_HOME');
  if isempty(subjdir)
    error('FREESURFER_HOME not defined as an environment variable (required for aparc_mask_flag)');
  end;
  subjdir = sprintf('%s/subjects',subjdir);
  subj = 'fsaverage';
  fullsubj = sprintf('%s/%s',subjdir,subj);
  if ~exist(fullsubj,'dir')
    error('average subject %s not found (required for aparc_mask_flag)',subj);
  end;
end;

fname_log = sprintf('%s/%s.log',parms.outdir,parms.outstem);
if ~exist(fname_log,'file') | parms.forceflag
  % create log file
  flog = fopen(fname_log,'wt');
  if flog==-1
    error('failed to open %s for writing',fname_log);
  end;

  % run glm calculations
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    % get names of surface stats files
    fprintf('%s: getting %s surface stats file names...\n',...
      mfilename,hemi);
    tic;
    [fname_surfstats,SubjIDs] = BOLD_MMIL_Get_SurfStats(RootDirs,...
      'StudyInfo',StudyInfo,...
      'infix',parms.infix,...
      'smoothsteps',parms.smoothsteps,...
      'sphsmoothsteps',parms.sphsmoothsteps,...
      'hemi',hemi);
    toc;
    if isempty(fname_surfstats)
      fprintf('%s: WARNING: %s surfstats files for specified subjects not found\n',...
        mfilename,hemi);
      fprintf(flog,'WARNING: %s surfstats files for specified subjects not found\n',...
        hemi);
      continue;
    end;

    % match SubjIDs and regressors
    [c,ia,ib] = intersect(SubjIDs,reg_subjIDs);
    fname_surfstats = fname_surfstats(ia,:);
    [nsubs,nconds] = size(fname_surfstats);
    if ~isempty(regressors)
      regressors = orig_regressors(ib,:);
    end;
    nregs = size(regressors,2);

    if ~isempty(parms.groupcols)
      % reduce regressors by 1
      ind_newregs = setdiff([1:nregs],parms.groupcols(1)-1);
      regressors = regressors(:,ind_newregs);
      regnames = orig_regnames(ind_newregs);
      ngroups = length(parms.groupcols);
      nvals = nregs-ngroups;
    else
      ngroups = 1;
      nvals = nregs;
    end;

    % create contrast vectors
    ncols = (nconds-1) + (ngroups-1) + nvals + 1;
    contrast_names = [];  contrast_vectors = [];
    k = 1; c = 1;

    % baseline/intercept is last column (but first contrast)
    tmp = zeros(1,ncols);
    tmp(ncols) = 1;
    contrast_vectors{k} = tmp;
    contrast_names{k} = 'baseline';
    k = k + 1;
    
    if nconds>1 & ~parms.multirunflag
      % individual conditions
      for i=2:nconds % if more than one cond, reduce by 1
        tmp = contrast_vectors{1};
        tmp(c) = 1;
        contrast_vectors{k} = tmp;
        contrast_names{k} = sprintf('cond%d',i);
        k = k + 1;
        c = c + 1;
      end;
    else
      c = nconds;
    end;

    % group variables and other regressors
    for i=1:nregs
      if i==parms.groupcols(1)-1, continue; end;
      tmp = zeros(1,ncols);
      tmp(c) = 1;
      if ismember(i,parms.groupcols-1)
        tmp(ncols) = 1;
      end;
      contrast_vectors{k} = tmp;
      contrast_names{k} = orig_regnames{i};
      k = k + 1;
      c = c + 1;
    end;

    % between group contrasts
    if ngroups>1
      for i=1:ngroups
        for j=i+1:ngroups
          tmp = zeros(1,ncols);
          if i~=1
            tmp(nconds-1+parms.groupcols(i)-2) = 1;
          end;
          tmp(nconds-1+parms.groupcols(j)-2) = -1;
          contrast_vectors{end+1} = tmp;
          contrast_names{end+1} = [orig_regnames{parms.groupcols(i)-1} '_' orig_regnames{parms.groupcols(j)-1}];
        end
      end
    end

    %% todo: between condition contrasts

    results = fs_groupavg_glm(fname_surfstats,...
      'regnames',regnames,...
      'regressors',regressors,...
      'groupregs',parms.groupcols-1,...
      'contrast_vectors',contrast_vectors,...
      'frame',parms.coef_frame,...
      'contrast_names',contrast_names);

    if isempty(results)
      fprintf('%s: WARNING: group glm for %s failed\n',...
        mfilename,hemi);
      fprintf(flog,'WARNING: group glm for %s failed\n',...
        hemi);
      continue;
    end;

    for i=1:length(results.contrasts)
      name = results.contrasts(i).name;
      fprintf(flog,'Contrast %d %s: n=%d, df=%d   (hemi %s)\n',...
        i,name,results.contrasts(i).n,results.contrasts(i).dof,hemi);
    end;
    fprintf(flog,'\n\n');

    for i=1:nsubs
      SubjID = SubjIDs{i};
      for j=1:nconds
        fname = fname_surfstats{i,j};
        if isempty(fname), continue; end;
        fprintf(flog,'subject %s, fname %s\n',...
          SubjID,fname);
      end;
    end;

    for i=1:length(results.contrasts)
      name = results.contrasts(i).name;
      for r=1:length(parms.statlist)
        statname = parms.statlist{r};
        if ~isfield(results.contrasts(i),statname)
          fprintf('%s: WARNING: results struct is missing field %s\n',...
            mfilename,statname);
          continue;
        end;
        vec = mmil_colvec(results.contrasts(i).(statname));
        if ~isempty(vec)
          fname = sprintf('%s/%s_%s_%s-%s.mgh',...
            parms.outdir,parms.outstem,name,statname,hemi);
          if parms.aparc_mask_flag
            vec = fs_mask_surfstats_with_aparc(vec,subj,hemi,subjdir);
          end;
          fs_save_mgh(vec,fname);
        end;
      end;
    end;
  end;
  fclose(flog);
end;

