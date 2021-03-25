function MMIL_SurfGroupAvg_BOLD(RootDirs,varargin)
%function MMIL_SurfGroupAvg_BOLD(RootDirs,[options])
%
% Required Input:
%  RootDirs:
%    a struct which must contain the following fields:
%         proc_bold, fsurf
%    and may contain the following fields:
%         orig, raw, proc, proc_dti, proc_bold
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
%    {default = [RootDirs.home '/MetaData/' ProjID '/BOLD_SurfGroupAvg']}
%  'outstem': output file stem
%    {default = 'BOLD_SurfGroupAvg'}
%  'fname_groups': full path name of csv file containing subject IDs
%     and group flags in additional columns
%     First column should have subject IDs, additional columns should contain
%       group info
%     First row should be column headers, number of additional rows
%       is number of subjects included in analysis
%  'aparc_mask_flag': [0|1] whether to mask non-cortical regions using aparc
%    {default = 0}
%  'coef_frame' - frame number for GLM coefficient in 3dDeconv results
%     {default = 2}
%  'forceflag' - [0|1] whether to run calculations even if output log file exist
%    {default = 0}
%
% Multiple conditions will be treated as runs and averaged together
%
% Created:  06/25/09 by Don Hagler
% Last Mod: 12/31/12 by Don Hagler
%

%% todo: accept ProjID as required input, RootDirs as optional
%% todo: accept QC_BOLD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin, 1)), return; end;
parms = mmil_args2parms(varargin, { ...
  'infix',[],[],...
  'StudyInfo',[],[],...
  'ProjID',[],[],...
  'fname_groups',[],[],...
  'outdir',[],[],...
  'outstem','BOLD_SurfGroupAvg',[],...
  'aparc_mask_flag',false,[false true],...
  'coef_frame',2,[1 Inf],...
  'forceflag',false,[false true],...
...
  'smoothsteps',0,[0 Inf],...
  'sphsmoothsteps',10,[0 Inf],...
...
  'statlist',{'mean','stdv','tstats','pvals'},{'mean','stdv','tstats','pvals'},...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'required_containers',{'proc_bold','fsurf'},[],...
  'modality','MRI',[],...
  'qcflag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter StudyInfo or generate StudyInfo struct if none supplied
args = MMIL_Args(parms,'MMIL_Check_StudyInfo');
[StudyInfo,RootDirs] = MMIL_Check_StudyInfo(parms.StudyInfo,RootDirs,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;

% set outputdir
if isempty(parms.outdir)
  parms.outdir = [RootDirs.home '/MetaData/' parms.ProjID '/' parms.outstem];
end;
if ~exist(parms.outdir,'dir')
  [s,m]=mkdir(parms.outdir);
  if ~s, error('failed to create output dir %s:\n%s',parms.outdir,m); end;
end;

% read regressor file if supplied
if ~isempty(parms.fname_groups)
  if ~exist(parms.fname_groups,'file')
    error('file %s not found',parms.fname_groups);
  end;
  % read fname_groups
  try
    input_info = mmil_readtext(parms.fname_groups);
    SubjIDs = {input_info{2:end,1}};
  catch
    error('failed to get info from groups file %s -- check format',...
      parms.fname_groups);
  end;
  groupnames = input_info(1,2:end);
  try
    groupflags = cell2mat(input_info(2:end,2:end));
  catch
    error('group info in file %s must all be integers',...
      parms.fname_groups);
  end

  reg_subjIDs =  {input_info{2:end,1}};
else
  groupflags = [];
end;
orig_groupflags = groupflags;
if isempty(groupflags)
  groupnames=[];
  gnamelist=[];
else
  [nsubs,ngroups] = size(groupflags);
  % make list of group names per subject
  gnamelist = cell(nsubs,1);
  for i=1:nsubs
    ind = find(squeeze(groupflags(i,:)));
    gnamelist{i} = groupnames{ind(1)};
  end;
end;
orig_gnamelist = gnamelist;

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

  % run avg calculations
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

    % match subjIDs in SubjIDs and groupflags
    [c,ia,ib] = intersect(SubjIDs,reg_subjIDs);
    fname_surfstats = fname_surfstats(ia,:);
    [nsubs,nconds] = size(fname_surfstats);
    if ~isempty(groupflags)
      groupflags = orig_groupflags(ib,:);
      gnamelist = orig_gnamelist(ib);
    end;

    condition_weights = ones(nconds,1)/nconds;
    results = fs_groupavg(fname_surfstats,...
      'gnamelist',gnamelist,...
      'frames',parms.coef_frame,...
      'verbose_flag',0,...
      'condition_weights',condition_weights);

    if isempty(results)
      fprintf('%s: WARNING: group avg for %s failed\n',...
        mfilename,hemi);
      fprintf(flog,'WARNING: group avg for %s failed\n',...
        hemi);
      continue;
    end;

    if isempty(results.contrasts), results.contrasts = results.groups; end;

    for i=1:length(results.contrasts)
      name = results.contrasts(i).name;
      fprintf(flog,'Contrast %d %s: n=%d, df=%d   (hemi %s)\n',...
        i,name,results.contrasts(i).n,results.contrasts(i).df,hemi);
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
          fprintf('%s: WARNING: results struct is missing field %s\n',mfilename,statname);
          continue;
          %error('results missing field %s',statname);
        end;
        vec = mmil_colvec(getfield(results.contrasts(i),statname));
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

