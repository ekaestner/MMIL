function ts_groupavg(studyinfo,varargin)
%function ts_groupavg(studyinfo,varargin)
%
%  Purpose: calculate group averages and t-stats from surface stats
%   in mgh format, resampled to common space of FreeSurfer's icosahedral sphere
%
%   Averages and t-tests will be calculated for each condition and group
%     specified in studyinfo, for the hemispheres specified by 'hemilist'
%   Paired t-tests between conditions or unpaired t-tests between groups
%     will be calculated only if 'cond_contrasts' or 'group_contrasts'
%     parameters are specified
%
% Usage:
%  ts_groupavg(studyinfo, 'key1', value1,...);
%
% Required Parameters:
%   studyinfo: structure with fields containing study info in cell arrays
%     studyinfo.subjnames: subject names (optional)
%     studyinfo.groupnames: group names (optional)
%     studyinfo.datafiles: file names for each subject and condition
%       (cell matrix)
%     studyinfo.condnames: condition names
%       (1 for each column in studyinfo.datafiles)
%     studyinfo.hemis: hemisphere labels ('lh' or 'rh')
%       (1 for each column in studyinfo.datafiles)
%    See ts_load_studyinfo.m for details
%
% Optional Parameters:
%   'prefix': prefix of all output files
%     {default = 'groupavg'}
%   'outdir': output directory
%     {default = pwd} (current working directory)
%   'hemilist': cell array of strings
%     specifying which hemis to include
%     {default = {'lh','rh'}}
%   'cond_contrasts': cell matrix of strings
%     specifying which conditions to contrast with paired t-tests
%     e.g. {'A','B'} or {'A','B';'A','C'}
%     if string 'all' is supplied, all possible combinations of conditions
%     will be contrasted
%     if empty, no condition contrasts will be calculated
%     {default = []}
%   'group_contrasts': cell matrix of strings
%     specifying which groups to contrast with unpaired t-tests
%     e.g. {'A','B'} or {'A','B';'A','C'}
%     if string 'all' is supplied, all possible combinations of groups
%     will be contrasted
%     if empty, no group contrasts will be calculated
%     {default = []}
%   'mbmask_flag': [0|1] whether to mask midbrain & corpus callosum
%     {default = 0}
%   'offset': value subtracted from all data points before calculating means and t-stats
%     For noise-normalized dSPMs, the baseline has an expected value of 1
%     {default = 1}
%   'forceflag': [0|1] whether to regenerate output files if they exist
%     {default = 0}
%   'verbose': [0|1] whether to print status information to the screen
%     {default = 1}
%
% See also: ts_load_studyinfo, fs_groupavg
%
%  Created:  01/13/2008 by Don Hagler
%  Modified: 06/03/2013 by Don Hagler

%% todo: log file?
%% todo: averages for frames from f1 to f2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse arguments
if (~mmil_check_nargs(nargin, 1)), return; end;
parms = mmil_args2parms( varargin, { ...
  'prefix','groupavg',[],...
  'outdir',pwd,[],...
  'forceflag',false,[false true],... 
  'hemilist',{'lh','rh'},[],...
  'cond_contrasts',{},[],...
  'group_contrasts',{},[],...
  'mbmask_flag',false,[false,true],...
  'nverts_ico7',163842,[],...
  'offset',1,[],...
  'subj','fsaverage',[],...
  'verbose',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parms
if ~iscell(parms.hemilist)
  parms.hemilist = {parms.hemilist};
end;
%parms.hemilist = unique(parms.hemilist);
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  if ~ismember(hemi,{'lh','rh'})
    error('entries in parms.hemilist must be ''lh'' or ''rh''');
  end;
end;

if ~exist(parms.outdir,'dir')
  [success,msg] = mkdir(parms.outdir);
  if ~success
    error('failed to create outdir %s: %s',parms.outdir,msg);
  end;
end;

if parms.mbmask_flag
  parms.subjdir = getenv('FREESURFER_HOME');
  if isempty(parms.subjdir)
    error('FREESURFER_HOME not defined as an environment variable (required if mbmask_flag=1)');
  end;
  parms.subjdir = sprintf('%s/subjects',parms.subjdir);
  fullsubj = sprintf('%s/%s',parms.subjdir,parms.subj);
  if ~exist(fullsubj,'dir')
    error('average subject %s not found in %s',parms.subj,parms.subjdir);
  end;
else
  parms.subjdir = [];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check studyinfo
if ~isfield(studyinfo,'datafiles') | isempty(studyinfo.datafiles)
  error('studyinfo.datafiles is empty');
elseif ~iscell(studyinfo.datafiles)
  studyinfo.datafiles = {studyinfo.datafiles};
end;
[nsubs,nconds] = size(studyinfo.datafiles);
if nsubs==1
  error('studyinfo.datafiles has data for only one subject');
end;
if isfield(studyinfo,'groupnames') & ~isempty(studyinfo.groupnames)
  if length(studyinfo.groupnames)~=nsubs
    error('length of studyinfo.groupnames (%d) does not match number of rows (subjects) of studyinfo.datafiles (%d)',...
      length(studyinfo.groupnames),nsubs);
  end;
  groups = unique(studyinfo.groupnames);
  ngroups = length(groups);
else
  groups = [];
  ngroups = 1;
end;
if length(studyinfo.condnames)~=nconds
  error('studyinfo.condnames must have %d entries to match studyinfo.datafiles',...
    nconds);
end;
if length(studyinfo.hemis)~=nconds
  error('studyinfo.hemis must have %d entries to match studyinfo.datafiles',...
    nconds);
end;
fprintf('%s: %d subjects, %d groups, %d conditions\n',...
  mfilename,nsubs,ngroups,nconds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up condition contrasts
if ~isempty(parms.cond_contrasts)
  if ~iscell(parms.cond_contrasts)
    if strcmp(parms.cond_contrasts,'all')
      tmp_conds = unique(studyinfo.condnames);
      if ~iscell(tmp_conds), tmp_conds={tmp_conds}; end;
      tmp_nconds = length(tmp_conds);
      parms.cond_contrasts = {};
      c = 1;
      for j=1:tmp_nconds-1
        for k=j+1:tmp_nconds
          parms.cond_contrasts{c,1} = tmp_conds{j};
          parms.cond_contrasts{c,2} = tmp_conds{k};
          c = c+1;
        end;
      end;
    else
      error('cond_contrasts is not a cell array and is not ''all''');
    end;
  end;
  if size(parms.cond_contrasts,2)~=2
    error('size of second dimension of cond_contrasts must be 2 (see help %s)',...
      mfilename);
  end;    
  ncond_contrasts = size(parms.cond_contrasts,1);
  if parms.verbose
    if ncond_contrasts==1
      fprintf('%s: %d condition contrast specified\n',mfilename,ncond_contrasts);
    else
      fprintf('%s: %d condition contrasts specified\n',mfilename,ncond_contrasts);
    end;
  end;
  % check that they are specified correctly
  for c=1:ncond_contrasts
    condA = parms.cond_contrasts{c,1};
    condB = parms.cond_contrasts{c,2};
    if strcmp(condA,condB)
      error('condition %s specified twice in condition contrast %d',...
        condA,c);
    end;
    for h=1:length(parms.hemilist)
      hemi = parms.hemilist{h};
      cA_ind = find(strcmp(studyinfo.condnames,condA) & strcmp(studyinfo.hemis,hemi));
      if isempty(cA_ind)
        error('condition %s, hemi %s from condition contrast %d not found in studyinfo',...
          condA,hemi,c);
      elseif length(cA_ind)>1
        error('condition %s, hemi %s from condition contrast %d duplicated in studyinfo',...
          condA,hemi,c);
      end;
      cB_ind = find(strcmp(studyinfo.condnames,condB) & strcmp(studyinfo.hemis,hemi));
      if isempty(cB_ind)
        error('condition %s, hemi %s from condition contrast %d not found in studyinfo',...
          condB,hemi,c);
      elseif length(cB_ind)>1
        error('condition %s, hemi %s from condition contrast %d duplicated in studyinfo',...
          condA,hemi,c);
      end;
    end;
  end;
else
  ncond_contrasts = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up group contrasts
if ngroups<=1, parms.group_contrasts=[]; end;
if ~isempty(parms.group_contrasts)
  if ~iscell(parms.group_contrasts)
    if strcmp(parms.group_contrasts,'all')
      parms.group_contrasts = {};
      c = 1;
      for j=1:ngroups
        for k=j+1:ngroups
          parms.group_contrasts{c,1} = groups{j};
          parms.group_contrasts{c,2} = groups{k};
          c = c+1;
        end;
      end;
    else
      error('group_contrasts is not a cell array and is not ''all''');
    end;
  end;
  % check that size is correct
  if size(parms.group_contrasts,2)~=2
    error('size of second dimension of group_contrasts must be 2 (see help %s)',...
      mfilename);
  end;    
  ngroup_contrasts = size(parms.group_contrasts,1);
  if parms.verbose
    if ngroup_contrasts==1
      fprintf('%s: %d group contrast specified\n',mfilename,ngroup_contrasts);
    else
      fprintf('%s: %d group contrasts specified\n',mfilename,ngroup_contrasts);
    end;
  end;
  % check that they are specified correctly
  for c=1:ngroup_contrasts
    groupA = parms.group_contrasts{c,1};
    groupB = parms.group_contrasts{c,2};
    if strcmp(groupA,groupB)
      error('group %s specified twice in group contrast %d',...
        groupA,c);
    end;
    gA_ind = find(strcmp(studyinfo.groupnames,groupA));
    if isempty(gA_ind)
      error('group %s from group contrast %d not found in studyinfo',...
        groupA,c);
    end;
    gB_ind = find(strcmp(studyinfo.groupnames,groupB));
    if isempty(gB_ind)
      error('group %s from group contrast %d not found in studyinfo',...
        groupB,c);
    end;
  end;
else
  ngroup_contrasts = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check data files
nframes = [];
verts = [];
datasize = [];
for i=1:nsubs
  for j=1:nconds
    fname = studyinfo.datafiles{i,j};
    if ~isempty(fname) & ~exist(fname,'file')
      error('data file %s not found',fname);
    else
      try
        [vol,M,mr_parms,volsz] = fs_load_mgh(fname,[],[],1);
      catch
        error('error loading %s as mgh file',fname);
      end;
      if isempty(datasize)
        datasize = volsz;
        nframes = datasize(4);
        nverts = datasize(1);
      elseif any(datasize~=volsz)
        error('file %s has mismatched data size ([%s] instead of [%s])',...
          fname,...
          regexprep(num2str(volsz),'\s+',','),...
          regexprep(num2str(datasize),'\s+',','));
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate averages and one-sample ttests for each condition and group
for c=1:nconds
  hemi = studyinfo.hemis{c};
  condname = studyinfo.condnames{c};
  for g=1:ngroups
    if isempty(groups)
      groupname = [];
      s_ind = [1:nsubs];
    else
      groupname = groups{g};
      s_ind = find(strcmp(groupname,studyinfo.groupnames));
    end;
    % construct output filenames, check if they exist
    tmp_prefix = sprintf('%s/%s',parms.outdir,parms.prefix);
    if ~isempty(groupname)
      tmp_prefix = sprintf('%s_group_%s',tmp_prefix,groupname);
    end;
    tmp_prefix = sprintf('%s_cond_%s',tmp_prefix,condname);
    outname_mean = sprintf('%s_mean-%s.mgh',tmp_prefix,hemi);
    outname_stdv = sprintf('%s_stdv-%s.mgh',tmp_prefix,hemi);
    outname_tstat = sprintf('%s_tstat-%s.mgh',tmp_prefix,hemi);
    outname_pval = sprintf('%s_pval-%s.mgh',tmp_prefix,hemi);
    outname_n = sprintf('%s_n-%s.mat',tmp_prefix,hemi);
    if ~exist(outname_mean,'file') | ~exist(outname_stdv,'file') |...
       ~exist(outname_tstat,'file') | ~exist(outname_pval,'file') |...
       ~exist(outname_n,'file') | parms.forceflag
      if parms.verbose
        if ~isempty(groupname)
          fprintf('%s: calculating group averages and ttests for group %s, condition %s, %s...\n',...
            mfilename,groupname,condname,hemi);
        else
          fprintf('%s: calculating group averages and ttests for condition %s %s...\n',...
            mfilename,condname,hemi);
        end;
      end;
      % initialize output
      outvol_mean = zeros(datasize);
      outvol_stdv = zeros(datasize);
      outvol_tstat = zeros(datasize);
      outvol_pval = zeros(datasize);
      for f=1:nframes
        if parms.verbose
          fprintf('%s: frame %d\n',mfilename,f);
        end;
        results = fs_groupavg(studyinfo.datafiles(s_ind,c),...
          'gnamelist',groupname,...
          'frames',f,...
          'offset',parms.offset,...
          'verbose_flag',parms.verbose);
        % even though fs_groupavg can do between group comparisons,
        %   we will do one group at a time here and do the unpaired ttests later
        %   (so that only one group's output data needs to be in memory at one time)
        if isempty(results)
          error('group average for group %s, condition %s failed',groupname,condname);
        end;
        % insert results for this frame into outvols
        outvol_mean(:,:,:,f) = results.groups.mean;
        outvol_stdv(:,:,:,f) = results.groups.stdv;
        outvol_tstat(:,:,:,f) = results.groups.tstats;
        outvol_pval(:,:,:,f) = results.groups.pvals;
      end;
      n = results.groups.n;
      df = results.groups.df;
      % save output
      if parms.verbose
        fprintf('%s: saving output...\n',mfilename);
      end;
      save_masked_mgh(outvol_mean,outname_mean,...
        hemi,parms);
      save_masked_mgh(outvol_stdv,outname_stdv,...
        hemi,parms);
      save_masked_mgh(outvol_tstat,outname_tstat,...
        hemi,parms);
      save_masked_mgh(outvol_pval,outname_pval,...
        hemi,parms);
      save(outname_n,'n','df');
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate average differences and paired ttests for each specified condition contrast
for c=1:ncond_contrasts
  condA = parms.cond_contrasts{c,1};
  condB = parms.cond_contrasts{c,2};
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    cA_ind = find(strcmp(studyinfo.condnames,condA) & strcmp(studyinfo.hemis,hemi));
    cB_ind = find(strcmp(studyinfo.condnames,condB) & strcmp(studyinfo.hemis,hemi));
    for g=1:ngroups
      if isempty(groups)
        groupname = [];
        s_ind = [1:nsubs];
      else
        groupname = groups{g};
        s_ind = find(strcmp(groupname,studyinfo.groupnames));
      end;
      % construct output filenames, check if they exist
      tmp_prefix = sprintf('%s/%s',parms.outdir,parms.prefix);
      if ~isempty(groupname)
        tmp_prefix = sprintf('%s_group_%s',tmp_prefix,groupname);
      end;
      tmp_prefix = sprintf('%s_cond_%s_VS_cond_%s',tmp_prefix,condA,condB);
      outname_mean = sprintf('%s_mean-%s.mgh',tmp_prefix,hemi);
      outname_stdv = sprintf('%s_stdv-%s.mgh',tmp_prefix,hemi);
      outname_tstat = sprintf('%s_tstat-%s.mgh',tmp_prefix,hemi);
      outname_pval = sprintf('%s_pval-%s.mgh',tmp_prefix,hemi);
      outname_n = sprintf('%s_n-%s.mat',tmp_prefix,hemi);
      if ~exist(outname_mean,'file') | ~exist(outname_stdv,'file') |...
         ~exist(outname_tstat,'file') | ~exist(outname_pval,'file') |...
         ~exist(outname_n,'file') | parms.forceflag
        if parms.verbose
          if ~isempty(groupname)
            fprintf('%s: calculating group averages and ttests for group %s, condition %s VS %s, %s...\n',...
              mfilename,groupname,condA,condB,hemi);
          else
            fprintf('%s: calculating group averages and ttests for condition %s VS %s, %s...\n',...
              mfilename,condA,condB,hemi);
          end;
        end;
        % initialize output
        outvol_mean = zeros(datasize);
        outvol_stdv = zeros(datasize);
        outvol_tstat = zeros(datasize);
        outvol_pval = zeros(datasize);
        for f=1:nframes
          if parms.verbose
            fprintf('%s: frame %d\n',mfilename,f);
          end;
          results = fs_groupavg(studyinfo.datafiles(s_ind,[cA_ind,cB_ind]),...
            'gnamelist',groupname,...
            'condition_weights',[1 -1],...
            'frames',f,...
            'offset',parms.offset,...
            'verbose_flag',parms.verbose);
          if isempty(results)
            error('group average for group %s, condition %s VS %s failed',groupname,condA,condB);
          end;
          % insert results for this frame into outvols
          outvol_mean(:,:,:,f) = results.groups.mean;
          outvol_stdv(:,:,:,f) = results.groups.stdv;
          outvol_tstat(:,:,:,f) = results.groups.tstats;
          outvol_pval(:,:,:,f) = results.groups.pvals;
        end;
        n = results.groups.n;
        df = results.groups.df;
        % save output
        if parms.verbose
          fprintf('%s: saving output...\n',mfilename);
        end;
        save_masked_mgh(outvol_mean,outname_mean,...
          hemi,parms);
        save_masked_mgh(outvol_stdv,outname_stdv,...
          hemi,parms);
        save_masked_mgh(outvol_tstat,outname_tstat,...
          hemi,parms);
        save_masked_mgh(outvol_pval,outname_pval,...
          hemi,parms);
        save(outname_n,'n','df');
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% difference of avgs and unpaired ttest for contrasts between groups
for g=1:ngroup_contrasts
  groupA = parms.group_contrasts{g,1};
  groupB = parms.group_contrasts{g,2};
  for c=1:nconds
    hemi = studyinfo.hemis{c};
    condname = studyinfo.condnames{c};

    % construct input filenames for groupA, check if they exist
    tmp_prefix = sprintf('%s/%s',parms.outdir,parms.prefix);
    tmp_prefix = sprintf('%s_group_%s',tmp_prefix,groupA);
    tmp_prefix = sprintf('%s_cond_%s',tmp_prefix,condname);
    innameA_mean = sprintf('%s_mean-%s.mgh',tmp_prefix,hemi);
    innameA_stdv = sprintf('%s_stdv-%s.mgh',tmp_prefix,hemi);
    innameA_n = sprintf('%s_n-%s.mat',tmp_prefix,hemi);
    if ~exist(innameA_mean,'file')
      error('file %s not found',innameA_mean);
    end;
    if ~exist(innameA_stdv,'file')
      error('file %s not found',innameA_stdv);
    end;
    if ~exist(innameA_n,'file')
      error('file %s not found',innameA_n);
    end;

    % construct input filenames for groupB, check if they exist
    tmp_prefix = sprintf('%s/%s',parms.outdir,parms.prefix);
    tmp_prefix = sprintf('%s_group_%s',tmp_prefix,groupB);
    tmp_prefix = sprintf('%s_cond_%s',tmp_prefix,condname);
    innameB_mean = sprintf('%s_mean-%s.mgh',tmp_prefix,hemi);
    innameB_stdv = sprintf('%s_stdv-%s.mgh',tmp_prefix,hemi);
    innameB_n = sprintf('%s_n-%s.mat',tmp_prefix,hemi);
    if ~exist(innameB_mean,'file')
      error('file %s not found',innameB_mean);
    end;
    if ~exist(innameB_stdv,'file')
      error('file %s not found',innameB_stdv);
    end;
    if ~exist(innameB_n,'file')
      error('file %s not found',innameB_n);
    end;

    % construct output filenames, check if they exist
    tmp_prefix = sprintf('%s/%s',parms.outdir,parms.prefix);
    tmp_prefix = sprintf('%s_group_%s_VS_group_%s',tmp_prefix,groupA,groupB);
    tmp_prefix = sprintf('%s_cond_%s',tmp_prefix,condname);
    outname_mean = sprintf('%s_mean-%s.mgh',tmp_prefix,hemi);
    outname_stdv = sprintf('%s_stdv-%s.mgh',tmp_prefix,hemi);
    outname_tstat = sprintf('%s_tstat-%s.mgh',tmp_prefix,hemi);
    outname_pval = sprintf('%s_pval-%s.mgh',tmp_prefix,hemi);
    outname_n = sprintf('%s_n-%s.mgh',tmp_prefix,hemi);
    if ~exist(outname_mean,'file') | ~exist(outname_stdv,'file') |...
       ~exist(outname_tstat,'file') | ~exist(outname_pval,'file') |...
       ~exist(outname_n,'file') | parms.forceflag
      if parms.verbose
        fprintf('%s: calculating group averages and ttests for group %s VS %s, condition %s, %s...\n',...
          mfilename,groupA,groupB,condname,hemi);
      end;
      % load precalculated group means and standard deviations
      volA_mean = fs_load_mgh(innameA_mean);
      volA_stdv = fs_load_mgh(innameA_stdv);
      volB_mean = fs_load_mgh(innameB_mean);
      volB_stdv = fs_load_mgh(innameB_stdv);
      % load num subjects
      clear n; load(innameA_n); nA = n;
      clear n; load(innameB_n); nB = n;
      n = nA + nB;
      df = n-2;
      % calculate difference of means
      outvol_mean = volA_mean - volB_mean;
      % calculate pooled standard deviations
      outvol_stdv = sqrt(((nA-1)*volA_stdv.^2+...
                          (nB-1)*volB_stdv.^2)/(n-2));
      % calculate tstats
      outvol_tstat = outvol_mean ./ (outvol_stdv ./ sqrt(eps+n) + eps);
      % calculate pvals
      outvol_pval = -log10(2*tcdf(-abs(outvol_tstat),df)); % 2-tailed
      % save output
      if parms.verbose
        fprintf('%s: saving output...\n',mfilename);
      end;
      save_masked_mgh(outvol_mean,outname_mean,...
        hemi,parms);
      save_masked_mgh(outvol_stdv,outname_stdv,...
        hemi,parms);
      save_masked_mgh(outvol_tstat,outname_tstat,...
        hemi,parms);
      save_masked_mgh(outvol_pval,outname_pval,...
        hemi,parms);
      save(outname_n,'n','df');
    end;
  end;
end;

return;

% function to save to mgh with masking
function save_masked_mgh(vol,fname,hemi,parms)
  if parms.mbmask_flag
    if size(vol,1)~=parms.nverts_ico7
      if parms.verbose
        fprintf('%s: WARNING: skipping mbmask because not ico7',...
          mfilename);
      end;
    else
      vec = squeeze(vol);
      vec = fs_mask_surfstats_with_aparc(vec,parms.subj,hemi,parms.subjdir);
      vol = reshape(vec,[size(vec,1),1,1,size(vec,2)]);
    end;
  end;
  fs_save_mgh(vol,fname);
return;
