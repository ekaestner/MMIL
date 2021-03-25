function results = fs_groupavg(fnamelist,varargin)
%function results = fs_groupavg(fnamelist,[options])
%
% Usage:
%  results = fs_groupavg(fnamelist,'key1', value1,...);
%
% Required Input:
%   fnamelist: cell array of input full path file names
%     All input files should be mgh format surface or volume files
%     registered to common space (e.g. spherical surface)
%     
%     fnamelist can also be a cell matrix of file names
%       Each row should correspond to a single subject
%       Each column should correspond to a condition
%     Specify a linear combination of these conditions with
%       condition_weights
%     For a paired t-test between two conditions, supply
%       fnamelist with two columns and condition_weights = [1 -1]
%
%  Optional Input:
%   'gnamelist' - cell array of group names, one for each row (subject)
%      of fnamelist
%      If ommitted, will calculate average of all files as single group
%   'frames' - vector of frame numbers (e.g. multiple time points)
%     If empty, will loop over all frames
%     {default = []}
%   'condition_weights: vector of linear weights for each
%     condition (each column of fnamelist)
%     {default = [1 1 1 ...]}
%   'offset' - value subtracted from all data points before calculating
%     means and t-stats
%     {default = 0}
%   'stats_flag' - [0|1] whether to calculate statistics other than means
%     (i.e. stdv, tstats, pvals)
%     {default = 1}
%   'contrasts_flag' - [0|1] whether to calculate contrasts between groups
%     {default = 1}
%   'verbose_flag' - [0|1] whether to print status messages
%     {default = 1}
%
% Output:
%   results: structure containing group averages, standard deviations
%     t-stats, and significance maps (-log10(p))
%
%   Separate maps will be generated for each group
%     and for each pairwise comparison between groups (if contrasts_flag=1).
%
% Created:  08/16/07 by Don Hagler
% Last Mod: 11/21/11 by Don Hagler
%

%% todo: allow user to specify which contrasts to do

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'gnamelist',[],[],...
  'frames',[],[],...
  'condition_weights',[],[],...
  'offset',0,[-Inf,Inf],...
  'verbose_flag',true,[false true],...
  'contrasts_flag',true,[false true],...
  'stats_flag',true,[false true],...
});
results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

if isempty(fnamelist)
  error('input file name list is empty');
elseif ~iscell(fnamelist)
  fnamelist = {fnamelist};
end;
[nsubs,nconds] = size(fnamelist);
if nsubs==1 && nconds==1
  error('input file name list contains only one file');
elseif nsubs==1
  fnamelist = reshape(fnamelist,[nconds,nsubs]);
  [nsubs,nconds] = size(fnamelist);
elseif nsubs<nconds
  if parms.verbose_flag
    fprintf('%s: WARNING: fnamelist has more conditions (# columns = %d) than subjects (# rows = %d)\n',...
      nconds,nsubs,mfilename);
  end;
end;

if isempty(parms.gnamelist)
  parms.gnamelist = {};
  for s=1:nsubs
    parms.gnamelist{s} = 'group';
  end;
else
  if ~iscell(parms.gnamelist)
    parms.gnamelist = {parms.gnamelist};
  end;
  if length(parms.gnamelist)==1
    for s=1:nsubs
      parms.gnamelist{s} = parms.gnamelist{1};
    end;
  end;
  if length(parms.gnamelist)~=nsubs
    error('length of gnamelist (%d) does not match number of rows (subjects) of fnamelist (%d)',...
      length(parms.gnamelist),nsubs);
  end;
end;

if isempty(parms.condition_weights)
  parms.condition_weights = ones(nconds,1);
end;
parms.condition_weights = reshape(parms.condition_weights,[prod(size(parms.condition_weights)),1]);
if size(parms.condition_weights,1) ~= nconds
  error('number of condition_weights (%d) does not match number of columns (conditions) of fnamelist (%d)',...
    size(parms.condition_weights,1),nconds);
end;

if ~isempty(parms.frames) && min(parms.frames)<1
  error('frame numbers must be > 1');
end;

if parms.verbose_flag
  tic
  fprintf('%s: checking input files...\n',mfilename);
end;
% check files
results.orig_volsz = [];
results.M = [];
rejectlist = zeros(nsubs,1);
for s=1:nsubs
  for c=1:nconds
    if isempty(fnamelist{s,c})
      rejectlist(s) = 1;
      continue;
    end;
    [vol,M,mrparms,volsz] = fs_load_mgh(fnamelist{s,c},[],[],1); % header-only
    if isempty(parms.frames), parms.frames = [1:volsz(4)]; end;
    if max(parms.frames)>volsz(4)
      error('only %d frames in file %s',volsz(4),fnamelist{s,c});
    end;
    if isempty(results.orig_volsz)
      results.orig_volsz = volsz;
      results.M = M;
    end;
    if any(results.orig_volsz~=volsz)
      error('size of input volume for\n%s (%d,%d,%d,%d)\ndoes not match that of\n%s (%d,%d,%d,%d)',...
        fnamelist{s,c},volsz(1),volsz(2),volsz(3),volsz(4),...
        fnamelist{1,1},results.orig_volsz(1),results.orig_volsz(2),results.orig_volsz(3),results.orig_volsz(4));
    end;
  end;
end;
if parms.verbose_flag, toc; end;

parms.gnamelist = {parms.gnamelist{~rejectlist}};
results.nvals = prod(results.orig_volsz(1:3));
results.nframes = length(parms.frames);
results.volsz = [results.orig_volsz(1:3) results.nframes];
results.frames = parms.frames;

nsubs_keep = length(find(~rejectlist));
tmp_fnamelist = cell(nsubs_keep,nconds);
k=1;
for s=1:nsubs
  if ~rejectlist(s)
    for c=1:nconds
      tmp_fnamelist{k,c} = fnamelist{s,c};
    end;
    k=k+1;
  end;
end;
fnamelist = tmp_fnamelist;
nsubs = nsubs_keep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare contrasts

groupNames = unique(parms.gnamelist);
numGroups = length(groupNames);
gnumlist = zeros(length(parms.gnamelist));
c=1;
results.groups = [];
if parms.contrasts_flag
  results.contrasts = [];
end;
for i=1:length(groupNames)
  gnumlist(find(strcmp(groupNames{i},parms.gnamelist)))=i;
  results.groups(i).name = groupNames{i};
  results.groups(i).n = 0;
  results.groups(i).mean = zeros(results.nvals,results.nframes,'single');
  if parms.stats_flag
    results.groups(i).stdv = zeros(results.nvals,results.nframes,'single');
    results.groups(i).df = [];
    results.groups(i).tstats = [];
    results.groups(i).pvals = [];
  end;
  if parms.contrasts_flag
    for j=1:length(groupNames)
      if i==j, continue; end;
      results.contrasts(c).name = sprintf('%s-VS-%s',groupNames{i},groupNames{j});
      results.contrasts(c).n = [];
      results.contrasts(c).mean = [];
      if parms.stats_flag
        results.contrasts(c).df = [];
        results.contrasts(c).stdv = [];
        results.contrasts(c).tstats = [];
        results.contrasts(c).pvals = [];
      end;
      c = c + 1;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data and calculate sums and sums of squares
for f=1:results.nframes
  frame = parms.frames(f);
  if parms.verbose_flag
    tic
    fprintf('%s: summing data for frame %d of %d...\n',mfilename,frame,results.orig_volsz(4));
  end;
  for s=1:nsubs
    data = zeros(results.nvals,nconds,'single');
    for c=1:nconds
      vec = fs_load_mgh(fnamelist{s,c},[],frame);
      vec = reshape(vec,[prod(size(vec)) 1]);
      data(:,c) = vec - parms.offset;
    end;
    vec = data*parms.condition_weights;
    gnum = gnumlist(s);
    if f==1 % add to the group only once
      results.groups(gnum).n = results.groups(gnum).n + 1;
    end;
    results.groups(gnum).mean(:,f) = results.groups(gnum).mean(:,f) + vec;
    if parms.stats_flag
      results.groups(gnum).stdv(:,f) = results.groups(gnum).stdv(:,f) + vec.^2;
    end;
  end;
  if parms.verbose_flag, toc; end;
end;

% calculate mean, stdv, tstats, pvals for single conditions
if parms.verbose_flag && parms.stats_flag
  tic
  fprintf('%s: calculating stats...\n',mfilename);
end;
for i=1:length(groupNames)
  n = results.groups(i).n;
  if n==0 % should not happen
    error('group %s has 0 subjects',groupNames{i});
  elseif n==1
    if parms.verbose_flag
      fprintf('%s: WARNING: group %s has 1 subject\n',mfilename,groupNames{i});
    end;
    if parms.stats_flag
      results.stdv = zeros(results.nvals,results.nframes,'single');
    end;
    continue;
  end;
  results.groups(i).mean = reshape(results.groups(i).mean,results.volsz);
  if parms.stats_flag
    results.groups(i).stdv = reshape(results.groups(i).stdv,results.volsz);
    results.groups(i).stdv = sqrt(max(0,(n*results.groups(i).stdv - results.groups(i).mean.^2))./(eps+n*(n-1)));
    results.groups(i).mean = results.groups(i).mean / (eps+n);
    results.groups(i).tstats = results.groups(i).mean ./ (results.groups(i).stdv ./ sqrt(eps+n));
    results.groups(i).tstats(results.groups(i).stdv==0)=0;
    results.groups(i).df = n-1;

    results.groups(i).pvals = (results.groups(i).tstats<0).*...
                                (log10(2*tcdf(results.groups(i).tstats,...
                                 results.groups(i).df))) +...
                              (results.groups(i).tstats>0).*...
                                (-log10(2*tcdf(-results.groups(i).tstats,...
                                 results.groups(i).df))); % 2-tailed
  else
    results.groups(i).mean = results.groups(i).mean / (eps+n);
  end;
end;
if parms.verbose_flag && parms.stats_flag, toc; end;

% calculate mean diffs, pooled stdv, two-sample tstats, pvals for contrasts
if parms.contrasts_flag
  c=1;
  for i=1:length(groupNames)
    for j=1:length(groupNames)
      if i==j, continue; end;
      n1 = results.groups(i).n;
      n2 = results.groups(j).n;
      n = n1+n2;
      results.contrasts(c).n = n;
      results.contrasts(c).mean = results.groups(i).mean - results.groups(j).mean;
      if parms.stats_flag
        results.contrasts(c).df = results.groups(i).df + results.groups(j).df;
        % assume equal variance for t-test
        results.contrasts(c).stdv = sqrt(((n1-1)*results.groups(i).stdv.^2+...
                                          (n2-1)*results.groups(j).stdv.^2)/(n-2));
        %% todo: t-test with unequal variance
        results.contrasts(c).tstats = results.contrasts(c).mean ./...
                                      (results.contrasts(c).stdv * sqrt(eps + n/(n1*n2)));
        results.contrasts(c).tstats(results.contrasts(c).stdv==0)=0;
        results.contrasts(c).pvals = (results.contrasts(c).tstats<0).*...
                                       (log10(2*tcdf(results.contrasts(c).tstats,...
                                        results.contrasts(c).df))) +...
                                     (results.contrasts(c).tstats>0).*...
                                      (-log10(2*tcdf(-results.contrasts(c).tstats,...
                                       results.contrasts(c).df))); % 2-tailed
      end;
      c = c + 1;
    end;
  end;
end;

