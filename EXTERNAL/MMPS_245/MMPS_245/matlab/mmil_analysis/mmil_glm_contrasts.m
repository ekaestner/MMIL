function output = mmil_glm_contrasts(regressors,varargin)
%function output = mmil_glm_contrasts(regressors,[options])
%
% Purpose: construct contrast vectors and contrast names
%
% Usage:
%  output = mmil_glm_contrasts(regressors,'key1', value1,...);
%
% Required Input:
%  regressors: design matrix for glm
%
% Optional Input:
%  'groupregs': vector of indices for group regressors (e.g. [1 2])
%     used for group contrasts
%     may be cell array of multiple sets of groups
%     if empty, no regressors treated as group variables
%     {default = []}
%  'regnames': cell array of regressor names
%    length should match number of columns in regressors
%    if ommitted, will assign names like reg1, reg2, reg3, etc.
%    {default = []}
%  'nconds': number of conditions
%    {default = 1}
%  'condnames': cell array of condition names
%    length should match nconds
%    if ommitted, will assign names like 'cond1', 'cond2', etc.
%    {default = []}
%  'group_contrast_flag': [0|1] construct contrasts for each pair-wise
%    difference between groupregs
%    {default = 0}
%  'combine_conds_flag': [0|1] whether to combine across conditions
%    {default = 0}
%  'baseflag': [0|1] include separate column for baseline/intercept
%    {default = 1}
%
% Output:
%   output: struct containing fields including
%     contrast_vectors, contrast_names, nsubs, nregs, etc.
%                                
% Created:  10/21/11 by Don Hagler
% Last Mod: 12/17/12 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = [];
if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms(varargin, { ...
  'groupregs',[],[],...
  'regnames',[],[],...
  'nconds',1,[],...
  'condnames',[],[],...
  'group_contrast_flag',false,[false true],...
  'combine_conds_flag',false,[false true],...
  'baseflag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters

if isempty(parms.groupregs)
  nsets = 0;
else
  if ~iscell(parms.groupregs)
    parms.groupregs = {parms.groupregs};
  end;
  nsets = length(parms.groupregs);
  for s=1:nsets
    % check that there are 0, 2, or more groups (not 1)
    ngroups = length(parms.groupregs{s});
    if ngroups==1
      error('number of groups in set %d equals 1',s);
    end;
  end;
end;

% if only one condition, consider this 0 conditions
if parms.nconds==1
  nconds = 0;
else
  %% todo: reduce cond regressors by 1?
  nconds = parms.nconds;
  % set condnames if not supplied
  if ~isempty(parms.condnames)
    if ~iscell(parms.condnames), parms.condnames = {parms.condnames}; end;
    if length(parms.condnames) ~= parms.nconds
      error('number of elements in condnames (%d) does not match nconds (%d)',...
        length(parms.condnames),parms.nconds);
    end;
  else
    parms.condnames = cell(parms.nconds,1);
    for i=1:parms.nconds
      parms.condnames{i} = sprintf('cond%d',i);
    end;
  end;
end;

nsubs = size(regressors,1);
nregs = size(regressors,2);

if nregs
  if ~isempty(parms.regnames)
    if ~iscell(parms.regnames), parms.regnames = {parms.regnames}; end;
    if length(parms.regnames) ~= nregs
      error('number of elements in regnames (%d) does not match nregs (%d)',...
        length(parms.regnames),nregs);
    end;
  else
    parms.regnames = cell(nregs,1);
    for i=1:nregs
      parms.regnames{i} = sprintf('reg%d',i);
    end;
  end;
end;

ind_groups = [];
groupnames = [];
if nsets>=1
  for s=1:nsets
    ind_tmp = parms.groupregs{s};
    ind_groups = cat(2,ind_groups,ind_tmp);
  end;
  groupnames = parms.regnames(ind_groups);
end;
ngroups = length(groupnames);

nvals = nregs - ngroups;
ncols = nconds + ngroups + nvals;
if parms.baseflag
  ncols = ncols + 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create contrast vectors

% initialize
contrast_names = [];
contrast_vectors = [];
k = 1; % index in contrast_vectors array
c = 1; % index in individual contrast vector

% contrast vector for baseline/intercept
zero_vec = zeros(1,ncols);
baseline_vec = zero_vec;
if parms.baseflag
  baseline_vec(1) = 1;
  c = c + 1;
else
  for s=1:nsets
    ind_tmp = parms.groupregs{s};
    ntmp = length(ind_tmp);
    baseline_vec(ind_tmp) = 1/ntmp;
  end;
end;
if parms.baseflag || ngroups>0
  contrast_vectors{k} = baseline_vec;
  contrast_names{k} = 'baseline';
  k = k + 1;
end;

% contrast vectors for each condition
if nconds & ~parms.combine_conds_flag
  for i=1:nconds %% todo: reduce cond regressors by 1?
    tmp = baseline_vec;
    tmp(c) = 1;
    contrast_vectors{k} = tmp;
    contrast_names{k} = parms.condnames{i};
    k = k + 1;
    c = c + 1;
  end;
else
  c = c + nconds;
end;

% contrast vectors for groups and other regressors
for i=1:nregs
  if ismember(i,ind_groups) && parms.baseflag
    tmp = baseline_vec;
  else
    tmp = zero_vec;
  end;
  tmp(c) = 1;
  contrast_vectors{k} = tmp;
  contrast_names{k} = parms.regnames{i};
  k = k + 1;
  c = c + 1;
end;

% between condition contrasts
if ~parms.combine_conds_flag
  for i=1:nconds
    for j=i:nconds
      if i==j, continue; end;
      tmp = zero_vec;
      tmp(1+i) = 1;
      tmp(1+j) = -1;
      contrast_vectors{k} = tmp;
      contrast_names{k} = [parms.condnames{i} '_vs_' parms.condnames{j}];
      k = k + 1;
    end
  end
end;

% between group contrasts
if parms.group_contrast_flag
  if parms.baseflag
    b = 1;
  else
    b = 0;
  end;
  for s=1:nsets
    % only compare groups within each set
    ind_tmp = parms.groupregs{s};
    ntmp = length(ind_tmp);
    tmp_groupnames = parms.regnames(ind_tmp);
    for i=1:ntmp
      for j=i:ntmp
        if i==j, continue; end;
        tmp = zero_vec;
        tmp(b+nconds+ind_tmp(i)) = 1;
        tmp(b+nconds+ind_tmp(j)) = -1;
        contrast_vectors{k} = tmp;
        contrast_names{k} = [tmp_groupnames{i} '_vs_' tmp_groupnames{j}];
        k = k + 1;
      end
    end
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% package output

output.regressors = regressors;
output.regnames = parms.regnames;
output.groupnames = groupnames;
output.condnames = parms.condnames;
output.contrast_vectors = contrast_vectors;
output.contrast_names = contrast_names;
output.ngroups = ngroups;
output.nsubs = nsubs;
output.nregs = nregs;
output.ncols = ncols;
output.nconds = parms.nconds;

