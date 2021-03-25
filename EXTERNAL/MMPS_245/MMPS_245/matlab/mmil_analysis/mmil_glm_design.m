function X = mmil_glm_design(regressors,varargin)
%function X = mmil_glm_design(regressors,[options])
%
% Purpose: create design matrix for glm
%
% Usage:
%  results = mmil_glm_design(regresors,'key1', value1,...);
%
% Required Input:
%  regressors: matrix containing regressors for each subject
%     each column should contain a different numerical regressor
%     each row should contain values for a different subject
%     size = [nsubs,nregs]
%
% Optional Input:
%  'nconds': number of conditions
%    {default = 1}
%  'nsubs': number of subjects
%    use this if regressors is empty
%    {default = []}
%  'baseflag': [0|1] include column for baseline/intercept
%    {default = 1}
%   
% Output:
%   X: design matrix  size = [nsubs,nparams]
%
% NOTE: this is a fixed effects, non-repeated measures design
%
% Created:  10/21/11 by Don Hagler
% Last Mod: 12/17/12 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms(varargin, { ...
  'nconds',1,[],...
  'nsubs',[],[],...
  'baseflag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nsubs,nregs] = size(regressors);
if nsubs==0
  if ~isempty(parms.nsubs)
    nsubs = parms.nsubs;
  else
    error('regressors is empty and nsubs not specified');
  end;
end;

%% todo: reduce cond regressors by 1?
% if only one condition, consider this 0 conditions
if parms.nconds==1
  nconds = 0;
else
  nconds = parms.nconds;
end;

ncols = nconds + nregs;
if parms.baseflag
  ncols = ncols + 1;
end;
nrows = nsubs*parms.nconds;

X = zeros(nrows,ncols);
if parms.baseflag
  b = 1;
  X(:,b) = ones(nrows,1); % baseline / intercept
else
  b = 0;
end;
if nconds > 1
  j = 1;
  for s=1:nsubs
    for i=1:nconds
      X(j,b+i) = 1;
      if nregs > 0
        X(j,b+nconds+[1:nregs]) = regressors(s,:);
      end;
      j = j + 1;
    end;
  end;
else
  for i=1:nregs
    X(:,b+i) = regressors(:,i);
  end;
end

