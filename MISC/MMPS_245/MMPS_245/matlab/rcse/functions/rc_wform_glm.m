function [results,data_fit] = rc_wform_glm(data,time,varargin)
%function [results,data_fit] = rc_wform_glm(data,time,[options])
%
% Required Input:
%   data: time series data; with size [ntpoints,nconds]
%   time: time vector; with size [ntpoints,1]
%
% Optional Parameters for general linear model:
%   'regressors': matrix of regressors; with size [nconds,nregs]
%     if not supplied,  will test model with intercept  only (mean)
%     {default = []}
%   'regnames': cell array of regressor names
%     {default = []}
%   'npoly': order of polynomial fit
%     {default = 1}
%
% Output:
%   results: struct containing  these fields:
%     data: input data matrix; size = [ntpoints,nconds]
%     time: input time vector; size = [ntpoints,1]
%     ntpoints: number of time points
%     nconds: number of conditions
%     glm: structure containing GLM betas, contrasts, etc.
%   data_fit: time series data; with size [ntpoints,nconds]
%
% Created:  04/05/16 by Don Hagler
% Last mod: 04/05/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
data_fit = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'regressors',[],[],...
  'regnames',[],[],...
  'npoly',1,[1:5],...
});

%% todo: input parm regressor names

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure ndims of data is 2
results.data = data;
if ndims(results.data) ~= 2
  error('data must be a 2-dimensional matrix');
end;

% make sure time is a vector
results.time = mmil_colvec(time);

% get dimensions of input data
[results.ntpoints,results.nconds] = size(data);

% check for mismatch between data and time
if results.ntpoints ~= length(results.time)
  error('number of columns in data matrix does not match number of elements in time vector');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% contstruct design matrix from regressors
X = ones(results.nconds,1);
if ~isempty(parms.regressors)
  nregs = size(parms.regressors,2);
  % check regnames
  if isempty(parms.regnames)
    for i=1:nregs
      parms.regnames{i} = sprintf('reg%d',i);
    end;
  end;
  % construct design matrix
  for i=1:nregs
    for j=1:parms.npoly
      X = cat(2,X,parms.regressors(:,i).^j);
    end;
  end;
else
  nregs = 0;
end;
nparms = size(X,2);

% create contrast vectors for glm
contrast_names = cell(nregs,1);
contrast_vectors = cell(nregs,1);
for i=1:nregs
  contrast_names{i} = parms.regnames{i};
  vec = zeros(1,nparms);
  j = 1 + parms.npoly*(i-1) + 1;
  k = j + parms.npoly - 1;
  vec(j:k) = 1;  
  contrast_vectors{i} = vec;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = data';
results.glm = mmil_glm_calc(X,D,...
  'contrast_vectors',contrast_vectors,...
  'contrast_names',contrast_names);

data_fit = (X*results.glm.betas)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

