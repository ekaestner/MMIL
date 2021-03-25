function [vol,beta] = mmil_regress_vol(vol,X,ind_valid)
%function [vol,beta] = mmil_regress_vol(vol,X,ind_valid)
%
% Purpose: use linear regression to project out
%   regressor time-series from each voxel of input time-series volume
%
% Required Input:
%   vol: timeseries volume (4D)
%     size should be [nx,ny,nz,nt]
%   X: matrix of regressor time-series
%     size should be [nt,nr]
%
% Optional Input:
%   ind_valid: vector of time point index numbers
%     to use in calculation of beta coefficients
%
% Output:
%   vol: post-regression timeseries volume (4D)
%     with size [nx,ny,nz,nt]
%   beta: vector of beta coefficients
%     with size [1,nr]
%
% Created:  07/13/12 Don Hagler
% Prev Mod: 11/01/12 Don Hagler
% Last Mod: 11/21/17 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('ind_valid','var'), ind_valid = []; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get size of volume
volsz = size(vol);
nvox = prod(volsz(1:3));
nframes = volsz(4);

% get size of tseries matrix
[nt,nr] = size(X);

if nt ~= nframes
  error('first dim of X must match last dim of vol');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reshape 4D volume into 2D matrix
V = reshape(vol,[nvox,nframes]);

% calculate beta coefficients
if isempty(ind_valid)
  beta = (pinv(X)*V')';
else
  beta = (pinv(X(ind_valid,:))*V(:,ind_valid)')';
end;

% calculate model fit
Y = beta*X';

% subtract fit from data
V = V - Y;

% reshape to original volume
vol = reshape(V,volsz);

