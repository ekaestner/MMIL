function R = ts_nbrhd_covar(surf,a0,a1,a2,tang_flag)
%function R = ts_nbrhd_covar(surf,[a0],[a1],[a2],[tang_flag])
%
% create neighborhood smoothness constraint covariance matrix
%
% Required Parameters:
%   surf: surface struct (see fs_read_surf.m)
%
% Optional Parameters:
%   a0: variance weight for each vertex in ROI
%     {default = 1.0}
%   a1: covariance weight for first degree neighbor
%     {default = 0.5}
%   a2: covariance weight for second degree neighbor
%     {default = 0.25}
%   tang_flag: [0|1] whether to include elements for tangential components
%     {default: 1}
%
% Output:
%   R: sparse source covariance matrix square in number of sources
%        if tang_flag=0, number of sources = num vertices
%        if tang_flag=1, number of sources = num vertices * 3
%
% Created:  06/16/15 by Don Hagler
% Last Mod: 05/31/16 by Don Hagler
%

R = [];

if ~exist('a0','var') | isempty(a0), a0 = 1.0; end;
if ~exist('a1','var') | isempty(a1), a1 = 0.5; end;
if ~exist('a2','var') | isempty(a2), a2 = 0.25; end;
if ~exist('tang_flag','var') | isempty(tang_flag), tang_flag = 1; end;
smf = 10^-5;

% smooth for each vertex
fprintf('%s: setting covariance values of neighbors for each vertex...\n',mfilename);
tic
surfQ = preprocessQ(surf);
tmpR = sparse(surf.nverts,surf.nverts);
for k=1:surf.nverts
  % first degree neighbors
  v1 = neighborsQ(surfQ,k,1);
  % second degree neighbors
  v2 = neighborsQ(surfQ,k,2);
  v2 = setdiff(v2,v1);
  % set covariance values
  tmpR(k,k) = a0;
  tmpR(k,v1) = a1;
  tmpR(k,v2) = a2;
end;
toc

% interleave tangential components
if tang_flag
  fprintf('%s: interleaving tangential components...\n',mfilename);
  tic
  ndips = 3*surf.nverts;
  R = speye(ndips,ndips);
  R(1:3:end,1:3:end) = tmpR;
  toc
else
  R = tmpR;
end;

