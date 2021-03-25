function R = ts_roi_nbrhd_covar(surf,roi_verts,a0,a1,a2,ax,az,tang_flag)
%function R = ts_roi_nbrhd_covar(surf,roi_verts,[a0],[a1],[a2],[ax],[az],[tang_flag])
%
% create ROI-neighborhood smoothness constraint covariance matrix
%
% Required Parameters:
%   surf: surface struct (see fs_read_surf.m)
%   roi_verts: cell array of vectors containing vertex numbers for each ROI
%
% Optional Parameters:
%   a0: variance weight for each vertex in ROI
%     {default = 1.0}
%   a1: covariance weight for first degree neighbor
%     {default = 0.5}
%   a2: covariance weight for second degree neighbor
%     {default = 0.25}
%   ax: covariance weight for all other vertices in ROIs
%     {default = 0}
%   az: variance weight for all other vertices in surface
%     {default = 0.1}
%   tang_flag: [0|1] whether to include elements for tangential components
%     {default = 1}
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

if ~mmil_check_nargs(nargin,2), return; end;

if ~exist('a0','var') | isempty(a0), a0 = 1.0; end;
if ~exist('a1','var') | isempty(a1), a1 = 0.5; end;
if ~exist('a2','var') | isempty(a2), a2 = 0.25; end;
if ~exist('ax','var') | isempty(ax), ax = 0.0; end;
if ~exist('az','var') | isempty(az), az = 0.1; end;
if ~exist('tang_flag','var') | isempty(tang_flag), tang_flag = 1; end;
smf = 10^-5;

if ~iscell(roi_verts)
  roi_verts = {roi_verts};
end;
nroi = length(roi_verts);

% smooth for each vertex
fprintf('%s: setting covariance values of neighbors for each vertex...\n',mfilename);
tic
surfQ = preprocessQ(surf);
R = az*speye(surf.nverts,surf.nverts);
for r=1:nroi
  verts = roi_verts{r};
  nverts = length(verts);
  tmpR = eye(nverts,nverts);
  for i=1:nverts
    v0 = verts(i);
    % set variance values
    tmpR(i,i) = a0;
    % set covariance values for first degree neighbors
    v1 = neighborsQ(surfQ,v0,1);
    [v1,j] = intersect(verts,v1);
    if ~isempty(v1)
      tmpR(i,j) = a1;
    end;
    % set covariance values for second degree neighbors
    v2 = neighborsQ(surfQ,v0,2);
    v2 = setdiff(v2,v1);
    [v2,j] = intersect(verts,v2);
    if ~isempty(v2)
      tmpR(i,j) = a2;
    end;
    % set covariance values for all other neighbors
    if ax>0
      vx = setdiff(verts,[v0;v1;v2]);
      [vx,j] = intersect(verts,vx);
      tmpR(i,j) = ax;
    end;
  end;
  % check for very small or negative eigenvalues
  [U,S,V] = svd(tmpR);
  if any(diag(S)<smf)
    S(S<smf) = smf;
    tmpR = U*S*V';
  end;
  R(verts,verts) = tmpR;
end;
toc

% interleave tangential components
if tang_flag
  fprintf('%s: interleaving tangential components...\n',mfilename);
  tic
  ndips = 3*surf.nverts;
  tmpR = R;
  R = speye(ndips,ndips);
  R(1:3:end,1:3:end) = tmpR;
  toc
end;

