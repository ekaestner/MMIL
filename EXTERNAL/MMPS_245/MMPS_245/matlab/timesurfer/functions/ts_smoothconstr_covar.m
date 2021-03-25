function R=ts_smoothconstr_covar(surf,smooth,tang_flag)
%function R=ts_smoothconstr_covar(surffile,[smooth],[tang_flag])
%
% create smoothness constraint covariance matrix
%
% Required Parameters:
%   surf: surface struct (see fs_read_surf.m)
%
% Optional Parameters:
%   smooth: number of smoothing steps
%     {default: 10}
%   tang_flag: [0|1] whether to include elements for tangential components
%     {default: 1}
%
% Output:
%   R: sparse source covariance matrix square in number of sources
%        if tang_flag=0, number of sources = num vertices
%        if tang_flag=1, number of sources = num vertices * 3
%
% Created:  01/20/08 by Don Hagler
% Last Mod: 06/11/09 by Don Hagler
%

R = [];

if ~exist('smooth','var') | isempty(smooth), smooth = 10; end;
if ~exist('tang_flag','var') | isempty(tang_flag), tang_flag = 1; end;
smf = 10^-5;

% smooth for each vertex
fprintf('%s: smoothing surface for each vertex...\n',mfilename);
tic
tmpR = sparse(surf.nverts,surf.nverts);
for k=1:surf.nverts
  vals = zeros(surf.nverts,1);
  vals(k) = 100;
  smvals = fs_smooth(surf,vals,smooth);
  smvals = smvals/max(smvals);
  smvals(smvals<smf) = 0;
  tmpR(k,:) = smvals;
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

