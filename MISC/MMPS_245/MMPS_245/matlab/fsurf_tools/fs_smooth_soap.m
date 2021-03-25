function [smoothedvals] = fs_smooth_soap(surf,vals,niter,nsmooth)
%function [smoothedvals] = fs_smooth_soap(surf,vals,niter,[nsmooth])
%
% Required Input:
%   surf is a structure containing:
%     vertices: x,y,z coordinates for each surface vertex
%     nbrs:   vertex numbers of neighbors for each vertex
%     nverts: number of vertices
%
%   vals: vector with nverts members
%   niter: number of iterations (soap bubble smoothing steps)
%
% Optional Input:
%   nsmooth: number of regular smoothing steps per soap bubble iteration
%
% Output:
%   smoothedvals: vector with nverts members
%
% created:  06/11/06 Don Hagler
% Last Mod: 02/24/10  by Don Hagler
%

smoothedvals = [];

if nargin < 3
  help(mfilename);
  return;
end;

if isempty(vals)
  error('vals is empty');
end;

if size(vals,2)~=1
  error('vals must have a single column (it has %d)',size(vals,2));
end;

if size(vals,1)~=surf.nverts
  error('number of vals (%d) does not match number of verts (%d)',...
    size(vals,1),surf.nverts);
end;

if ~exist('nsmooth','var')
  nsmooth = 1;
end;
if isempty(nverts_only)
  nsmooth = 1;
end;

%fprintf('%s(%d,%d): soap bubble smoothing',mfilename,niter,nsmooth);
vals = [0;vals]; % create dummy vertex with index 1, value 0
orig_vals = vals;
orig_v = ind(orig_vals);
maxnbrs = size(surf.nbrs,2);
surf.nbrs = [ones(maxnbrs,1)';surf.nbrs + 1]; % nbrs contains zeros, change to 1
num_nbrs = sum(surf.nbrs>1,2)+1; % including vertex and its neighbors
for i=1:niter
  vals(orig_v)=orig_vals;
%  fprintf('*');
  for j=1:nsmooth
%    fprintf('.');
    Y=[vals, vals(surf.nbrs(:,:))]; % sum values from vertex and its neighbors
    vals=sum(Y,2)./num_nbrs;
  end;
end;
smoothedvals = vals(2:end);

%fprintf('\n');

