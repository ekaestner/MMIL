function [smoothedvals] = fs_smooth_sparse(surf,vals,niter)
%function [smoothedvals] = fs_smooth_sparse(surf,vals,niter)
%
% surf is a structure containing:
%   vertices: x,y,z coordinates for each surface vertex
%   nbrs:   vertex numbers of neighbors for each vertex
%   nverts: number of vertices
%
% vals is vector with nverts members
%
% niter is number of iterations (smoothing steps)
%
% Created:  01/20/08  by Don Hagler
% Last Mod: 10/24/10  by Don Hagler
%

smoothedvals = [];

if ~mmil_check_nargs(nargin,3), return; end;

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

%fprintf('%s(%d): smoothing',mfilename,niter);
vals = [0;vals]; % create dummy vertex with index 1, value 0
fixedvals = find(vals);
maxnbrs = size(surf.nbrs,2);
nbrs = [ones(maxnbrs,1)';surf.nbrs + 1]; % nbrs contains zeros, change to 1
nbrs(fixedvals,:)=1; % don't let fixedvals get changed by neighbors

for iter=1:niter
%  fprintf('.');
  Y=[vals, vals(nbrs(:,:))]; % sum values from vertex and neighbors
  nzvec = find(max(abs(Y),[],2)); % only include non-zero values
  vals(nzvec)=sum(Y(nzvec,:),2)./sum(Y(nzvec,:)~=0,2);
end;
smoothedvals = vals(2:end);

%fprintf('\n');

