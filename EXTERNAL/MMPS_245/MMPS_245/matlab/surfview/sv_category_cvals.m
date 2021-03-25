function cvals = sv_category_cvals(vals,varargin)
%function cvals = sv_category_cvals(vals,[options])
%
% Purpose: create matrix of rgb color values (0 - 1) from vals
%
% Required Input:
%   vals: vector of values with size [nverts,ncategory]
%     each column of vals specifies weights for a category
%     e.g. fuzzy clusters
%
% Optional Input:
%   'cmap': color map matrix with size [n,3]
%     alternatively, may be any string that is a valid input
%       for colormap function (e.g. 'hot', 'gray', 'jet', etc.)
%     {default = 'jet'}
%   'curvvals': vector of cortical surface curvature values
%     to modulate gray level
%     {default = []}
%   'curvfact': gray-scale contrast of curvature values
%     {default = 0.2}
%
% Output:
%   cvals: matrix of rgb color values (0 - 1) with size [nverts,3]
%
% Created:  09/21/12 by Don Hagler (based on code by Anders Dale)
% Last Mod: 06/29/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
cvals = [];

parms = mmil_args2parms(varargin,{...
  'cmap','jet',[],...
  'curvvals',[],[],...
  'curvfact',0.2,[0,1],...
...
  'fmax',1,[],...
  'fmin',0,[],...
  'fmid',0.5,[],...
});

nvals = size(vals,1);
ncategory = size(vals,2);

if isempty(parms.curvvals)
  curvvals = -ones(nvals,1);
else
  curvvals = mmil_colvec(parms.curvvals);
  if length(curvvals) ~= nvals
    error('length of curvvals (%d) does not match vals (%d)',...
      length(curvvals),nvals);
  end;
end;

if ischar(parms.cmap)
  figure;
  set(gcf,'visible','off');
  cmap = colormap(parms.cmap);
  close(gcf);  
else
  cmap = parms.cmap;
end;
if size(cmap,2) ~= 3, error('cmap must have size = [n,3]'); end;
ncols = size(cmap,1);

cmax = ncategory;
cmin = 1;
cmid = (cmin+cmax)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate gray-scale colors from curvature
cvals_curv = parms.curvfact*(1-parms.curvfact*sign(curvvals))*ones(1,3);

% mapping from category number to scaling factor
cvec = ([1:cmax]-cmin)/(cmax-cmin);

% index to colormap from categories
civec = 1+cvec*(size(cmap,1)-1);

% set colors for each category
category_colors = interp1(cmap,civec,'linear*','extrap');

% ignore negative values and scale to max
avals = vals;
avals(vals<0) = 0;
avals = bsxfun(@rdivide,avals,max(avals,1)); 

% set colors for each vertex
cvals_vals = avals*category_colors;

% calculate maximum value for each vertex
mvals = max(avals,[],2);

% initialize vector of values based on mvals
fvals = zeros(size(mvals));

% set fvals for avals between fmin and fmid
ivec = find(mvals>parms.fmin & mvals<=parms.fmid);
fvals(ivec) = (0.5*(mvals(ivec)-parms.fmin)/(parms.fmid-parms.fmin));

% set fvals for avals greater than fmid
ivec = find(mvals>parms.fmid);
fvals(ivec) = (0.5+0.5*(mvals(ivec)-parms.fmid)/(parms.fmax-parms.fmid));

% clip fvals > 1
fvals = min(1,fvals);

% calculate weighting for each vertex
wvec = 2*min(fvals,0.5);

% combine color vals for curv and vals
wmat = repmat(wvec,[1 3]);
cvals = (1-wmat).*cvals_curv + wmat.*cvals_vals;

% clip values greater than one (overlap)
cvals(cvals>1) = 1;

