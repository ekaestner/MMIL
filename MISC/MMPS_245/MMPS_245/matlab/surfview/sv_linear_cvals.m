function cvals = sv_linear_cvals(vals,varargin)
%function cvals = sv_linear_cvals(vals,[options])
%
% Purpose: create matrix of rgb color values (0 - 1) from vals
%   using linear scaling between fmin, fmid, and fmax values
%
% Required Input:
%   vals: vector of values with size [nverts,1]
%
% Optional Input:
%   'fmax': value corresponding to end of color scale
%     if empty, will be set to max(vals)
%     {default = []}
%   'fmin': value corresponding to start of color scale
%     if empty, will be set to fmax/10
%     {default = []}
%   'fmid': value corresponding to middle of color scale
%     if empty, will be set to (fmin+fmax)/2
%     {default = []}
%   'cmap': color map matrix with size [n,3]
%     alternatively, may be any string that is a valid input
%       for colormap function (e.g. 'hot', 'gray', 'jet', etc.)
%     {default = 'mmil_cmap_blueblackred'}
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
% Last Mod: 09/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
cvals = [];

parms = mmil_args2parms(varargin,{,...
  'fmax',[],[],...
  'fmin',[],[],...
  'fmid',[],[],...
  'cmap','mmil_cmap_blueblackred',[],...
  'curvvals',[],[],...
  'curvfact',0.2,[0,1],...
});

vals = mmil_colvec(vals);
nvals = length(vals);
avals = abs(vals);

if isempty(parms.curvvals)
  curvvals = -ones(nvals,1);
else
  curvvals = mmil_colvec(parms.curvvals);
  if length(curvvals) ~= nvals
    error('length of curvvals (%d) does not match vals (%d)',...
      length(curvvals),nvals);
  end;
end;

if isempty(parms.fmax), parms.fmax = max(avals); end;
if isempty(parms.fmin), parms.fmin = parms.fmax/10; end;
if isempty(parms.fmid), parms.fmid = (parms.fmin+parms.fmax)/2; end;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate gray-scale colors from curvature
cvals_curv = parms.curvfact*(1-parms.curvfact*sign(curvvals))*ones(1,3);

% initialize vector of values based on vals
fvals = zeros(size(vals));

% set fvals for avals between fmin and fmid
ivec = find(avals>parms.fmin & avals<=parms.fmid);
fvals(ivec) = sign(vals(ivec)).*...
             (0.5*(avals(ivec)-parms.fmin)/(parms.fmid-parms.fmin));

% set fvals for avals greater than fmid
ivec = find(avals>parms.fmid);
fvals(ivec) = sign(vals(ivec)).*...
             (0.5+0.5*(avals(ivec)-parms.fmid)/(parms.fmax-parms.fmid));

% clip abs(fvals) > 1
fvals = max(-1,min(1,fvals));

% find index to color map for each fval
civec = round(1+ncols/2) +0.5*fvals*(ncols-1);
civec = max(1,min(ncols,civec));

% linearly interpolate from color map to color values
cvals_vals = interp1(cmap,civec,'linear*','extrap');

% combine color vals for curv and vals
wvec = 2*min(abs(fvals),0.5);
wmat = repmat(wvec,[1 3]);
cvals = (1-wmat).*cvals_curv + wmat.*cvals_vals;

