function [xvals,yvals] = mmil_estimate_bsbc(volA,volB,varargin)
%function [xvals,yvals] = mmil_estimate_bsbc(volA,volB,[options])
%
% Purpose: estimate the intensity transformation required to match
%     images between scanners / scan parameters
%   based on matching intensity histograms
%
% Required Input:
%   volA: reference volume (3D matrix)
%   volB: volume to be matched to reference (3D matrix)
%
% Optional Input:
%   'volmaskA': brain mask for volA
%     if not supplied, will correct values for entire image
%     {default = []}
%   'volmaskB': brain mask for volB
%     if not supplied, will correct values for entire image
%     {default = []}
%   'xvals': vector of pre-transformation intensity values
%     {default = linspace(0,150,1000)}
%   'voxmatch_flag': [0|1] assume voxel-by-voxel match between images
%     {default = 0}
%   'smooth_flag': smooth histograms before calculating transformation
%     {default = 1}
%   'lambda_sm': regularization term to smooth histograms
%     {default = 1e3}
%   'clip_flag': limit transformation to values between p and 1-p
%     {default = 1}
%   'prctile': percentile margin used for clipping intensity range
%     {default = 1e-3}
%
% Output:
%   xvals: vector of original intensity values
%   yvals: vector of transformed intensity values
%
% Created:  01/08/15 by Don Hagler
% Last Mod: 11/11/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xvals = []; yvals = [];
if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'volmaskA',[],[],...
  'volmaskB',[],[],...
  'xvals',linspace(0,150,1000),[],...
  'voxmatch_flag',false,[false true],...
  'smooth_flag',true,[false true],...
  'lambda_sm',1e3,[],...
  'clip_flag',true,[false true],...
  'prctile',1e-3,[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xvals = parms.xvals;

if ~isempty(parms.volmaskA)
  ind_maskA = find(parms.volmaskA>0);
  vecA = volA(ind_maskA);
else
  vecA = volA(:);
end;

if ~isempty(parms.volmaskB)
  ind_maskB = find(parms.volmaskB>0);
  vecB = volB(ind_maskB);
else
  vecB = volB(:);
end;

if parms.voxmatch_flag
  ind_validA = isfinite(vecA+vecB);
  ind_validB = ind_validA;
else
  ind_validA = isfinite(vecA);
  ind_validB = isfinite(vecB);
end;

hc = max(0.01,hist(vecA(ind_validA),xvals));
chcA = cumsum(hc)/sum(hc);
if parms.smooth_flag
  tmp = mmil_logit(chcA,1);
  [dummy tmp] = mmil_smoothfun1d(xvals,tmp,xvals,parms.lambda_sm);
  chcA = mmil_logit(tmp,0);
end;

hc = max(0.01,hist(vecB(ind_validB),xvals));
chcB = cumsum(hc)/sum(hc);
if parms.smooth_flag
  tmp = mmil_logit(chcB,1);
  [dummy tmp] = mmil_smoothfun1d(xvals,tmp,xvals,parms.lambda_sm);
  chcB = mmil_logit(tmp,0);
end;

if ~parms.clip_flag
  yvals = interp1(chcA,xvals,chcB,'linear','extrap');
else
  % indices defining range of values within prctile range
  iA1 = max(find(chcA<=parms.prctile))+1;
  iA2 = min(find(chcA>=1-parms.prctile))-1;
  iB1 = max(find(chcB<=parms.prctile))+1;
  iB2 = min(find(chcB>=1-parms.prctile))-1;
  % values and cumulative histograms within prctile range
  chcA_clipA = chcA(iA1:iA2);
  xvals_clipA = xvals(iA1:iA2);
  chcB_clipB = chcB(iB1:iB2);
  xvals_clipB = xvals(iB1:iB2);
  % transformation values within prctile range
  yvals_clipB = interp1(chcA_clipA,xvals_clipA,chcB_clipB,...
                        'linear','extrap');
  % transformation values for entire range
  yvals = zeros(size(xvals));
  yvals(1:iB1) = ...
    linspace(min(parms.xvals),yvals_clipB(1),iB1);
  yvals(iB2:end) = ...
    linspace(yvals_clipB(end),max(parms.xvals),length(yvals)-iB2+1);
  yvals(iB1:iB2) = yvals_clipB;
end;

