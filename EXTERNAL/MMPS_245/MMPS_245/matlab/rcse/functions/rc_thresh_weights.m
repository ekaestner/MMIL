function [tw,tv]=rc_thresh_weights(w,v,thresh,threshabsflag)
%function [tw,tv]=rc_thresh_weights(w,v,thresh,threshabsflag)
%
% Required Input:
%  w: vector of weights
%  v: vector of corresponding vertex numbers
%  thresh: threshold to be applied to weights
%
% Optional Input:
%   threshabsflag: [0|1] whether to take absolute value of w
%     before applying threshold
%
% Output:
%   tw: vector of weights, excluding those below thresh
%   tv: vector of vertex numbers, excluding those below thresh
%
% Early Mod: 12/22/06 by Don Hagler
% Last Mod:  02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('threshabsflag','var') || isempty(threshabsflag), threshabsflag = 1; end;

tw = w;
tv = v;

if isempty(w), return; end

if threshabsflag
  iw = find(abs(w)>thresh);
else
  iw = find(w>thresh);
end

tw = w(iw);
tv = v(iw);

