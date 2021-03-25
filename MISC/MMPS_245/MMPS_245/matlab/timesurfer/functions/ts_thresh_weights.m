function [tw,tv]=ts_thresh_weights(w,v,thresh,absflag)
%function [tw,tv]=ts_thresh_weights(w,v,thresh,[absflag])
%
% Required Input:
%   w: vector of weights
%
% Optional Input:
%   v: vector of vertex numbers
%     If empty, assume all vertices are in w
%     {default = []}
%   thresh: threshold applied to values in w
%     {default = 0.01}
%   absflag: [0|1] whether to apply threshold to absolute value
%     {default = 1}
%
% Last Mod: 03/01/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('v','var'), v = []; end;
if ~exist('thresh','var') || isempty(thresh), thresh = 0.01; end;
if ~exist('absflag','var') || isempty(absflag), absflag = 1; end;

tw = w;
tv = v;

if isempty(w), return; end
if isempty(v), v = 1:length(w); end;

if absflag
  iw = find(abs(w)>thresh);
else
  iw = find(w>thresh);
end

if nargout>1
  tw = w(iw);
  tv = v(iw);
else
  tw(~iw) = 0;
end;  

