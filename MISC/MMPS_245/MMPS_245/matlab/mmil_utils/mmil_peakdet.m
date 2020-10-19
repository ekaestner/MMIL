function [maxtab, mintab]=mmil_peakdet(v,delta,x)
%function [maxtab, mintab]=mmil_peakdet(v,delta,x)
%
% Purpose: detect peaks in a vector
%
% Required Input:
%   v: vector of values
%   delta: minimum difference between adjacent minima and maxima
%     a point is considered a maximum peak if it has the maximal
%     value, and was preceded (to the left) by a value lower by delta.
%
% Optional Input:
%   x: vector of X-values corresponding to v
%     if not supplied, maxtab and mintab container indices
%
% Output:
%   maxtab: matrix with two columns containing maxima
%   mintab: matrix with two columns containing minima
%    column 1 contains indices in V, and column 2 the found values.
%    if x is supplied, column 1 contains corresponding X-values
%
% Created:  by Eli Billauer
% Last Mod: 04/07/14 by Don Hagler
%

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.

maxtab = [];
mintab = [];

v = v(:); % Just in case this wasn't a proper vector

if nargin < 3
  x = (1:length(v))';
else 
  x = x(:);
  if length(v)~= length(x)
    error('Input vectors v and x must have same length');
  end
end
  
if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta <= 0
  error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1:length(v)
  this = v(i);
  if this > mx, mx = this; mxpos = x(i); end
  if this < mn, mn = this; mnpos = x(i); end
  
  if lookformax
    if this < mx-delta
      maxtab = [maxtab ; mxpos mx];
      mn = this; mnpos = x(i);
      lookformax = 0;
    end  
  else
    if this > mn+delta
      mintab = [mintab ; mnpos mn];
      mx = this; mxpos = x(i);
      lookformax = 1;
    end
  end
end
