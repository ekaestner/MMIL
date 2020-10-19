function h = plot3d(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14);
% PLOT3D  'plot3' the given three column data
% function h = plot3d(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14);
% x is n*three columns wide, parameters p1-p14 are optional standard plot3
%  parameters.  h is the handle returned by plot3
% The first three columns of x are plotted as one set, then the next three,
%  etc.

% Copyright (c) 1995 by John C. Mosher
% Los Alamos National Laboratory
% Group ESA-MT, MS J580
% Los Alamos, NM 87545
% email: mosher@LANL.Gov
%
% Permission is granted to modify and re-distribute this code in any manner
%  as long as this notice is preserved.  All standard disclaimers apply.

% March 3, 1995 author

cols = size(x,2);		% number of columns

s = 'hp = plot3(x(:,1:3:cols),x(:,2:3:cols),x(:,3:3:cols)'; % basic string

for i = 1:(nargin-1),		% if extra parameters given
  s = [s sprintf(',p%.0f',i)];  % add next parameter
end

s = [s ');'];			% close off plot3 command

eval(s);			% do it!
if(nargout),			% user wanted handle
  h = hp;
end

return
