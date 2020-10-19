function mmil_subplot_tight(n,m,i)
%function mmil_subplot_tight(n,m,i)
%
% Purpose: create a borderless subplot
%
% Required Parameters:
%   n: number of rows
%   n: number of columns
%   i: index of current subplot
%    
% from: http://www.briandalessandro.com/blog/how-to-make-a-borderless-subplot-of-images-in-matlab/
%
% Created:  11/25/13 by Don Hagler
% Last Mod: 11/25/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

[c,r] = ind2sub([m n], i);
subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 1/n])

return;


