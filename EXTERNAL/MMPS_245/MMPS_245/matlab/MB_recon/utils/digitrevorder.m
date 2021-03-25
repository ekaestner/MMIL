% Copyright (C) 2013 Mike Miller
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.

function [y, i] = digitrevorder (x, r)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (~ isvector (x))
    error ('digitrevorder: X must be a vector');
  elseif (~(isscalar (r) && r == fix (r) && r >= 2 && r <= 36))
    error ('digitrevorder: R must be an integer between 2 and 36');
  else
    tmp = log (numel (x)) / log (r);
    if (fix (tmp) ~= tmp)
      error ('digitrevorder: X must have length equal to an integer power of %d', r);
    end
  end

  old_ind = 0:numel (x) - 1;
  new_ind = base2dec (fliplr (dec2base (old_ind, r)), r);

  i = new_ind + 1;
  y(old_ind + 1) = x(i);

  if (size(x,2) == 1)
    y = y';
  else
    i = i';
  end

