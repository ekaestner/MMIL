function m = mmil_cell2num(c)
%function m = mmil_cell2num(c)
%
% This function takes a cell matrix of strings and converts
% each cell into a number.
%
% Last Mod:  03/23/11 by Don Hagler
%

sz = size(c);
m = zeros(sz);
for i=1:prod(size(c))
  m(i) = sscanf(c{i},'%f');
end;

