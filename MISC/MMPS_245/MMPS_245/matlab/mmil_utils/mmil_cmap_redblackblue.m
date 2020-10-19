function map = mmil_cmap_redblackblue()
%function map = mmil_cmap_redblackblue()
%
% Created:  01/25/07 by Don Hagler
% Last Mod: 09/27/12 by Don Hagler
%

map = mmil_cmap_blueblackred();
map=map(end:-1:1,:);

