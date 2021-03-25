function cmap = bluehot(N);
% BLUEHOT make blue hot to red hot dual intensity colorscale
% function cmap = bluehot(N);
% returns N=64, if N not given.  
% Negative values are blue, positive values are red, with intuitive
%  convention that positive values are "warmer" than negative.  

if(exist('N')~=1),
  N = 64;			% default size
end

cmap1 = hot(ceil(N/2));		% in case its odd
cmap2 = hot(floor(N/2));

cmap = [flipud(cmap1(:,[3 2 1])); cmap2];
% cmap = flipud([flipud(cmap);cmap(:,[3 2 1])]); % near is "hot", neg is "cold"

return
