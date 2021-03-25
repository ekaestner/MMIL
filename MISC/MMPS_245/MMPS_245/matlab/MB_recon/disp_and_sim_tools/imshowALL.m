function imshowALL(im, dispRange, showSize, newFig)
%
% function imshowALL(im [, dispRange, showSize, newFig])
%
% On the same figure, show all 2D grayscale images in the image matrix 'im'.
%
% Inputs:
%   im        - The image matrix.
%   dispRange - The grayscale range for display.
%   showSize  - The 'Size' argument for matlab built-in function 'montage'.
%   newFig    - 1: open a new figure, 0: no new figure. Default:0.
%
% (c) Kangrong Zhu,     Stanford University     2011

if ~exist('dispRange', 'var') || isempty(dispRange)
    dispRange = [0, max(abs(im(:)))];
end
if ~exist('newFig', 'var')
    newFig = 0;
end

sz = size(im);
im = reshape(im, sz(1), sz(2), 1, prod(sz(3:end)));

if newFig
    figure;
end
if exist('showSize','var')
    montage(im, dispRange, 'Size', showSize); 
else
    montage(im, dispRange);
end
colorbar;