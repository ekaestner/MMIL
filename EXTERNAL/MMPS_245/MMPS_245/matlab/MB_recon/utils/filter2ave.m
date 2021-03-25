function im = filter2ave(im, n)
%
% function im = filter2ave(im, [n=3])
%
% Two-dimensional averaging digital filter.
%
% Inputs
%   im - Input 2-D images. Dimension: [X, Y, 3rd Dimension, 4th Dimension, ...].
%   n  - Number of points for the averaging filter. Dimension: a
%        scalar (same for X and Y) or a vector with 2 elements.
% Output
%   im - Average filtered 2-D images. Dimension same as the input 'im'.
%
% (c) Kangrong Zhu      Stanford University     Oct 2014

if ~exist('n', 'var') || isempty(n)
    n = 3;
end

if isscalar(n)
    n = [n, n];
end

if (n(1) == 1) && (n(2) == 1)
    return;
end

h = ones(n(1), n(2)) ./ ( n(1) * n(2) );

sz = size(im);
nim = prod( sz(3 : end) );
for im_idx = 1 : nim
    im(:, :, im_idx) = filter2(h, im(:, :, im_idx));
end

return
