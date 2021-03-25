function im = elimiso(im, area_thresh)
%
% function im = elimiso(im, area_thresh)
%
% Eliminate isolated pixels in the input images.
%
% Inputs
%   im          - Input images, background is zero.
%   area_thresh - Regions with area < area_thresh will be eliminated.
%
% Output
%   im          - Images with isolated pixels eliminated, having the same
%                 matrix size as the input 'im'.
%
% (c) Kangrong Zhu      Stanford University     2011

mask = sign(abs(im));
sz = size(im);
nim = prod(sz(3 : end));

for im_idx = 1 : nim
    L = bwlabel(mask(:, :, im_idx));
    stat = regionprops(L, 'Area');
    idx = find([stat.Area] > area_thresh);
    mask(:, :, im_idx) = ismember(L, idx);
end
im = mask .* im;

return