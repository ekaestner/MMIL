function im = get_montage_mtx(im, msz)
%
% Get a matrix corresponding to a montage of the input images.
%
% Input
%   im  - Input images. The first 2 dimensions correspond to the 2 dimensions in an image.
%   msz - Montage size, [number of rows, number of columns] in montage.
%
% Output
%   im  - The matrix corresponding to a montage of the input images.
%
% (c) Kangrong Zhu  April 2015

sz = size(im);
nim = prod(sz(3:end));
im = reshape(im, [sz(1), sz(2), nim]);

if ~exist('msz', 'var') || isempty(msz)
    [nrow, ncol] = get_montage_nrow_ncol(nim);
else
    nrow = msz(1);
    ncol = msz(2);
end

im = cat(3, im, zeros(sz(1), sz(2), nrow*ncol - nim));
im = permute(im, [2,1,3]);
im = reshape(im, [sz(2), sz(1), ncol, nrow]);
im = permute(im, [1,3,2,4]);
im = reshape(im, [sz(2)*ncol, sz(1)*nrow]);
im = permute(im, [2,1]);

return