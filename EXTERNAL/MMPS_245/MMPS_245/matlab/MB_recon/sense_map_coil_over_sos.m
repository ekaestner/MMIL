function [smap, im] = sense_map_coil_over_sos(smap_ksp, smap_sz, coil_dim, crop_smap)
%
% function [smap, im] = sense_map_coil_over_sos(smap_ksp, smap_sz [, coil_dim, crop_smap])
%
% Calculate the sensitivity maps using the coil images divided by the 
% square root of sum of square(SOS) image.
%
% Inputs
%    smap_ksp  - K-space data for sensitivity map calculation. Dim: [Kx, Ky, ...].
%    smap_sz   - The size for the 2D sensitivity maps. The input rawdata
%                will be zero-padded to this size in the first 2 dimensions. 
%                Dim: [FE, PE, ...](Only first 2 elements will be used).
%    coil_dim  - Dimension for the coil in the input 'smap_ksp'. Default:
%                last dimension in 'smap_ksp'.
%    crop_smap - True: crop the sensitivity maps based on an intensity threshold.
%
% Outputs
%    smap      - The sensitivity maps, size: [smap_sz(1), smap_sz(2), ...].
%                The 3rd onward dimensions are the same as those in the 
%                input 'smap_ksp'.
%    im        - The square root of sum of square(SOS) image.
%
% (c) Kangrong Zhu,     Stanford University     April 2013

% Input data size
in_dat_sz = size(smap_ksp);

if ~exist('coil_dim', 'var') || isempty(coil_dim)
    coil_dim = length(in_dat_sz);
end

% Output data size
out_dat_sz = in_dat_sz;
out_dat_sz(1 : 2) = smap_sz(1 : 2);

% Pad the raw data to a size of [smap_sz(1), smap_sz(2)] in the first 2 dimensions
padded = zpad(smap_ksp, out_dat_sz);

% Sensitivity maps - coil images divided by sos image
ic = ifft2c(padded);
im = sos(ic, coil_dim);

repmat_sz = ones(1, length(out_dat_sz));
repmat_sz(coil_dim) = out_dat_sz(coil_dim);
im_rep = repmat(im, repmat_sz);
smap = ic ./ im_rep;

% Crop the sensitivity maps
if crop_smap
    imsz = size(im);
    msk = zeros(imsz);
    nims = prod(imsz(3:end));
    se = strel('disk', 2);
    for im_idx = 1 : nims
        im_temp = im(:, :, im_idx);
        im_temp = im_temp ./ max(abs(im_temp(:)));
        % thresh = graythresh(im_temp);
        thresh = 0.1;
        msk(:, :, im_idx) = im2bw(im_temp, thresh);
        msk(:, :, im_idx) = elimiso(msk(:, :, im_idx), 200);
        msk(:, :, im_idx) = imdilate(msk(:, :, im_idx), se);
    end
    msk = repmat(msk, repmat_sz);
    smap = smap .* msk;
end

return