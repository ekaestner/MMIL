function [snr_map, sig_map, noi_map] = imsnr(im, rep_dim, fit_order, debug)
%
% function [snr_map, sig_map, noi_map] = imsnr(im [, rep_dim, fit_order, debug])
%
% Calculates the SNR of the input images which have repeated measurements.
%
% Input:
%   im        - Input image matrix.
%   rep_dim   - The dimension for the repeated measurements in 'im'. This
%               should be the temporal dimension if a temporal series is used.
%               Default: The last dimension in 'im'.
%   fit_order - The order for the signal drift fitting in the repeatedly
%               measured dimension. Default: 4.
%   debug     - true: Print out the current image index that is being calculated;
%               false: No printing. Default: false.
%
% Output:
%   snr_map   - The SNR maps. Dim: the same as 'im', except that the
%               repeatedly measured dimension now has a size of 1.
%   sig_map   - The signal maps. Dim: the same as 'snr_map'.
%   noi_map   - The noise maps. Dim: the same as 'snr_map'.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

% Set defaults if not passed.
dat_siz = size(im);

if ~exist('rep_dim', 'var') || isempty(rep_dim)
    rep_dim = length(dat_siz);
end

if ~exist('fit_order', 'var') || isempty(fit_order)
    fit_order = 4;
end

if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end

% Constants and matrix sizes.
USE_DEFAULT_NORMALIZATION = 0; % For function 'std', use the default normalization by N-1.
REP_DIM_RESHAPED          = 4; % The repeatedly measured dimension after reshaping the input image matrix.

dat_dim     = 1 : length(dat_siz);
non_rep_dim = dat_dim(dat_dim ~= rep_dim);          % Dimensions not corresponding to the repeated measurements.
nx          = dat_siz(1);
ny          = dat_siz(2);
nim         = prod( dat_siz(non_rep_dim(3 : end)) ); % # of images to calculate the SNR for.
nt          = dat_siz(rep_dim);                      % # of repeated measurements for each image.
t           = (1 : nt).';                            % Indices of the repeated measurements. Used in removing the drifting in a time series.

% Reshape the image matrix, so that the last dimension corresponds to the repeated measurements.
im = reshape( permute(im, [non_rep_dim, rep_dim]), [nx, ny, nim, nt]);

% Fit and remove a polynomial to each pixel time course.
src_mtx = zeros(nt, fit_order+1);                  % Matrix for indices of the temporal phases.
for order = 0 : fit_order
    src_mtx(:, order+1) = t.^order;
end

sig_map = zeros(nx, ny, nim);
for im_ind = 1 : nim
    
    if debug
        fprintf('Calculating SNR for image %d/%d...\n', im_ind, nim);
    end
    
    for x = 1 : nx
        for y = 1 : ny
            tgt_mtx = reshape( im(x, y, im_ind, :), [nt, 1]);
            poly_coe = pinv(src_mtx) * tgt_mtx;    % The polynomial coefficients.
            
            polynomial = src_mtx * poly_coe;
            im(x, y, im_ind, :) = im(x, y, im_ind, :) - ...
                reshape(polynomial, [1, 1, 1, nt]); % Only signal fluctuation is left in 'im'.
            
            sig_map(x, y, im_ind) = mean(polynomial);
            
        end
    end
end

% Calculate the SNR maps.
noi_map = std(im, USE_DEFAULT_NORMALIZATION, REP_DIM_RESHAPED);
snr_map = sig_map ./ noi_map;
snr_map(sig_map < 10^(-12)) = 0;

% Reshape the SNR maps to the original matrix size.
dat_siz(rep_dim) = 1;
snr_map = reshape(snr_map, dat_siz);

return;
