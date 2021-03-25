function [snr_map, noise_map] = pmr_calc_snr(orig, stack, rep_dim)
%
% function snr_map = pmr_calc_snr(orig, stack, [rep_dim = ndims(stack)])
%
% Pseudo Multiple Replica: Calculate SNR maps from image stack.
% Reference: Philip M. Robson, et al. MRM 60:895-907 (2008).
%
% Inputs
%   orig    - The original reconstructed image, that is, that without
%             additional noise added. Its size in dimension 'rep_dim' is 1.
%   stack   - Stack of image replicas. Dimensions are the same as the input
%             'orig', except in the 'rep_dim' dimension.
%   rep_dim - The dimension of replica in the input 'stack'. Size of 'orig'
%             and 'stack' is different only in this dimension.
%
% Output
%   snr_map - SNR maps. Dimension the same as the input 'orig'.
%   noise_map - Noise maps. Dimension the same as the input 'orig'.
%
% (c) Kangrong Zhu      Stanford University     Nov 2014

if ~exist('rep_dim', 'var') || isempty(rep_dim)
    rep_dim = ndims(stack);
end

USE_DEFAULT_NORMALIZATION = 0; % For function 'std', use the default normalization by N-1.

signal_map = abs(orig);
noise_map = std(abs(stack), USE_DEFAULT_NORMALIZATION, rep_dim);
snr_map = signal_map ./ noise_map;

return