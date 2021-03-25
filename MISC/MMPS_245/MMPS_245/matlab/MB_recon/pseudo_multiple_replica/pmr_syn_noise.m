function noise_mtx = pmr_syn_noise(psi_mtx, noise_mtx_sz, coil_dim)
%
% function noise_mtx = pmr_syn_noise(psi_mtx, noise_mtx_sz, [coil_dim=length(noise_mtx_sz)])
%
% Pseudo Multiple Replica: Synthesize scaled and correlated noise records
% for each channel of the phased-array from uncorrelated
% Gaussian-distributed white noise with unity standard deviation.
% Reference: Philip M. Robson, et al. MRM 60:895-907 (2008).
%
% Inputs
%   psi_mtx      - Coil noise covariance matrix. Dim: [nc, nc]
%   noise_mtx_sz - Size of the desired noise matrix.
%   coil_dim     - Coil dimension in the desired noise matrix.
%
% Output
%   noise_mtx    - Noise matrix.
%
% Example
%   noise_mtx = pmr_syn_noise(eye(32), [80, 80, 32, 10], 3);
%   returns matrix with IID noise with a size of [80, 80, 32, 10] for 32 coils.
%
% (c) Kangrong Zhu      Stanford University     Oct 2014

nd = length(noise_mtx_sz); % Number of dimensions in the desired noise matrix
if ~exist('coil_dim', 'var') || isempty(coil_dim)
    coil_dim = nd;
end
nc = noise_mtx_sz(coil_dim); % Number of coils
non_coil_dim_sz = noise_mtx_sz([1 : coil_dim-1, coil_dim+1 : nd]); % Size for non-coil dimensions
np = prod(non_coil_dim_sz); % Total number of noise data points to generate
res_permute_order = [2 : coil_dim, 1, coil_dim+1 : nd];

[V, D] = eig(psi_mtx);
matrix_square_root = V * sqrt(D) * inv(V); % Dim: [nc, nc]

raw_mtx_sz = [nc, np];
noise_mtx = randn(raw_mtx_sz) + 1i * randn(raw_mtx_sz); % Uncorrelated Gaussian-distributed white noise with unity standard deviation. Dim: [nc, np]
noise_mtx = matrix_square_root * noise_mtx; % Scaled and correlated noise records for each channel of the phased-array. Dim: [nc, np]
noise_mtx = reshape(noise_mtx, [nc, non_coil_dim_sz]); % Dim: noise_mtx_sz, except the coil dimension is moved to the 1st dimension
noise_mtx = permute(noise_mtx, res_permute_order); % Dim: noise_mtx_sz

return