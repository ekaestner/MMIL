function m = calc_coil_noise_cov(d, c)
%
% function m = calc_coil_noise_cov(d, [c=ndims(d)])
%
% Calculates coil noise covariance matrix from noise data matrix.
%
% Inputs
%   d - Noise data matrix, with the c-th dimension being the coil dimension.
%   c - Coil dimension in input 'd'.
%
% Output
%   m - Coil noise covariance matrix. Dim: [nc, nc], where nc = size(d, c).
%
% (c) Kangrong Zhu  Stanford University     June 2014

if ~exist('c', 'var') || isempty(c)
    c = ndims(d);
end
s = size(d);

nc = s(c);
p = [1:c-1, c+1:length(s), c]; % Order for matrix permutation
d = permute(d, p); % Permute to make coil dimension the last dimension
d = reshape(d, [], nc);
m = calc_cov(d);

return