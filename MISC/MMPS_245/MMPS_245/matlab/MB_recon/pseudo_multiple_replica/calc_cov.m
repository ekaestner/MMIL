function cov_mtx = calc_cov(dat)
%
% function cov_mtx = calc_cov(dat)
%
% Calculates covariance matrix of complex signals.
% Note this is different from the Matlab built-in cov function in that
% this function normalizes by 2N, whereas the Matlab built-in cov function
% normalizes by N-1 or by N.
% Reference: Philip M. Robson, et al. MRM 60:895-907 (2008).
%
% Input
%   dat     - Each column is a variable, each row is an observation. Dim: [N, nv]
%
% Output
%   cov_mtx - Covariance matrix. Dim: [nv, nv]
%
% (c) Kangrong Zhu      Stanford University     Oct 2014

[N, nv] = size(dat);
cov_mtx = (1/2/N) * ( dat.' * conj(dat) ); % The 1/2 scaling factor accounts for both real and imaginary parts

return