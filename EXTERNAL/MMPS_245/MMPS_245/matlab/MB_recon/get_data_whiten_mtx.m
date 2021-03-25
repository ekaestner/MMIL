function mtx = get_data_whiten_mtx(d, d_type)
%
% Get the data whitening matrix.
%
% Inputs
%   d      - The data to use to form the data whitening matrix.
%            If 'd_type' is 'coil_noise_std', this is a vector containing the coil noise standard deviations. length(d) = nc (i.e. Number of Coils);
%            If 'd_type' is 'psi_mtx', this is the coil noise covariance matrix. Dim: [nc, nc].
%   d_type - Type of input data 'd'.
%            'coil_noise_std': 'd' contains coil noise standard deviation, i.e. sqrt of the diagonal elements in the coil noise covariance matrix.
%            'psi_mtx': 'd' contains coil noise covariance matrix.
%            default: if 'd' is a row or column vector, default to 'coil_noise_std', if 'd' is a square matrix, default to 'psi_mtx'.
%
% Output
%   mtx    - Data whitening matrix. whitened_data = mtx * raw_data. Dim: [nc, nc].
%
% (c) Kangrong Zhu  Stanford University     March 2015

if ~exist('d_type', 'var') || isempty(d_type)
    if isrow(d) || iscolumn(d)
        d_type = 'coil_noise_std';
    else
        if size(d, 1) == size(d, 2)
            d_type = 'psi_mtx';
        else
            error('Input data dimension wrong.');
        end
    end
end

switch d_type
    case 'coil_noise_std' % Approximation using coil noise standard deviation
        mtx = diag(1 ./ d); % Dim: [nc, nc]
        
    case 'psi_mtx' % Use full coil noise covariance matrix
        [V, D] = eig(d);
        mtx = V * diag(1./sqrt(diag(D))) / V; % whitening matrix is (psi_mtx)^(-1/2) = V * (1./sqrt(D)) * inv(V). Dim: [nc, nc]
end

return