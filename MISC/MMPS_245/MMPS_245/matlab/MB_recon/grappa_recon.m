function dat = grappa_recon(us_dat, ker, method)
%
% function dat = grappa_recon(us_dat, ker [, method])
%
% GRAPPA reconstruction.
%
% Inputs:
%   us_dat - Undersampled k-space data, with unacquired lines set to 0.
%            Dim: [n_fe(full), n_pe(full), n_coil].
%   ker    - Interpolation kernel, i.e. the output of the 'grappa_kernel' function.
%            If 'method' is 'imSpace', this is the multiplication kernel in the image space.
%            If 'method' is 'kSpace', this is the convolution kernel in the k-space. 
%            Dim: [FE(source), PE(source), n_coil(source), n_coil(target)].
%   method - The same as the input 'method' in function 'grappa_kernel'.
%            'imSpace'(Default): the reconstruction will be carried out
%                                by multiplication in the image space.
%            'kSpace': the reconstruction will be carried out by 
%                      convolution in the k-space. 
%
% Output:
%   dat    - Reconstructed data, having the same size as the input 'us_dat'.
%            This is the k-space data if 'method' is 'kSpace', and is the
%            image space data if 'method' is 'imSpace'.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

% -- Inputs
if ~exist('method', 'var') || isempty(method)
    method = 'imSpace';
end

if ~(strcmp(method, 'kSpace') || strcmp(method, 'imSpace'))
    error('''method'' must be either ''kSpace'' or ''imSpace''.');
end

% -- Data size
[n_fe, n_pe, n_coil] = size(us_dat);

% -- GRAPPA reconstruction
dat = complex(zeros(n_fe, n_pe, n_coil));

if strcmp(method, 'kSpace')
    for tgt_coil = 1:n_coil
        for src_coil = 1:n_coil
            dat(:, :, tgt_coil) = dat(:, :, tgt_coil) + conv2( ...
                us_dat(:, :, src_coil), ker(:, :, src_coil, tgt_coil), 'same');
        end
    end
end

if strcmp(method, 'imSpace')
    if (size(ker, 1) ~= n_fe) || (size(ker,2) ~= n_pe)
        error('The input kernel size is not consistent with the undersampled k-space data size.');
    end
    
    us_dat = ifft2c(us_dat);                                    % undersampled coil images. The variable name 'us_dat' is kept unchanged to reduce memory use.
    for coil = 1 : n_coil
        dat(:, :, coil) = sum(us_dat .* ker(:, :, :, coil), 3); % reconstructed coil images
    end
end  

return
