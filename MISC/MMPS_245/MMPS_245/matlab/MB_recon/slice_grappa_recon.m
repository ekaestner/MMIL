function dat = slice_grappa_recon(dat_mux, ker, us_msk, method)
%
% function dat = slice_grappa_recon(dat_mux, ker, us_msk, [method='imSpace'])
%
% Slice-GRAPPA reconstruction.
% Reference: Kawin Setsompop, et al. MRM 2012;67(5):1210-24.
%
% Inputs
%   dat_mux - Slice-multiplexed k-space data. Dim: [Kx(=nx), Ky(=ny), Coil(=nc)].
%   ker     - Interpolation kernel, i.e. the output of the 'slice_grappa_kernel' function.
%             If 'method' is 'imSpace', this is the multiplication kernel in the image space.
%             If 'method' is 'kSpace', this is the convolution kernel in the k-space.
%             Dim: [FE(source), PE(source), nc(source), nc(target), nz(target, z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%   us_msk  - Undersample mask with fields 'ky', 'kz', 'omegaz'.
%             See get_ky_omegaz_us_msk.m for details.
%   method  - How the reconstruction will be carried out.
%             'imSpace': Multiplication in image space.
%             'kSpace':  Convolution in k-space.
%
% Output
%   dat     - Reconstructed single-slice k-space data. Dim: [Kx(=nx), Ky(=ny), Coil(=nc), SimultaneousSlice Z(=nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2014

SPECIFIED_OMEGAZ_ENC = 0;

% -- Inputs
if ~exist('method', 'var') || isempty(method)
    method = 'imSpace';
end

if ~(strcmp(method, 'kSpace') || strcmp(method, 'imSpace'))
    error('''method'' must be either ''kSpace'' or ''imSpace''.');
end

% -- Data size
[nx, ny, nc] = size(dat_mux);
nz = size(ker, 5);

% -- Slice-GRAPPA reconstruction
dat = zeros(nx, ny, nc, nz);

if strcmp(method, 'kSpace')
    for z = 1 : nz
        for tgt_coil = 1:nc
            for src_coil = 1:nc
                dat(:, :, tgt_coil, z) = dat(:, :, tgt_coil, z) + conv2(dat_mux(:, :, src_coil), ker(:, :, src_coil, tgt_coil, z), 'same');
            end
        end
    end
end

if strcmp(method, 'imSpace')
    if (size(ker, 1) ~= nx) || (size(ker,2) ~= ny)
        error('The input kernel size is not consistent with the slice-multiplexed k-space data size.');
    end
    
    dat_mux = ifft2c(dat_mux);                                                     % Slice-multiplexed coil images. The variable name 'dat_mux' is kept unchanged to reduce memory use.
    for z = 1 : nz
        for tgt_coil = 1 : nc
            dat(:, :, tgt_coil, z) = sum(dat_mux .* ker(:, :, :, tgt_coil, z), 3); % Reconstructed single-slice coil images
        end
    end
end

% -- Remove FTz encoding phase from each slice
us_msk = encode_ftz_pha(us_msk, nz, SPECIFIED_OMEGAZ_ENC);                         % us_msk.ftz_pha is the encoding phase added to each individual slice for each sample. Dim: [nsamp(=ny), nz(z indices -floor(nz/2):1:(ceil(nz/2)-1))]
if strcmp(method, 'imSpace')
    dat = fft2c(dat);
end
for z = 1 : nz
    pha_mtx = repmat(us_msk.ftz_pha(:, z).', [nx, 1]);
    for coil = 1 : nc
        dat(:, :, coil, z) = dat(:, :, coil, z) .* conj(pha_mtx);
    end
end

return