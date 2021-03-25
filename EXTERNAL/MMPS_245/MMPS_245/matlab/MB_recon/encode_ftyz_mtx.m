function ft_mtx = encode_ftyz_mtx(msk, ny, nz, add_vcc)
%
% function ft_mtx = encode_ftyz_mtx(msk, ny, nz, add_vcc)
% 
% Calculate Fourier transform (FT) encoding matrix for the y-z plane.
% In y: DFT encoding.
% In z: DFT encoding(CAIPI) or DTFT encoding(MICA).
%
% Inputs
%   msk    - Structure for sampling mask on the ky-omegaz plane. Fields 'add_vcc', 'kz',
%            'omegaz', 'ftz_pha (if not exist, will be calculated using 'kz' and 'omegaz')'
%            and  'ky (used when ny > 1)' are used in this function. See get_ky_omegaz_us_msk.m for details.
%   ny     - Full matrix size in y. If ny <= 1, the DFTy part will not be included in the FT encoding matrix.
%   nz     - Number of simultaneous slices.
%   add_vcc- 0: No virtual coils, use only actual coils; 1: Use both actual and virtual coils; 2: Use only virtual coils.
%
% Output
%   ft_mtx - The FT encoding matix for the y-z plane.
%            Dim: [sample(nsamp*(1+extra_dat_vcc)), Y->Z(=ny*nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))].
%            if add_vcc == 0: extra_dat_vcc = 0, Order in 1st dimension is sample ordering in actual coils.
%            if add_vcc == 1: extra_dat_vcc = 1, Order in 1st dimension is sample ordering in actual coils -> sample ordering in virtual coils.
%            if add_vcc == 2: extra_dat_vcc = 0, Order in 1st dimension is sample ordering in virtual coils.
%
% (c) Kangrong Zhu,     Stanford University     Oct 2013

KEEP_ORIG_SZ             = 1;
NEG_OMEGAZ_ENC           = 1;
SPECIFIED_OMEGAZ_ENC     = 0;
ONLY_ACTUAL_COILS        = 0;
ACTUAL_AND_VIRTUAL_COILS = 1;
ONLY_VIRTUAL_COILS       = 2;

if ~exist('add_vcc', 'var') || isempty(add_vcc)
    add_vcc = ONLY_ACTUAL_COILS;
end

if (ny > 1) && (length(msk(1).ky) ~= length(msk(1).kz))
    error('Number of samples mismatch between msk.ky and msk.kz.');
end
nsamp = length(msk(1).kz);                           % Number of samples on the ky-omegaz plane (i.e. number of echos)

% Include the FTz(DFTz or randomly sampled DTFTz) part in the FT encoding matrix
if ~isfield(msk, 'ftz_pha')
    msk = encode_ftz_pha(msk, nz, SPECIFIED_OMEGAZ_ENC); % msk.ftz_pha is the encoding phase added to each individual slice for each sample. Dim: [nsamp, nz(z indices (-floor(nz/2):1:(ceil(nz/2)-1))])
end
switch add_vcc
    case ONLY_ACTUAL_COILS
        extra_dat_vcc = 0;
    case ACTUAL_AND_VIRTUAL_COILS
        extra_dat_vcc = 1;
        tmp_msk = encode_ftz_pha(msk, nz, NEG_OMEGAZ_ENC);
        msk.ftz_pha = cat(1, msk.ftz_pha, tmp_msk.ftz_pha);
    case ONLY_VIRTUAL_COILS
        extra_dat_vcc = 0;
        msk = encode_ftz_pha(msk, nz, NEG_OMEGAZ_ENC);
end
msk.ftz_pha = ifftshift(msk.ftz_pha, 2);             % z indices become ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))
ftz_mtx = zeros(nsamp*(1+extra_dat_vcc), ny*nz);     % The FTz part in the encoding matrix. Dim: [sample(=nsamp*(1+extra_dat_vcc)), Y->Z(=ny*nz)(z indices are ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))]
for slice = 1 : nz
    ftz_mtx(:, (slice-1)*ny+1 : slice*ny) = repmat( msk.ftz_pha(:, slice), [KEEP_ORIG_SZ, ny]); % Do not divide by sqrt(nz)(Not using orthonormal basis but using the actual encoding phase) so that the signal scaling is consistent between mux=1 and mux>1 in simulation.
end

ft_mtx = ftz_mtx;

% If ny > 1, include the DFTy part in the FT encoding matrix
if ny > 1
    ky_indices = -ny/2 : 1 : ny/2-1;                 % Raw data matrix always arranged top->bottom(-ky->+ky) in ky no matter whether it was acquired top->bottom or bottom->top.
    ky_indices = ky_indices(msk.ky(:)).';
    y_indices = fftshift(0 : ny-1);                  % fftshift is essential because the 0-th indexed PE pixel is in the center along y.
    dfty_mtx = exp(-1i*2*pi*ky_indices*y_indices/ny) ./ sqrt(ny); % Encoding using orthonormal DFT, to be consistent with all the (i)fftc and (i)fft2c used in the recon. Using orthonormal basis here doesn't affect simulation since the ky encoding is the same for mux=1 and mux>1.
    switch add_vcc
        case ONLY_ACTUAL_COILS
            % DO NOTHING
        case ACTUAL_AND_VIRTUAL_COILS
            dfty_mtx = cat(1, dfty_mtx, exp(-1i*2*pi*(-ky_indices)*y_indices/ny) ./ sqrt(ny));
        case ONLY_VIRTUAL_COILS
            dfty_mtx = exp(-1i*2*pi*(-ky_indices)*y_indices/ny) ./ sqrt(ny);
    end
    dfty_mtx = repmat(dfty_mtx, [KEEP_ORIG_SZ, nz]); % Dim: [nsamp*(1+extra_dat_vcc), ny*nz]
    
    ft_mtx = ft_mtx .* dfty_mtx;
end

return