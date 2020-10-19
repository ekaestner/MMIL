function [dat, nt] = mux_encode(dat, msk, nt, noise, eddy)
% function [dat, nt] = mux_encode(dat, msk, nt, [noise, eddy])
%
% Encode single-slice k-space data into slice-multiplexed k-space data.
%
% Inputs
%   dat    - Fully sampled single-slice K-space data. Dim: [Kx(=nx), Ky(=ny, full matrix size in y), Echo(=nec), Slice(=nsl),
%            Coil(=nc), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))]
%   msk    - Undersampling mask on the ky-omegaz plane. Fields 'ky','kz','omegaz'
%            are used in this function. See get_ky_omegaz_us_msk.m for details.
%   nt     - Number of time points to synthesize. Can be any positive integer
%            if length(msk)==1; Will be set equal to length(msk) if length(msk)>1.
%   noise  - If not empty, add this noise to the encoded data. Dim: length(noise(:)) = length(dat(:)) = nx*nsamp*nec*nsl*nc*nt.
%   eddy   - If not empty: Contains phase terms the eddy current effects cause.
%              Eddy current effect will be included in the encoding process.
%              xKyOmegaz-data-with-eddy-current = xKyOmegaz-data-without-eddy-current .* eddy.
%              If funtcion 'epi_pha_correct' is used to calculate the phase terms,
%              'eddy' should be reshaped from the output 'pha_flt'.
%              Dim: [X(=nx), PE(=nsamp), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))), Coil(=nc)].
%            If empty: eddy current effect will not be included.
%
% Outputs
%   dat    - Encoded slice-multiplexed k-space data. Dim: [Kx(=nx), Ky(=nsamp, number samples on ky-omegaz plane),
%            Echo(=nec), Slice(=nsl), Coil(=nc), Time(=nt)].
%   nt     - Number of time points in the synthesized data set.
%
% (c) Kangrong Zhu      Stanford University     Oct 2013

add_noise = (exist('noise', 'var') && ~isempty(noise)); % True: Add noise to encoded slice-multiplexed k-space data
include_eddy = (exist('eddy', 'var') && ~isempty(eddy)); % True: Include the eddy current effect in the encoding process

KEEP_ORIG_SZ = 1;
FE_DIM = 1;

nmsk = length(msk);
if (nmsk > 1) && (nt ~= nmsk)
    nt = nmsk; % nt equal to nmsk if nmsk>1; nt can be any number if nmsk==1
end
nsamp = length(msk(1).kz);

[nx, ny, nec, nsl, nc, nz] = size(dat);
dat = ifft2c(dat); % Dim: [X(=nx), Y(=ny), Echo, Slice, Coil, SimultaneousSlice Z(=nz)]
add_vcc = false; % Need to be false to synthesize acquired data

res = zeros(nsamp, nx, nec, nsl, nc, nmsk); % Dim: [Sample(=nsamp), X, Echo, Slice, Coil, Undersampling Mask]
for msk_idx = 1 : nmsk
    ftyz_mtx = encode_ftyz_mtx(msk(msk_idx), ny, nz, add_vcc); % Dim: [Sample(=nsamp*(1+add_vcc)), Y->Z(=ny*nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))))]
    for x = 1 : nx
        for ec = 1 : nec
            for sl = 1 : nsl
                for c = 1 : nc
                    encode_mtx = ftyz_mtx;
                    if include_eddy
                        eddy_mtx = eddy(x, :, ec, sl, :, c); % Dim: [X(=1), Sample(=nsamp), Echo(=1), Slice(=1), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))), Coil(=1)]
                        eddy_mtx = permute(eddy_mtx, [2,7,5,1,3,4,6]); % Dim: [Sample(=nsamp), Y(=1), SimultaneousSlice Z(=nz)]
                        eddy_mtx = repmat(eddy_mtx, [1, ny, 1]); % Dim: [Sample(=nsamp), Y(=ny), SimultaneousSlice Z(=nz)]
                        eddy_mtx = reshape(eddy_mtx, [nsamp, ny*nz]); % Dim: [Sample(=nsamp), Y->Z(=ny*nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))))]
                        
                        encode_mtx = encode_mtx .* eddy_mtx;
                    end
                    
                    dat_raw = reshape(dat(x, :, ec, sl, c, :), [ny*nz, 1]); % Dim: [Y->Z(=ny*nz), 1]
                    res(:, x, ec, sl, c, msk_idx) = encode_mtx * dat_raw; % Right side dim: [Sample(=nsamp), 1]
                end
            end
        end
    end
end
dat = permute(res, [2,1,3,4,5,6]); % Dim: [X(=nx), Sample(=nsamp), Echo, Slice, Coil, Undersampling Mask(=nmsk)]
dat = fftc(dat, FE_DIM); % Dim: [Kx(=nx), Sample(=nsamp), Echo, Slice, Coil, Undersampling Mask(=nmsk)]
if (nmsk == 1) && (nt > 1)
    dat = repmat(dat, [KEEP_ORIG_SZ, KEEP_ORIG_SZ, KEEP_ORIG_SZ, KEEP_ORIG_SZ, KEEP_ORIG_SZ, nt]); % Dim: [Kx(=nx), Sample(=nsamp), Echo, Slice, Coil, Time]
end

if add_noise
    dat = dat + reshape(noise, [nx, nsamp, nec, nsl, nc, nt]);
end

return