function [psf, max_intr] = us_msk_psf(us_msk, ny, nz, kydir, pdf)
%
% function [psf, max_intr] = us_msk_psf(us_msk, ny, [nz=number of kz lines], [kydir=BOTTOM_UP], [pdf=(1/R)*ones(ny, nkz)])
%
% Calculates spatial point spread function and maximum interference for slice-multiplexed acquisition (Assuming a delta function is placed at the center of the middle slice in the spatial domain)
%
% Inputs
%   us_msk   - Ky-omegaz undersample mask with fields 'ky', 'kz' and 'omegaz'. See get_ky_omegaz_us_msk.m for details.
%   ny       - Full prescribed matrix size in y (phase encoding direction).
%   nz       - Number of simultaneous slices.
%   kydir    - Ky traversal direction. Value equal to TOP_DOWN, CENTER_OUT or BOTTOM_UP.
%   pdf      - Probably density function on ky-omegaz plane. Dim: [ny, nkz(=length(us_msk(1).omegaz))].
%
% Output
%   psf      - Point spread function, assuming a delta function is placed
%              at the center of the middle slice. If a delta function is
%              placed at the center of another slice, the point spread function will be circularly shifted.
%   max_intr - Maximum interference in the spatial point spread function.
%
% (c) Kangrong Zhu  Stanford University     June 2014

BOTTOM_UP  = 2;
if ~exist('kydir', 'var')
    kydir = BOTTOM_UP;
end

% Matrix representation of undersample mask
mtx = us_msk_get_mtx(us_msk, ny, kydir, []);
nkz = length(mtx(1).omegaz); % Number of kz lines

if ~exist('nz', 'var') || isempty(nz)
    nz = nkz;
end
if ~exist('pdf', 'var') || isempty(pdf)
    R = (ny/length(us_msk(1).ky)) * nkz;
    pdf = (1/R) * ones(ny, nkz);
end

% PSF
nmsk = length(us_msk);
max_intr = zeros(1, nmsk);
psf = zeros(ny, nz, nmsk);
for msk_idx = 1 : nmsk
    
    % Raw PSF with nkz lines along z
    msk = (mtx(msk_idx).mtx >= 1);
    psf_raw = fftshift(ifft2(ifftshift(msk ./ pdf))); % Point Spread Function when the middle slice contains a delta function signal in image domain
    
    % Crop or pad PSF to nz lines along z
    if nz == nkz
        psf_nz = psf_raw;
    end
    if nz < nkz
        psf_nz = crop(psf_raw, [ny, nz]);
    end
    if nz > nkz
        npad_left = ceil((nz-nkz)/2);
        npad_right = (nz-nkz) - npad_left;
        psf_nz = psf_raw;
        % Repeat PSF to pad left side
        for pad_idx = 1 : ceil(npad_left/nkz)
            if npad_left >= pad_idx*nkz               % Need to use all nkz lines in z for this padding
                psf_nz = [psf_raw, psf_nz];
            else                                      % Last padding, will only need some of the lines in z
                npad_this_time = npad_left - (pad_idx-1)*nkz;
                psf_nz = [psf_raw(:, nkz-npad_this_time+1:nkz), psf_nz];
            end
        end
        % Repeat PSF to pad right side
        for pad_idx = 1 : ceil(npad_right/nkz)
            if npad_right >= pad_idx*nkz              % Need all nkz lines in z
                psf_nz = [psf_nz, psf_raw];
            else                                      % Need only some of the lines in z
                npad_this_time = npad_right - (pad_idx-1)*nkz;
                psf_nz = [psf_nz, psf_raw(:, 1:npad_this_time)];
            end
        end
    end
    
    % Output
    psf(:, :, msk_idx) = psf_nz;
    
    if nz <= nkz                                      % The spatial domain PSF was not repeated
        psf_alias = psf_nz;
    else                                              % nz > nkz, the spatial domain PSF was repeated
        psf_alias = crop(psf_nz, [ny, nkz]);          % Only look at one instance
    end
    psf_alias = ifftshift(psf_alias);
    if nkz > 1
        psf_alias = psf_alias(:, 2:end);              % Exclude the slice which contains the delta function in image domain so that any inplane aliasing within this slice is excluded.
    else % nkz = 1, size(psf_alias) = [ny, 1]
        psf_alias = psf_alias(2 : ny);
    end
    max_intr(msk_idx) = max(abs(psf_alias(:)));
end

return