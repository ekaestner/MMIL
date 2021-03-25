function [sense_rsnr, sense_noise, sense_noise_ref, nmsk, nsamp_ref] = calc_sense_rsnr(p, nt, dat_cal)
%
% function [sense_rsnr, sense_noise, sense_noise_ref, nmsk, nsamp_ref] = calc_sense_rsnr(p, nt, dat_cal)
%
% Calculates retained SNR and noise maps using the hybrid-space SENSE reconstruction model, by calling function mux_calc_im_noise.m.
%
% Inputs
%   p               - Parameter structure. See mux_epi_params.m for details.
%                     In addition to fields defined in mux_epi_params.m, here p may also contain fields:
%                     sense_noise_ref - Reference noise maps for calculating the retained SNR maps. Dimension same as the output 'sense_noise_ref'.
%                                       If this field doesn't exist or is empty, will calculate the reference noise maps inside this function.
%                     nsamp_ref       - Number of acquired points on the ky-omegaz plane in the reference scan. Must exist and be nonempty if p.sense_noise_ref exists and is nonempty.
%                     us_msk          - Undersampling mask. See get_ky_omegaz_us_msk.m for details.
%                                       If this field doesn't exist or is empty, will calculate the undersampling mask inside this function.
%   nt              - Number of time points in the slice-accelerated acquisition, used to determine how many undersampling masks there are.
%   dat_cal         - Calibration data. Dim: [Kx(=nxi, before ramp sampling correction), Ky(=nky, size of raw k-space data loaded from p-file, inplane acceleration interleaved with zeros, partial ky not zero-padded), Echo(=nec), Slice(=nsl), Coil(=nci, before coil compression), Time(i.e. Kz, = p.mux)]
%
% Outputs
%   sense_rsnr      - Retained SNR maps. Retained SNR maps are accurate and not scaled.
%                     Dim: [X(=p.nx_pres, if p.zpad_image=0; =p.zpad_size(1), if p.zpad_image=1), Y(=p.ny_pres, if p.zpad_image=0; =p.zpad_size(2), if p.zpad_image=1), Echo(=nec), Slice(=nsl*nz), SigularDimension(=1), UndersamplingMask(=nmsk)].
%   sense_noise     - Noise maps. The noise maps may be scaled when a noise data file is used to calculate coil noise covariance because the data in the noise data file are scaled.
%                     Dim: Same as sense_rsnr.
%   sense_noise_ref - Noise maps of a reference scan (a mux1 scan (if p.sense_rsnr_ref_type is 'mux1') or a slice-DFT-encoded calibration scan (if p.sense_rsnr_ref_type is 'muxcal')) that has the same inplane acceleration as the mux scan.
%                     The noise maps may be scaled when a noise data file is used to calculate coil noise covariance.
%                     Dim: [X(=size(sense_rsnr,1)), Y(=size(sense_rsnr,2)), Echo(=nec), Slice(=nsl*nz), SigularDimension(=1), UndersamplingMask(=1)].
%   nmsk            - Number of unique undersampling masks simulated.
%   nsamp_ref       - Number of acquired points on the ky-omegaz plane in the reference scan (a mux1 scan or a slice-DFT-encoded calibration scan).
%
% (c) Kangrong Zhu  Stanford University     March 2015

%% Prepare inputs to function mux_calc_im_noise.m
[nxi, nky, nec, nsl, nci, nkz] = size(dat_cal); % nky - Size in Ky of the acquired k-space data (raw k-space data loaded from p-file. Inplane acceleration interleaved with zeros, partial ky not zero-padded); nci: Number of coils before coil compression

if nec ~= 1
    error('Number of echoes must be 1 in the calibration data.');
end

% Undersampling mask
if isfield(p, 'us_msk') && ~isempty(p.us_msk)
    us_msk = p.us_msk;
else
    us_msk = get_ky_omegaz_us_msk(p, nky, nt, true);
end

% Sensitivity maps
if ~isfield(p, 'smap') || isempty(p.smap)
    error('No sensitivity maps');
end

% Coil noise covariance matrix
if ~isfield(p, 'psi_mtx') || isempty(p.psi_mtx)
    error('No coil noise covariance matrix.');
end

% Data whitening type
if ~isfield(p, 'whiten_type') || isempty(p.whiten_type)
    error('No data whitening type.');
end

% Ramp sampling correction filter
if p.vrgf
    if ~isfield(p, 'pmr_vrgf') || isempty(p.pmr_vrgf)
        error('No vrgf file.');
    end
    ramp_flt = rawload_vrgf(nxi, p.nx_pres, p.pmr_vrgf);
else
    ramp_flt = [];
end

% Coil compression matrix
if p.coil_compress
    if ~isfield(p, 'ccmtx') || isempty(p.ccmtx)
        error('No coil compression matrix.');
    end
else
    p.ccmtx = [];
end

% Image noise type
if ~isfield(p, 'sense_rsnr_im_noise_type') || isempty(p.sense_rsnr_im_noise_type)
    error('No image noise type specified.');
end

% Coil noise standard deviation
if ~isfield(p, 'coil_noise_std') || isempty(p.coil_noise_std) || any(isnan(p.coil_noise_std))
    if strcmp(p.whiten_type, 'coil_noise_std')
        error('No coil noise standard deviation.');
    else
        p.coil_noise_std = [];
    end
end

% Debugging flag
if ~isfield(p, 'debug') || isempty(p.debug)
    p.debug = false;
end

%% Calculate image noise of the SENSE-reconstructed slice-multiplexed images
sense_noise = mux_calc_im_noise(us_msk, p.smap, p.psi_mtx, p.whiten_type, p.pha_coe_default, ramp_flt, p.ccmtx, p.sense_rsnr_im_noise_type, p.sense_rsnr_propagate_type, p.coil_noise_std, p.debug); % Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec=1), Slice*SimultaneousSlice Z(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)),  MaskType(=nmsk)]

%% Output results
nmsk = length(us_msk);
nz = size(p.smap, 5);
sense_noise = reshape(sense_noise, [p.nx_pres, p.ny_pres, nec, nsl*nz, 1, nmsk]);
if p.zpad_image
    sense_noise = ifft2c(zpad(fft2c(sense_noise), [p.zpad_size(1), p.zpad_size(2), nec, nsl*nz, 1, nmsk]));
end

%% Calculate image noise of the reference scan
nsamp = length(us_msk(1).ky); % Number of k-space samples in the slice-multiplexed acquisition
if isfield(p, 'sense_noise_ref') && ~isempty(p.sense_noise_ref)
    if ~isfield(p, 'nsamp_ref') || isempty(p.nsamp_ref)
        error('Must input nsamp_ref with sense_noise_ref.');
    else
        nsamp_ref = p.nsamp_ref;
        sense_noise_ref = p.sense_noise_ref;
    end
else
    switch p.sense_rsnr_ref_type
        case 'muxcal' % Use DFT encoded calibration acquisition as reference
            if p.debug
                fprintf(' Calculating image noise of mux calibration with matching inplane acceleration...\n');
            end
            
            % - Undersampling mask
            us_msk_muxcal.omegaz = encode_dftz_omegaz(p.cap_fov_shift_cal);
            nkz = length(us_msk_muxcal.omegaz);
            us_msk_muxcal.ky = repmat(us_msk(1).ky(:), [nkz, 1]);
            us_msk_muxcal.kz = reshape( repmat(1 : nkz, [nsamp, 1]), [nsamp*nkz, 1]);
            
            % - Calcualte image noise
            sense_noise_ref = mux_calc_im_noise(us_msk_muxcal, p.smap, p.psi_mtx, p.whiten_type, p.pha_coe_default, ramp_flt, p.ccmtx, p.sense_rsnr_im_noise_type, p.sense_rsnr_propagate_type, p.coil_noise_std, p.debug); % Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec=1), Slice*SimultaneousSlice Z(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))]
            
            % - Number of samples acquired in the reference scan
            nsamp_ref = length(us_msk_muxcal.ky);
            
        case 'mux1' % Use MUX1 acquisition as reference
            if p.debug
                fprintf(' Calculating image noise of mux1 acquisition with matching inplane_acceleration...\n');
            end
            
            % - Set up inputs to function mux_calc_im_noise
            
            % Parameter
            p_mux1 = pmr_params_set_mux1(p);
            p_mux1.pmr_slices = repmat(p_mux1.pmr_slices(:), [nz, 1]); % Dim: [nsl*nz, 1]. Let the simulation reload the same data from ref.dat for Nyquist ghosting correction for each slice
            p_mux1.num_coils = nci;
            p_mux1.num_mux_cycle = 1;
            
            % Undersmpling mask
            us_msk_mux1 = get_ky_omegaz_us_msk(p_mux1, nky, nt, true); % The CAIPI-type us_msk_mux1 only contains one mask, length(us_msk_mux1) = 1
            
            % Coil compression matrix
            if p_mux1.coil_compress
                ccmtx_mux1 = repmat(p.ccmtx, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nz]);
            else
                ccmtx_mux1 = [];
            end
            
            % Sensitivity maps
            sz = size(p.smap);
            smap_mux1 = reshape(fftshift(p.smap, 5), [sz(1), sz(2), sz(3), sz(4)*sz(5), 1, sz(6)]);
            
            % Nyquist ghosting correction coefficients
            pha_coe_mux1 = repmat(p.pha_coe_default, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nz, p.KEEP_ORIG_SZ]);
            
            % - Calculate image noise
            sense_noise_ref = mux_calc_im_noise(us_msk_mux1, smap_mux1, p.psi_mtx, p.whiten_type, pha_coe_mux1, ramp_flt, ccmtx_mux1, p.sense_rsnr_im_noise_type, p.sense_rsnr_propagate_type, p.coil_noise_std, p.debug); % Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec=1), Slice*SimultaneousSlice Z(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))]
            
            % - Number of samples acquired in the reference scan
            nsamp_ref = length(us_msk_mux1.ky);
    end
    
    % Zero pad reconstructed images
    if p.zpad_image
        sense_noise_ref = ifft2c(zpad(fft2c(sense_noise_ref), [p.zpad_size(1), p.zpad_size(2), nec, nsl*nz]));
    end
end

%% Calculate and output retained SNR maps

% Acceleration factor w.r.t. reference scan
R = nsamp_ref / nsamp;

% Retained SNR
sz = size(sense_noise);
sense_rsnr = zeros(sz(1), sz(2), nec, nsl*nz, 1, nmsk);
for msk_idx = 1 : nmsk
    sense_rsnr(:, :, :, :, 1, msk_idx) = sense_noise_ref ./ sense_noise(:, :, :, :, 1, msk_idx) .* sqrt(R);
end

%% If not in debugging mode, don't output reference noise maps
if ~p.debug
    sense_noise_ref = [];
end

return