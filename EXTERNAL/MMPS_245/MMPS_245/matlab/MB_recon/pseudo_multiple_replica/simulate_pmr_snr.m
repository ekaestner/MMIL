function [pmr_rsnr, pmr_noise, pmr_snr, pmr_noise_ref, nmsk, dat_mux, ccmtx, smap, noise, im_sing, pmr_res] = simulate_pmr_snr(p, nt, dat_cal)
%
% function [pmr_rsnr, pmr_noise, pmr_snr, pmr_noise_ref, nmsk] = simulate_pmr_snr(p, nt, dat_cal)
%
% Calculates retained SNR, noise and SNR maps using pseudo-multiple replica simulation, by calling function pmr_sim_mux_seq.m.
%
% Inputs
%   p           - Parameter structure. See mux_epi_params.m for details.
%                 In addition to fields defined in mux_epi_params.m, here p may also contain fields:
%                 pmr_noise_ref - Reference noise maps for calculating the retained SNR maps. Dimension same as the output 'pmr_noise_ref'.
%                                 If this field doesn't exist or is empty, will calculate the reference noise maps inside this function.
%                 us_msk        - Undersampling mask. See get_ky_omegaz_us_msk.m for details.
%                                 If this field doesn't exist or is empty, will calculate the undersampling mask inside this function.
%                 pmr_recon_sing - True: Reconstruct original single-slice images.
%   nt          - Number of time points in the slice-accelerated acquisition, used to determine how many undersampling masks there are.
%   dat_cal     - Calibration data. Dim: [Kx(=nxi, before ramp sampling correction), Ky(=nky, size of raw k-space data loaded from p-file, inplane acceleration interleaved with zeros, partial ky not zero-padded), Echo(=nec), Slice(=nsl), Coil(=nc), Time(i.e. Kz, = p.mux)]
%
% Outputs
%   pmr_rsnr    - Retained SNR maps. Retained SNR maps are accurate and not scaled.
%                 Dim: [X(=p.nx_pres, if p.zpad_image=0; =p.zpad_size(1), if p.zpad_image=1), Y(=p.ny_pres, if p.zpad_image=0; =p.zpad_size(2), if p.zpad_image=1), Echo(=nec), Slice(=nsl*nz), SigularDimension(=1), UndersamplingMask(=nmsk)].
%   pmr_noise   - Noise maps. Noise maps may be scaled when a noise data file is used to calculate coil noise covariance because the data in the noise data file are scaled.
%                 Dim: Same as pmr_rsnr.
%   pmr_snr     - SNR maps. SNR maps may be scaled when a noise data file is used to calculate coil noise covariance.
%                 Dim: Same as pmr_rsnr.
%   pmr_noise_ref - Noise maps for a MUX1 scan (if p.pmr_ref_type is 'mux1') or for a slice-DFT-encoded calibration scan (if p.pmr_ref_type is 'muxcal') that has the same inplane acceleration as the mux scan.
%                 This was used as the reference noise maps for calculating retained SNR maps.
%                 Noise maps may be scaled when a noise data file is used to calculate coil noise covariance.
%                 Dim: [X(=size(pmr_rsnr,1)), Y(=size(pmr_rsnr,2)), Echo(=nec), Slice(=nsl*nz), SigularDimension(=1), UndersamplingMask(=1)].
%   nmsk        - Number of unique undersampling masks simulated.
%   Output the following to assist the simulation of different undersampling patterns.
%     dat_mux   - Synthesized slice-multiplexed data. Nonempty only when p.debug is true. Fields include:
%                 dat_mux(i).dat - Corresponding to the i-th undersampling mask. The 1st time point contains no added noise, the 2nd time point contains added noise that conforms to the input coil noise covariance matrix. Dim: [Kx(=p.nx_pres), PE(=nsamp), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=2)].
%     ccmtx     - Coil compression matrices. Dim: [Coil(=p.num_coils, before compression), Coil(=p.num_vcoils, after compression), X(=1 if p.cc_method is 'single'; =p.nx_pres if p.cc_method is 'geometric'), Echo(=nec), Slice(=nsl)]
%     smap      - Sensitivity maps. Nonempty only when reconstruction method is SENSE. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Coil(=p.num_coils, after coil compression when applicable), SetOfSensitivityMaps(=p.espirit_nmap if p.smap_type is 'espirit'; =1 if p.smap_type is 'coil_over_sos')].
%     noise     - Noise added to the synthesized slice-multiplexed data. Dim: [Kx(=nxi), PE(=nsamp), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=1+p.pmr_num_rep)].
%     im_sing   - SOS-coil-combined reconstructed original single-slice images. Nonempty only when p.pmr_recon_sing is true. Dim: [X(=p.nx_pres), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))].
%     pmr_res   - Result structure of the PMR simulation, corresponding to the output sense/oneD_grappa/slice_grappa/split_slice_grappa of the function pmr_sim_mux_seq.
%
% (c) Kangrong Zhu      Stanford University     Dec 2014

OUTPUT_IDX_DAT_MUX = 6;
OUTPUT_IDX_PMR_RES = 11;

%% Prepare inputs to function pmr_sim_mux_seq.m
% Size
[nxi, nky, nec, nsl, nc, nkz] = size(dat_cal); % nxi: Size in Kx before ramp sampling correction; nky: Size in Ky of the acquired k-space data (raw k-space data loaded from p-file. Inplane acceleration interleaved with zeros, partial ky not zero-padded); nc: Number of physical coils.
nz = p.mux;

% Parameter structure
p.sim_sense = false;
p.sim_1d_grappa = false;
p.sim_slice_grappa = false;
p.sim_split_slice_grappa = false;
switch p.mux_recon_method
    case 'sense'
        p.sim_sense = true;
    case '1Dgrappa'
        p.sim_1d_grappa = true;
    case 'slice-grappa'
        p.sim_slice_grappa = true;
    case 'split-slice-grappa'
        p.sim_split_slice_grappa = true;
end
if strcmp(p.pmr_ref_type, 'muxcal')
    p.pmr_sim_cal = true;
else
    p.pmr_sim_cal = false;
end
if ~isfield(p, 'pmr_recon_sing')
    p.pmr_recon_sing = false;
end

% Undersampling mask
if isfield(p, 'us_msk') && ~isempty(p.us_msk)
    us_msk = p.us_msk;
    p = rmfield(p, 'us_msk');
else
    us_msk = get_ky_omegaz_us_msk(p, nky, nt, true);
end

% Eddy current effects
eddy = []; % TODO: Not considering eddy current effects for now

% Noise to add to the synthesized slice-multiplexed data
if isfield(p, 'noise')
    noise = p.noise;
    p = rmfield(p, 'noise');
else
    noise = [];
end

% Original single-slice data
if isfield(p, 'im_sing')
    im_sing = p.im_sing;
    p = rmfield(p, 'im_sing');
else
    im_sing = [];
end

% When using muxcal as reference, do not need to simulate repeated calibration data if reference noise map is already input
if isfield(p, 'pmr_noise_ref') && ~isempty(p.pmr_noise_ref) && strcmp(p.pmr_ref_type, 'muxcal') && p.pmr_sim_cal
    p.pmr_sim_cal = false;
end

%% Conduct pseudo-multiple replica simulation
if nargout >= OUTPUT_IDX_DAT_MUX % For simulating different undersampling patterns
    [sense, oneD_grappa, slice_grappa, split_slice_grappa, pmr_cal, dat_mux, ccmtx, smap, noise, im_sing] = pmr_sim_mux_seq(p, us_msk, dat_cal, eddy, noise, im_sing);
    p.ccmtx = ccmtx;
    p.smap = smap;
else
    [sense, oneD_grappa, slice_grappa, split_slice_grappa, pmr_cal] = pmr_sim_mux_seq(p, us_msk, dat_cal, eddy, noise, im_sing);
end

%% Output results
nmsk = length(us_msk);
switch p.mux_recon_method
    case 'sense'
        pmr_res = sense;
    case '1Dgrappa'
        pmr_res = oneD_grappa;
    case 'slice-grappa'
        pmr_res = slice_grappa;
    case 'split-slice-grappa'
        pmr_res = split_slice_grappa;
end
sz = size(pmr_res(1).snr);
if length(sz) < 2
    sz(end+1 : 2) = 1; % Make sure sz(1), sz(2)) all exist
end
pmr_snr = zeros(sz(1), sz(2), nec, nsl*nz, 1, nmsk); % Dim: X(=p.nx_pres, if p.zpad_image=0; =p.zpad_size(1), if p.zpad_image=1), Y(=p.ny_pres, if p.zpad_image=0; =p.zpad_size(2), if p.zpad_image=1), Echo(=nec), Slice(=nsl*nz), SigularDimension(=1), UndersamplingMask(=nmsk)].
pmr_noise = zeros(sz(1), sz(2), nec, nsl*nz, 1, nmsk); % Dim: X(=p.nx_pres, if p.zpad_image=0; =p.zpad_size(1), if p.zpad_image=1), Y(=p.ny_pres, if p.zpad_image=0; =p.zpad_size(2), if p.zpad_image=1), Echo(=nec), Slice(=nsl*nz), SigularDimension(=1), UndersamplingMask(=nmsk)].
for msk_idx = 1 : nmsk
    pmr_snr(:, :, :, :, 1, msk_idx) = pmr_res(msk_idx).snr;
    pmr_noise(:, :, :, :, 1, msk_idx) = pmr_res(msk_idx).noise;
end
if nargout < OUTPUT_IDX_PMR_RES
    clear pmr_res;
end

%% Pseudo-multiple replica simulation for reference scan
if isfield(p, 'pmr_noise_ref') && ~isempty(p.pmr_noise_ref)
    pmr_noise_ref = p.pmr_noise_ref;
else
    switch p.pmr_ref_type
        case 'muxcal' % Use DFT encoded calibration acquisition as reference
            pmr_noise_ref = pmr_cal.noise;
            
        case 'mux1' % Use MUX1 acquisition as reference
            if p.debug
                fprintf(' Simulating multiple replica of mux1 acquisition...\n');
            end
            
            % - Set up inputs to function pmr_sim_mux_seq
            % Parameter
            p_mux1 = pmr_params_set_mux1(p);
            p_mux1.pmr_slices = repmat(p_mux1.pmr_slices(:), [nz, 1]); % Dim: [nsl*nz, 1]. Let the simulation reload the same data from ref.dat for Nyquist ghosting correction for each slice
            
            % Undersampling mask
            us_msk_mux1 = get_ky_omegaz_us_msk(p_mux1, nky, nt, true); % The CAIPI-type us_msk_mux1 only contains one mask, length(us_msk_mux1) = 1
            
            % Coil compression matrix
            if p_mux1.coil_compress
                p_mux1.ccmtx = repmat(p.ccmtx, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, nz]);
            end
            
            % Sensitivity maps
            if p_mux1.sim_sense
                sz = size(p.smap);
                p_mux1.smap = reshape(fftshift(p.smap, 5), [sz(1), sz(2), sz(3), sz(4)*sz(5), 1, sz(6)]);
            end
            
            % Single-slice data
            dat_cal_mux1 = mux_dftz(dat_cal, p.T_DIM, p.cap_fov_shift_cal, nz, 'decode'); % Dim: [Kx(=nxi), Ky(=nky), Echo(=nec), Slice(=nsl), Coil(=nc), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))))]
            dat_cal_mux1 = reshape(permute(fftshift(dat_cal_mux1, 6), [1,2,3,4,6,5]), [nxi, nky, nec, nsl*nz, nc, 1]); % Dim: [Kx(=nxi), Ky(=nky), Echo(=nec), Slice(=nsl*nz), Coil(=nc), SimultaneousSlice Z(=1)], Change the dimensions so that the number of simultaneous slices becomes 1
            
            % Eddy current effects
            eddy = []; % TODO: Not considering eddy current effects for now
            
            % - Pseudo-multiple replica simulation
            [sense_mux1, oneD_grappa_mux1, slice_grappa_mux1, split_slice_grappa_mux1] = pmr_sim_mux_seq(p_mux1, us_msk_mux1, dat_cal_mux1, eddy);
            switch p_mux1.mux_recon_method
                case 'sense'
                    pmr_noise_ref = sense_mux1.noise;
                case '1Dgrappa'
                    pmr_noise_ref = oneD_grappa_mux1.noise;
                case 'slice-grappa'
                    pmr_noise_ref = slice_grappa_mux1.noise;
                case 'split-slice-grappa'
                    pmr_noise_ref = split_slice_grappa_mux1.noise;
            end
    end
end

% Acceleration factor w.r.t. reference scan
switch p.pmr_ref_type
    case 'muxcal'
        R = p.mux;
    case 'mux1'
        R = 1;
end

%% Calculate and output retained SNR maps
pmr_rsnr = zeros(sz(1), sz(2), nec, nsl*nz, 1, nmsk);
for msk_idx = 1 : nmsk
    pmr_rsnr(:, :, :, :, 1, msk_idx) = pmr_noise_ref ./ pmr_noise(:, :, :, :, 1, msk_idx) .* sqrt(R);
end

%% If not in debugging mode, don't output SNR maps and reference noise maps
if ~p.debug
    pmr_snr = [];
    pmr_noise_ref = [];
end

return