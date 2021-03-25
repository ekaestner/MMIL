function [dat, p] = mux_epi_process_data_sense(dat, p)
%
% function [dat, p] = mux_epi_process_data_sense(dat, p)
%
% Process slice-multiplexed EPI data. Use SENSE for reconstruction.
%
% Inputs
%   dat - MUX EPI data, already phase and ramp-sample corrected. The first
%         p.mux_encoded*p.num_mux_cycle time points are mux phase
%         cycling time points(totally p.num_mux_cycle mux phase cycles).
%         The rest time points are accelerated MUX EPI data. Dim: [Kx, Ky,
%         Echo, Slice, Coil, Time].
%   p   - Parameter structure. See mux_epi_params.m for details. The following
%         fields are used in this function: internal_cal, partial_ky, ny_pres,
%         smap_acssz, debug, mux, num_mux_cycle, cap_fov_shift_cal, cap_fov_shift,
%         caipi, mica_br, mica_rand, cap_seed_shift, cap_blip_start, cap_blip_inc, kydir,
%         BOTTOM_UP, PE_DIM, C_DIM, T_DIM, KEEP_ORIG_SZ, inplane_R, return_sos_im,
%         smap, smap_type, espirit_nmap, espirit_ksize, espirit_eigThresh_1,
%         espirit_eigThresh_2, crop_smap, pha_coe(Set in function epi_load_tseries).
%
% Output
%   dat - The reconstructed images.
%         If p.return_sos_im == true, these are the square-root-of-sum-of-squares coil-combined magnitude images.
%         If p.return_sos_im == false, these are the reconstructed complex images.
%         Dim: [X, Y, Echo, Slice*SimultaneousSlice, Coil(Always=1), Time].
%   p   - Output parameter structure. Fields that might have been changed or added:
%         smap_acssz - Might have been changed when calculating sensitivity maps by calling function sense_dat_cal_to_smap.
%         smap       - If ((p.calc_snr_pmr || p.calc_sense_rsnr) && (~isfield(p, 'smap') || isempty(p.smap))), this has been added to save the sensitivity maps for future use.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2013

%% Parameters
sz_dat = get_dat_sz(dat, p);                                % Input data size
nky_input = sz_dat.y;                                       % nky in the input data
if p.partial_ky
    ny_part = nky_input;                                    % Number of acquired ky lines in a partial ky acquisition
end

if isempty(p.ny_pres)
    if p.partial_ky
        error('Must specify ''p.ny_pres'' when partial ky acquisition was used.');
    else
        p.ny_pres = sz_dat.y;
    end
end

if ~isfield(p, 'debug') || isempty(p.debug)
    p.debug = false;
end

use_previous_smap = (isfield(p, 'smap') && ~isempty(p.smap)); % True: Use previously saved sensitivity maps; False: will calculate sensitivity maps from calibration data
if ~use_previous_smap && (p.num_mux_cycle == 0)
    error('No previously saved sensitivity maps and no calibration data to calculate sensitivity maps.');
end

%% Sensitivity maps
% Fully sampled calibration data
if p.num_mux_cycle > 0
    dat_cal = dat(:, :, :, :, :, 1 : p.mux_encoded*p.num_mux_cycle); % Fully sampled slice-multiplexed data, Dim: [Kx(=nx), Ky(=nky_input), Echo(=nec), Slice(=nsl), Coil(=nc), Time(=p.mux_encoded*p.num_mux_cycle)]
    dat_cal = reshape(dat_cal, [sz_dat.x, sz_dat.y, sz_dat.ec, sz_dat.sl, sz_dat.c, p.mux_encoded, p.num_mux_cycle]); % Dim: [Kx(=nx), Ky(=nky_input), Echo(=nec), Slice(=nsl), Coil(=nc), Kz(=p.mux_encoded), Time(=p.num_mux_cycle)]
    dat_cal = mux_dftz(dat_cal, p.T_DIM, p.cap_fov_shift_cal, p.mux, 'decode'); % Dim: [Kx(=nx), Ky(=nky_input), Echo(=nec), Slice(=nsl), Coil(=nc), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Time(=p.num_mux_cycle)]
end

% Sensitivity maps
if use_previous_smap
    if p.debug
        fprintf('  Using previously saved sensitivity maps...\n');
    end
    smap = p.smap; % Load previously saved sensitivity maps
else
    if p.debug
        fprintf('  Calculating sensitivity maps...\n');
    end

    [smap, p] = sense_dat_cal_to_smap(dat_cal(:, :, :, :, :, :, p.num_mux_cycle), p); % Use the last group of mux phase cycling as calibration data. smap Dim: [X(=nx), Y(=p.ny_pres), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Coil(=nc), SetOfSensitivityMaps(=p.espirit_nmap if p.smap_type is 'espirit'; =1 if p.smap_type is 'coil_over_sos')]
    
    if p.calc_snr_pmr || p.calc_sense_rsnr
        p.smap = smap; % Save sensitivity maps for future use
    end
end

%% Coil-combine fully sampled data
if p.num_mux_cycle > 0
    if p.debug
        fprintf('  Coil-combining fully sampled calibration data...\n');
        if p.return_sos_im % Use SOS to combine the coil images of the fully sampled calibration data
            fprintf('   Using SOS...\n');
        else % Use SENSE to combine the coil images of the fully sampled calibration data
            fprintf('   Using SENSE...\n');
        end
    end
    
    dat_cal = fftshift(permute(dat_cal, [1,2,3,4,6,5,7]), 5); % Dim: [Kx(=nx), Ky(=nky_input), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=nc), Time(=p.num_mux_cycle)]
    dat_cal = reshape(dat_cal, [sz_dat.x, nky_input, sz_dat.ec, sz_dat.sl*p.mux, sz_dat.c, p.num_mux_cycle]); % Dim: [Kx(=nx), Ky(=nky_input), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=nc), Time(=p.num_mux_cycle)]
    
    % Use SENSE
    if ~p.return_sos_im
        dat_cal = sense_coil_comb(dat_cal, p); % Output dat_cal Dim: [X(=nx), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=1), Time(=p.num_mux_cycle)]. Input dat_cal Dim: [Kx(=nx), Ky(=nky_input), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=nc), Time(=p.num_mux_cycle)].
    end
    
    % Partial ky
    if p.partial_ky
        if p.return_sos_im
            datin_domain = 'kxky';
            dat_cal = cat(p.PE_DIM, dat_cal, zeros(sz_dat.x, p.ny_pres-nky_input, sz_dat.ec, sz_dat.sl*p.mux, sz_dat.c, p.num_mux_cycle));
        else
            datin_domain = 'xy';
        end
        dat_cal = homodyne(dat_cal, ny_part, p.homodyne_ntran, p.homodyne_niter, datin_domain, 'xy'); % Dim: [X(=nx), Y(=p.ny_pres), Echo(=nec), Slice(=nsl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=1 if ~p.return_sos_im; =nc if p.return_sos_im), Time(=p.num_mux_cycle)]
    else
        if p.return_sos_im
            dat_cal = ifft2c(dat_cal);
        end
    end
    
    % Use SOS
    if p.return_sos_im
        dat_cal = sos(dat_cal, p.C_DIM);
    end
else
    dat_cal = [];
end

%% Recon accelerated data
if p.debug
    fprintf('  Reconstructing accelerated data using SENSE...\n');
end

if size(dat, p.T_DIM) > p.mux_encoded*p.num_mux_cycle
    % ky-kz undersampling mask for the accelerated data
    if p.calc_snr_pmr && isfield(p, 'pmr_us_msk') && ~isempty(p.pmr_us_msk)
        us_msk = p.pmr_us_msk;
    else
        us_msk = get_ky_omegaz_us_msk(p, sz_dat.y, size(dat, p.T_DIM) - p.mux_encoded*p.num_mux_cycle, true);
    end
    
    % Data
    dat = dat(:, us_msk(1).ky, :, :, :, p.mux_encoded*p.num_mux_cycle+1:end); % dat contains accelerated data only. If partial ky was used, dat doesn't contain the unacquired k-space portion. Dim: [Kx, Ky(=nsamp), Echo, Slice, Coil, Time(accelerated time points)]

    % SENSE
    if isfield(p, 'pha_coe')
        pha_coe = p.pha_coe(:, us_msk(1).ky, :, :, :, :);
        if size(pha_coe, 6) == 1 % One set of coefficients for all coils
            pha_coe = repmat(pha_coe, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, sz_dat.c]);
        end
        eddy = get_pha_flt( - pha_coe, sz_dat.x);
    else
        eddy = [];
    end
    [dat, ranks] = mux_recon_sense(dat, us_msk, p, smap, eddy); % dat Dim: [X, Y, Echo, Slice*SimultaneousSlice, Coil(=1), Time, SetOfSensitivityMaps]
    
    % Partial ky
    if p.partial_ky
        if p.debug
            fprintf('  Solving partial ky acquisition...\n');
        end
        if p.homodyne_niter == 0
            fprintf('   Zero-filling...\n');
        else
            fprintf('   Homodyne...\n');
        end

        dat = homodyne(dat, ny_part, p.homodyne_ntran, p.homodyne_niter, 'xy', 'xy'); % Dim: [X, Y, Echo, Slice*SimultaneousSlice, Coil(=1), Time, SetOfSensitivityMaps]
        sz_dat.y = p.ny_pres;
    end
    
    % Combine images for different sets of sensitivity maps, if multiple sensitivity maps were used in ESPIRiT
    if strcmp(p.smap_type, 'espirit') && p.espirit_nmap>1
        dat(isnan(dat)) = 0;
        dat = sum(dat, 7);
    end
else
    dat = [];
end

%% Concatenate fully sampled and accelerated data
dat = cat(p.T_DIM, dat_cal, dat);

%% Zero pad reconstructed images
if p.zpad_image
    sz_dat = get_dat_sz(dat, p);
    dat = ifft2c(zpad(fft2c(dat), [p.zpad_size(1), p.zpad_size(2), sz_dat.ec, sz_dat.sl, sz_dat.c, sz_dat.t, sz_dat.kz]));
end

%% Final images
if p.return_sos_im
    dat = abs(dat);
end

return