function dat = mux_epi_process_data_slice_grappa(dat, p)
%
% function dat = mux_epi_process_data_slice_grappa(dat, p)
%
% Process slice-multiplexed EPI data. Use slice-GRAPPA for reconstruction.
% Reference: Kawin Setsompop, et al. MRM 2012;67(5):1210-24.
%
% Inputs
%   dat - MUX EPI data, already phase and ramp-sample corrected. The first
%         p.mux_encoded*p.num_mux_cycle time points are mux phase
%         cycling time points(totally p.num_mux_cycle mux phase cycles).
%         The rest time points are accelerated MUX EPI data. Dim: [Kx, Ky,
%         Echo, Slice, Coil, Time].
%   p   - Parameter structure. See mux_epi_params.m for details. The following
%         fields are used in this function: FE_DIM, PE_DIM, EC_DIM, SL_DIM,
%         C_DIM, T_DIM, KZ_DIM, mux, ny_pres, partial_ky, debug, mux_kersz_slice_grappa,
%         mux_acssz_slice_grappa, inplane_R, inplane_kersz, inplane_acssz, cap_fov_shift_cal,
%         num_mux_cycle, grappa_domain, homodyne_niter, homodyne_ntran,
%         return_sos_im, use_gzblips, caipi, mica_br, mica_rand, cap_seed_shift,
%         cap_blip_start, cap_blip_inc, kydir, BOTTOM_UP, cap_fov_shift.
%
% Output
%   dat - Reconstructed images. Dim: [X, Y, Echo, Slice*SimultaneousSlice, Coil(=1, if
%         (p.return_sos_im == true) || (p.return_sos_im == false && p.sense_coil_comb == true);
%         =nc, otherwise), Time].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2014

%% Parameters
datsz = get_dat_sz(dat, p);
nz = p.mux;                                                              % Number of simultaneous slices

if isempty(p.ny_pres)
    if p.partial_ky
        error('Must specify ''p.ny_pres'' when partial ky acquisition was used.');
    else
        p.ny_pres = datsz.y;
    end
end

if ~isfield(p, 'debug') || isempty(p.debug)
    p.debug = false;
end

% Default parameters for through-plane slice-GRAPPA recon
if nz > 1
    if isempty(p.mux_kersz_slice_grappa.x)
        p.mux_kersz_slice_grappa.x = 7;
    end
    if isempty(p.mux_kersz_slice_grappa.y)
        p.mux_kersz_slice_grappa.y = 7;
    end
    if isempty(p.mux_acssz_slice_grappa.x)
        p.mux_acssz_slice_grappa.x = datsz.x;
    end
    if isempty(p.mux_acssz_slice_grappa.y)
        p.mux_acssz_slice_grappa.y = p.ny_pres/p.inplane_R; % Use all available k-space data. For partial ky acquisitions, including the zeros in the kernel fitting process doesn't affect the results.
    end
end

% Default parameters for in-plane GRAPPA recon
if isempty(p.inplane_kersz.x)
    p.inplane_kersz.x = 7;
end
if isempty(p.inplane_kersz.y)
    p.inplane_kersz.y = 4;
end
if isempty(p.inplane_acssz.x)
    p.inplane_acssz.x = datsz.x;
end
if isempty(p.inplane_acssz.y)
    p.inplane_acssz.y = p.ny_pres;                                       % Always use all available calibration data.
end

%% Separate calibration and accelerated data
% If partial ky acquisition was used, zero-pad k-space to the prescribed size.
if p.partial_ky
    dat = cat(p.PE_DIM, dat, zeros(datsz.x, p.ny_pres-datsz.y, datsz.ec, datsz.sl, datsz.c, datsz.t));

    ny_part = datsz.y;                                                   % # of acquired ky lines in partial ky acquisition
    datsz.y = p.ny_pres;
end

dat_ref = dat(:, :, :, :, :, 1:abs(p.cap_fov_shift_cal)*p.num_mux_cycle);
dat_ref = reshape(dat_ref, [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, abs(p.cap_fov_shift_cal), p.num_mux_cycle]);
dat_ref = mux_dftz(dat_ref, p.T_DIM, p.cap_fov_shift_cal, nz, 'decode'); % Solve individual slices. Dim: [Kx(=datsz.x), Ky(=datsz.y), Echo(=datsz.ec), Slice(=datsz.sl), Coil(=datsz.c), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Time(=p.num_mux_cycle)]
dat_ref = permute(dat_ref, [1,2,3,4,6,5,7]);                             % Dim: [Kx(=datsz.x), Ky(=datsz.y), Echo(=datsz.ec), Slice(=datsz.sl), SimultaneousSlice Z(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Coil(=datsz.c), Time(=p.num_mux_cycle)]
dat_ref = fftshift(dat_ref, 5);                                          % Dim: [Kx(=datsz.x), Ky(=datsz.y), Echo(=datsz.ec), Slice(=datsz.sl), SimultaneousSlice Z(=nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=datsz.c), Time(=p.num_mux_cycle)]

dat_cal = dat_ref(:, :, :, :, :, :, p.num_mux_cycle);                    % Last group of mux phase cycling, Dim: [Kx(=datsz.x), Ky(=datsz.y), Echo(=datsz.ec), Slice(=datsz.sl), SimultaneousSlice Z(=nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=datsz.c)]
dat = dat(:, :, :, :, :, abs(p.cap_fov_shift_cal)*p.num_mux_cycle+1:end);

% Now dat_ref, dat_cal, and dat are in k-space

%% Through-plane slice-GRAPPA recon
if p.mux > 1
    if p.debug
        if p.use_split_slice_grappa
            str = ' split-';
        else
            str = ' ';
        end
        fprintf(['  Conducting through-plane' str 'slice-GRAPPA recon...\n']);
    end
    
    % ky-kz undersampling mask for the accelerated data
    if p.calc_snr_pmr && isfield(p, 'pmr_us_msk') && ~isempty(p.pmr_us_msk) % Enable the option to pass the undersampling mask via the parameter structure, so that this routine can be used for simulation purpose.
        us_msk = p.pmr_us_msk;
    else
        us_msk = get_ky_omegaz_us_msk(p, datsz.y, size(dat, p.T_DIM), true); % Undersampling mask on the ky-omegaz plane, with inplane acceleration
    end
    nmsk = length(us_msk);
    if nmsk > 1 % Undersampling pattern changes for each time point
        dat_reconed = zeros(datsz.x, datsz.y, datsz.ec, datsz.sl*p.mux, datsz.c, size(dat, p.T_DIM)); % Dim: [Kx(=datsz.x), Ky(=datsz.y), Echo(=datsz.ec), Slice(=datsz.sl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=datsz.c), Time(=Number of accelerated mux acquisiitons)]
        for msk_idx = 1 : nmsk
            dat_reconed(:, :, :, :, :, msk_idx) = mux_recon_slice_grappa(dat(:, :, :, :, :, msk_idx), dat_cal, us_msk(msk_idx), p);
        end
        dat = dat_reconed;
        clear dat_reconed;
    else
        dat = mux_recon_slice_grappa(dat, dat_cal, us_msk(1), p); % Dim: [Kx(=datsz.x), Ky(=datsz.y), Echo(=datsz.ec), Slice(=datsz.sl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1)), Coil(=datsz.c), Time(=Number of accelerated mux acquisiitons)]
    end
end

dat_cal = reshape(dat_cal, [datsz.x, datsz.y, datsz.ec, datsz.sl*nz, datsz.c]); % Dim: [Kx(=datsz.x), Ky(=datsz.y), Echo(=datsz.ec), Slice(=datsz.sl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))), Coil(=datsz.c)]

% Now dat_ref, dat_cal and dat are in k-space.

%% In-plane GRAPPA recon
if p.inplane_R > 1
    if p.debug
        fprintf('  Conducting inplane GRAPPA recon...\n');
    end

    % Position for ACS signal
    [inplane_acspos, p.inplane_acssz] = grappa_acspos( struct('x',{datsz.x}, 'y',{p.ny_pres}), ...
        p.inplane_acssz, struct('x',{[1, datsz.x]}, 'y',{[1, datsz.y]}), p.debug);
    
    % In-plane GRAPPA recon
    for echo = 1 : datsz.ec
        for slice = 1 : datsz.sl*nz
            acs_dat = reshape( dat_cal(inplane_acspos.x, inplane_acspos.y, echo, slice, :), ...
                [p.inplane_acssz.x, p.inplane_acssz.y, datsz.c]);        % Dim: [Kx, Ky, Coil]
            ker = grappa_kernel(acs_dat, [p.inplane_kersz.x, p.inplane_kersz.y], ...
                p.inplane_R, p.grappa_domain, [datsz.x, datsz.y]);       % Dim: [FE, PE, Coil(=datsz.c, source), Coil(=datsz.c, target)]
            
            for time = 1 : size(dat, p.T_DIM)
                us_dat = reshape( dat(:, :, echo, slice, :, time), [datsz.x, datsz.y, datsz.c]); % Undersampled data, dim: [Kx, Ky, Coil]
                dat(:, :, echo, slice, :, time) = reshape(grappa_recon(us_dat, ker, p.grappa_domain), [datsz.x, datsz.y, 1, 1, datsz.c]);
            end
        end
    end
    
    if ~p.return_sos_im && p.sense_coil_comb % In this case, transform dat into k-space
        if strcmp(p.grappa_domain, 'imSpace')
            dat = fft2c(dat);
        end
    else % In this case, transform dat into image space
        if strcmp(p.grappa_domain, 'kSpace')
            dat = ifft2c(dat);
        end
    end
else
    if ~p.return_sos_im && p.sense_coil_comb
        % Do nothing
    else
        dat = ifft2c(dat);
    end
end

% Now dat is in k-space if (~p.return_sos_im && p.sense_coil_comb), is in image space otherwise.

% Concatenate fully sampled reference data and accelerated data
dat_ref = reshape(dat_ref, [datsz.x, datsz.y, datsz.ec, datsz.sl*nz, datsz.c, p.num_mux_cycle]); % Dim: [Kx(=datsz.x), Ky(=datsz.y), Echo(=datsz.ec), Slice(=datsz.sl*nz, z indices -floor(nz/2):1:(ceil(nz/2)-1))), Coil(=datsz.c), Time(=p.num_mux_cycle)]
if p.slice_grappa_odd_even_fit
    dat_ref = slice_grappa_odd_even_fit_residual_ghost_correct(dat_ref, p.pha_coe, p);
end
if ~p.return_sos_im && p.sense_coil_comb
    % Do nothing
else
    dat_ref = ifft2c(dat_ref);
end
dat = cat(p.T_DIM, dat_ref, dat);

% Now dat is in k-space if (~p.return_sos_im && p.sense_coil_comb), is in image space otherwise.

%% Combine coil images using SENSE, if needed
if ~p.return_sos_im && p.sense_coil_comb
    if p.partial_ky
        ny_to_use = ny_part;
    else
        ny_to_use = datsz.y;
    end
    dat = sense_coil_comb(dat(:, 1 : ny_to_use, :, :, :, :), p); % Dim: [X, Y(=p.ny_pres), Echo, Slice, Coil, Time]
    datsz.c = 1;
end

% Now dat is in image space.

%% Homodyne
if p.partial_ky
    if p.debug
        fprintf('  Solving partial ky acquisition...\n');
    end
    if p.homodyne_niter == 0
        fprintf('   Zero-filling...\n');
    else
        fprintf('   Homodyne...\n');
    end
    dat = homodyne(dat, ny_part, p.homodyne_ntran, p.homodyne_niter, 'xy', 'xy'); % Data input to and output from homodyne are all image space data
end

%% Zero pad reconstructed images
if p.zpad_image
    datsz = get_dat_sz(dat, p);
    dat = ifft2c(zpad(fft2c(dat), [p.zpad_size(1), p.zpad_size(2), datsz.ec, datsz.sl, datsz.c, datsz.t, datsz.kz]));
end

%% Final images
if p.return_sos_im
    dat = sos(dat, p.C_DIM);
end

return