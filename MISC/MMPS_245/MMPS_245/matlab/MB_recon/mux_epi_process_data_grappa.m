function dat = mux_epi_process_data_grappa(dat, p)
%
% function dat = mux_epi_process_data_grappa(dat, p)
%
% Process slice-multiplexed EPI data. Use GRAPPA for reconstruction.
%
% Inputs
%   dat - MUX EPI data, already phase and ramp-sample corrected. The first
%         p.mux_encoded*p.num_mux_cycle time points are mux phase
%         cycling time points(totally p.num_mux_cycle mux phase cycles).
%         The rest time points are accelerated MUX EPI data. Dim: [Kx, Ky,
%         Echo, Slice, Coil, Time].
%   p   - Parameter structure. See mux_epi_params.m for details. The following
%         fields are used in this function: inplane_R, mux, partial_ky,
%         ny_pres, num_mux_cycle, debug, inplane_kersz, inplane_acssz,
%         mux_kersz_1d_grappa, mux_acssz_1d_grappa, PE_DIM, C_DIM, T_DIM, caipi, cal_dat_tpoints,
%         cap_fov_shift_cal, cap_fov_shift, grappa_domain, return_sos_im.
%
% Output
%   dat - Reconstructed images. Dim: [FE, PE, Echo, Slice*SimultaneousSlice, Coil(=1, if
%         (p.return_sos_im == true) || (p.return_sos_im == false && p.sense_coil_comb == true);
%         =nc, otherwise), Time].
%
% Calibration data for GRAPPA
%   (A) Inplane
%         time points in dat to recon                                         -- time points in dat used as calibration
%         1 : abs(cap_fov_shift_cal)*num_mux_cycling(Not accelerated inplane) -- NONE
%         abs(cap_fov_shift_cal)*num_mux_cycling+1 : end                      -- abs(cap_fov_shift_cal)*(num_mux_cycling-1)+1 (Zblip Off) or synthesized from abs(cap_fov_shift_cal)*(num_mux_cycling-1)+1 : abs(cap_fov_shift_cal)*num_mux_cycling (Zblip On)
%   (B) Through-plane
%         time points in dat to recon                                         -- time points in dat used as calibration
%         1 : abs(cap_fov_shift_cal)*num_mux_cycling                          -- NONE
%         abs(cap_fov_shift_cal)*num_mux_cycling+1 : end                      -- abs(cap_fov_shift_cal)*(num_mux_cycling-1)+1 : abs(cap_fov_shift_cal)*num_mux_cycling
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

%% -- Parameters
FIRST = 1;                                                               % Used to choose the first time point in mux phase cycling

sz_dat = get_dat_sz(dat, p);

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

% -- Undersample mask on the ky-omegaz plane, without inplane acceleration
if p.use_gzblips && p.mux>1 && p.inplane_R>1
    us_msk = get_ky_omegaz_us_msk(p, sz_dat.y, [], false);
end

% -- For inplane accelerated scans
if p.inplane_R > 1
    % - Calibration data for solving the inplane aliasing
    if p.use_gzblips && p.mux>1                                          % Gz blip On acquisition
        dat_cal_raw = dat(:, :, :, :, :, p.cal_dat_tpoints);             % Last group of mux phase cycling
        if p.cap_fov_shift_cal ~= p.cap_fov_shift                        % CAIPI FOV shift different for calibration vs. accelerated data
            dat_cal_raw = mux_dftz(dat_cal_raw, p.T_DIM, p.cap_fov_shift_cal, p.mux, 'decode'); % Solve individual slices. Dim: [Kx(=sz_dat.x), Ky(=sz_dat.y), Echo(=sz_dat.ec), Slice(=sz_dat.sl), Coil(=sz_dat.c), SimultaneousSlice Z(=p.mux)]
            dat_cal_raw = mux_dftz(dat_cal_raw, p.T_DIM, p.cap_fov_shift, p.mux, 'encode');     % Synthesize fully sampled calibration data with a CAIPI FOV shift of p.cap_fov_shift. Dim: [Kx(=sz_dat.x), Ky(=sz_dat.y), Echo(=sz_dat.ec), Slice(=sz_dat.sl), Coil(=sz_dat.c), Time(=abs(p.cap_fov_shift)]
        end
        dat_cal = zeros(sz_dat.x, sz_dat.y, sz_dat.ec, sz_dat.sl, sz_dat.c);
        for idx = 1 : sz_dat.y
            dat_cal(:, idx, :, :, :) = dat_cal_raw(:, us_msk.ky(idx), :, :, :, us_msk.kz(idx)); % Synthesize a 2D slice-accelerated image that is fully sampled in-plane.
        end
    else                                                                 % Direct aliasing of simultaneous slices
        dat_cal = dat(:, :, :, :, :, p.mux*(p.num_mux_cycle-1) + FIRST); % Use the 1st time point in the last group of mux phase cycling.
    end
    sz_cal = get_dat_sz(dat_cal, p);
    
    % - Default parameters for inplane GRAPPA recon
    if isempty(p.inplane_kersz.x)
        p.inplane_kersz.x = 7;
    end
    if isempty(p.inplane_kersz.y)
        p.inplane_kersz.y = 4;
    end
    if isempty(p.inplane_acssz.x)
        p.inplane_acssz.x = sz_cal.x;
    end
    if isempty(p.inplane_acssz.y)
        if p.partial_ky
            p.inplane_acssz.y = p.ny_pres;                               % Always use all available calibration data.
        else
            p.inplane_acssz.y = sz_cal.y;
        end
    end
end

% -- For through-plane accelerated scans
if p.mux >1
    % - Default parameters for through-plane GRAPPA recon
    if isempty(p.mux_kersz_1d_grappa.x)
        p.mux_kersz_1d_grappa.x = 7;
    end
    if isempty(p.mux_kersz_1d_grappa.y)
        p.mux_kersz_1d_grappa.y = 4;
    end
    if isempty(p.mux_acssz_1d_grappa.x)
        p.mux_acssz_1d_grappa.x = 32;
    end
    if isempty(p.mux_acssz_1d_grappa.y)
        p.mux_acssz_1d_grappa.y = p.ny_pres * (p.mux-1);
    end
end

%% -- Calculate inplane GRAPPA interpolation kernel.
if p.inplane_R > 1
    if p.debug                                                           % For debugging
        im_cal = sos(ifft2c(dat_cal), p.C_DIM);
        fprintf('  Calculating inplane GRAPPA kernel...\n');
    end

    % - Position for ACS signal
    if p.partial_ky
        total_cal_ny = p.ny_pres;                                        % Calculate ACS positions using full matrix size.
    else
        total_cal_ny = sz_cal.y;
    end
    [inplane_acspos, p.inplane_acssz] = grappa_acspos( struct('x',{sz_cal.x}, 'y',{total_cal_ny}), ...
        p.inplane_acssz, struct('x',{[1, sz_cal.x]}, 'y',{[1, sz_cal.y]}), p.debug);

    % - Inplane kernel
    switch p.grappa_domain
        case 'kSpace'
            inplane_ker = zeros(p.inplane_kersz.x, p.inplane_R*ceil(p.inplane_kersz.y/2)*2-1, ...
                sz_cal.c, sz_cal.c, sz_cal.ec, sz_cal.sl, sz_cal.t);
        case 'imSpace'
            inplane_ker = zeros(sz_dat.x, sz_dat.y, sz_cal.c, sz_cal.c, sz_cal.ec, sz_cal.sl, sz_cal.t);
    end

    for time = 1 : sz_cal.t
        for slice = 1 : sz_cal.sl
            for echo = 1 : sz_cal.ec

                acs_dat = reshape( dat_cal(inplane_acspos.x, inplane_acspos.y, echo, slice, :, time), ...
                    [p.inplane_acssz.x, p.inplane_acssz.y, sz_cal.c]);   % Dim: [FE, PE, Coil]

                inplane_ker(:,:,:,:, echo, slice, time) = ...
                    grappa_kernel(acs_dat, [p.inplane_kersz.x, p.inplane_kersz.y], ...
                    p.inplane_R, p.grappa_domain, [sz_dat.x, sz_dat.y]);
            end
        end
    end

end

%% -- Process accelerated time points.

if p.debug
    im_raw = sos(ifft2c(dat), p.C_DIM);
end

% -- Inplane GRAPPA recon
if p.inplane_R > 1
    if p.debug
        fprintf('  Solving inplane aliasing...\n');
    end

    for time = 1 : sz_dat.t
        if time <= p.mux_encoded*p.num_mux_cycle                         % Skip inplane interpolation for the first 'p.mux_encoded*num_mux_cycle' time points.
            continue;
        end

        for slice = 1 : sz_dat.sl
            for echo = 1 : sz_dat.ec

                us_dat = reshape( dat(:, :, echo, slice, :, time), ...
                    [sz_dat.x, sz_dat.y, sz_dat.c]);                     % Undersampled data, dim: [FE, PE, Coil]
                
                ker = inplane_ker(:, :, :, :, echo, slice, FIRST);       % Currently sz_cal.t is always 1, so always use the FIRST time point.
                
                dat(:, :, echo, slice, :, time) = reshape( ...
                    grappa_recon(us_dat, ker, p.grappa_domain), [sz_dat.x, sz_dat.y, 1, 1, sz_dat.c]);
                
                if strcmp(p.grappa_domain, 'imSpace')
                    dat(:, :, echo, slice, :, time) = fft2c(dat(:, :, echo, slice, :, time));
                end
                
            end
        end
    end
    
    if p.debug
        im_solved_inplane = sos(ifft2c(dat), p.C_DIM);
    end
end

% -- If partial ky acquisition was used, zero-pad k-space to the prescribed size.
if p.partial_ky
    dat = cat(p.PE_DIM, dat, zeros(sz_dat.x, ...
        p.ny_pres-sz_dat.y, sz_dat.ec, sz_dat.sl, sz_dat.c, sz_dat.t));

    ny_part = sz_dat.y;                                                  % # of acquired ky lines in partial ky acquisition
    sz_dat.y = p.ny_pres;
end

% -- Through-plane GRAPPA recon.
if p.mux > 1
    if p.debug
        fprintf('  Solving through plane aliasing...\n');
    end

    [dat, sz_dat] = mux_recon_1Dgrappa(dat, p);                          % Dim:  [FE, PE, Echo, Slice, Coil, TemporalPhase]
    
    if p.debug
        if strcmp(p.grappa_domain, 'imSpace')
            im_solved_thrplane = sos(dat, p.C_DIM);
        else                                                             % if p.grappa_domain is 'kSpace'
            im_solved_thrplane = sos(ifft2c(dat), p.C_DIM);
        end
    end
end

% -- Change to image space. Combine coil images using SENSE if needed.
if ~p.return_sos_im && p.sense_coil_comb                                 % Combine coil images using SENSE
    if p.debug
        fprintf('  SENSE coil combination...\n');
    end
    if strcmp(p.grappa_domain, 'imSpace') && (p.mux > 1)                 % Change to k-space
        dat = fft2c(dat);
    end
    if p.partial_ky
        ny_to_use = ny_part;
    else
        ny_to_use = sz_dat.y;
    end
    dat = sense_coil_comb(dat(:, 1 : ny_to_use, :, :, :, :), p);         % Dim: [X, Y(=p.ny_pres), Echo, Slice, Coil, Time]
    sz_dat.c = 1;
else                                                                     % Don't combine coil images
    if strcmp(p.grappa_domain, 'kSpace') || (p.mux == 1)                 % Change to image space if p.mux=1 or if through-plane GRAPPA returned kspace data
        dat = ifft2c(dat);
    end
end

% -- If partial k-space acquisition was used, recover missing k-space data.
if p.partial_ky
    if p.debug
        fprintf('  Solving partial ky acquisition...\n');
        if p.homodyne_niter == 0
            fprintf('   Zero-filling...\n');
        else
            fprintf('   Homodyne...\n');
        end
    end
    % PWW add loop to reduce memory usage
    % dat = homodyne(dat, ny_part, p.homodyne_ntran, p.homodyne_niter, 'xy', 'xy'); % Data input to and output from homodyne are all image space data
    for t_idx = 1:size(dat,6)
        dat(:,:,:,:,:,t_idx) = homodyne(dat(:,:,:,:,:,t_idx), ny_part, p.homodyne_ntran, p.homodyne_niter, 'xy', 'xy'); % Data input to and output from homodyne are all image space data
    end
end

% -- Zero pad reconstructed images
if p.zpad_image
    dat = ifft2c(zpad(fft2c(dat), [p.zpad_size(1), p.zpad_size(2), sz_dat.ec, sz_dat.sl, sz_dat.c, sz_dat.t, sz_dat.kz]));
end

% -- Final images
if p.return_sos_im
    dat = sos(dat, p.C_DIM);
end

return
