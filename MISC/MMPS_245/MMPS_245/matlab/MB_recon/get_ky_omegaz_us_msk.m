function us_msk = get_ky_omegaz_us_msk(p, nky, nt, inc_inp)
%
% function us_msk = get_ky_omegaz_us_msk(p, nky, [nt=1], [inc_inp=false])
%
% Returns the undersample mask on the ky-omegaz plane for accelerated undersampled MUX EPI acquisitions.
%
% Inputs
%   p       - Parameter structure. See mux_epi_params.m for details.
%             The following fields are used in this function: mux, use_gzblips, caipi,
%             mica_br, mica_rand, cap_seed_shift, cap_blip_start, cap_blip_inc,
%             inplane_R, kydir, BOTTOM_UP, cap_fov_shift.
%   nky     - Size of the acquired k-space data (raw k-space data loaded
%             from p-file. Inplane acceleration interleaved with zeros, 
%             partial ky not zero-padded) in the PE dimension.
%   nt      - Only needed for MICA acquisitions. Number of acquired time points.
%   inc_inp - Only needed for CAIPI-type and direct slice-aliasing acquisitions.
%             True: Include inplane acceleration in the undersample mask.
%             False: Don't include inplane acceleration.
%
% Output
%   us_msk  - A structure describing the undersample pattern on the ky-omegaz 
%             plane (kx dimension is always fully sampled). If the undersample
%             pattern changes with time, us_msk(t) is the undersample mask
%             for the t-th time point. Each us_msk contains fields:
%             ky: The ky index(1,2...nky) each sample corresponds to.
%                 Dim: a vector with length 'nsamp(Number of sampled points
%                 on the ky-omegaz plane)'.
%             kz: The omegaz index(1,2...length(omegaz)) each sample
%                 corresponds to. Dim: a vector with length 'nsamp'.
%             omegaz: The sample candidates in the frequency spectrum of 
%                 the simultaneous slices, in range [-pi, pi]. The order
%                 of omegaz corresponds to the order used in the mux cycling
%                 calibration time points. Dim: [1, abs(p.cap_fov_shift)].
%                 The n-th sample acquires the point (ky(n), omegaz(kz(n)))
%                 on the ky-omegaz plane.
%             ftz_pha: The FTz encoding phase added to each individual slice
%                 for each sample on the ky-omegaz plane. Dim: [Sample(=nsamp),
%                 SimultaneousSliceZ(=nz, z indices (-floor(nz/2):1:(ceil(nz/2)-1))].
% 
% (c) Kangrong Zhu,     Stanford University     Aug 2013

if ~exist('nt', 'var') || isempty(nt)
    nt = 1;
end

if ~exist('inc_inp', 'var') || isempty(inc_inp)
    inc_inp = false;
end

CAP_MAX_PE = 128 * 64;
RAND_MOD = hex2dec('7fffffff');

% Undersample mask
if p.use_gzblips && p.mux > 1                 % Gz blip on && mux>1
    if p.caipi                                % CAIPI-type, one us_msk will be calculated
        us_msk = get_caipi_us_msk(p, nky, inc_inp);
    end
    if p.mica_br || p.mica_rand               % MICA, one or nt us_msk will be calculated
        % Ky indices
        ky_indices = 1 : 1 : nky;
        if p.inplane_R > 1
            ky_indices = bypass_lower_ileaves(ky_indices, p); % Account for inplane acceleration. Raw data matrix is arranged top->bottom(-ky->+ky) in ky no matter whether it was acquired top->bottom or bottom->top.
        end
        
        % Sampled frequency indices
        etl = length(ky_indices);             % Echo train length, i.e. number of samples in the frequency spectrum
        if p.mica_br                          % Bit-reversed MICA. Indices in bit-reversed order
            idx_mono = 0 : 1 : 2^nextpow2(etl)-1; % Indices in monotonically increasing order, padded to length of power of 2
            idx_bitrev = digitrevorder(idx_mono, 2);
            idx_bitrev = idx_bitrev((idx_bitrev>=0) & (idx_bitrev<etl)); % Desired indices are 0,1,...etl-1
        end
        if p.mica_rand
            cap_max_random = min(floor(CAP_MAX_PE/etl), nt);
            mica_pe_ordering = zeros(etl, cap_max_random);
            for pe_ordering_idx = 1 : cap_max_random
                if (pe_ordering_idx == 1) || (p.cap_seed_shift > 0)
                    mica_pe_ordering(:, pe_ordering_idx) = shuffle(etl, pe_ordering_idx-1);
                end
            end
        end
        
        % Number of undersample masks
        if p.cap_seed_shift > 0               % Each time point applies a unique sampling scheme
            nmsk = nt;
        else
            nmsk = 1;                         % All time points used the same sampling scheme
        end
        
        % Undersample masks
        us_msk(nmsk) = struct('ky', {[]}, 'kz', {[]}, 'omegaz', {[]});
        for msk_idx = 1 : nmsk
            us_msk(msk_idx).ky = ky_indices;
            us_msk(msk_idx).omegaz = 2*pi * ((0:1:etl-1) - (etl-1)/2) /etl; % Sampled frequencies, omegaz
            if p.mica_br                      % Bit-reversed MICA
                us_msk(msk_idx).kz = circshift(idx_bitrev+1, [0, -p.cap_seed_shift*(msk_idx-1)]); % Order for sampling omegaz
            else if p.mica_rand               % Random MICA
                    blipstart = mod(p.cap_seed_shift*etl*(msk_idx-1), cap_max_random*etl);
                    us_msk(msk_idx).kz = mica_pe_ordering(mod(blipstart:1:blipstart+etl-1, cap_max_random*etl)+1) + 1;
                end
            end
            if p.kydir == p.BOTTOM_UP
                us_msk(msk_idx).kz = us_msk(msk_idx).kz(end:-1:1);
            end
        end
    end
    if p.mica_perturbed_caipi
        % First get a CAIPI undersmaple mask
        inc_inp = true;
        us_msk_caipi = get_caipi_us_msk(p, nky, inc_inp);
        
        % Ky indices
        us_msk.ky = us_msk_caipi.ky;
        
        % Random integer series
        etl = length(us_msk.ky);
        cap_max_random = min(floor(CAP_MAX_PE/etl), nt);
        mica_pe_ordering = zeros(etl, cap_max_random);
        for pe_ordering_idx = 1 : cap_max_random
            if (pe_ordering_idx == 1) || (p.cap_seed_shift > 0)
                mica_pe_ordering(:, pe_ordering_idx) = randomint(etl, pe_ordering_idx-1, RAND_MOD);
            end
        end
        
        % Random kz perturbation
        blipstart = p.cap_blip_start;
        kz_off_rand = (mica_pe_ordering(mod(blipstart:1:blipstart+etl-1, cap_max_random*etl)+1) - (RAND_MOD-1)/2) / (RAND_MOD-1); % in [-1/2, 1/2]
        kz_off_rand = kz_off_rand * p.cap_kz_rand_pert_range; % in [-1/2*cap_kz_rand_pert_range, 1/2*cap_kz_rand_pert_range]
        
        % Add random kz perturbation to normal CAIPI kz encoding
        blip_polarity = sign(p.cap_fov_shift);
        omegaz_single_blip = - 2*pi/abs(p.cap_fov_shift);
        if p.kydir == p.BOTTOM_UP
            kz_off_rand = kz_off_rand(end: -1 : 1);           % Must correspond to us_msk_caipi.omegaz(us_msk_caipi.kz), which is already flipped if ky is acquired bottom->up
        end
        us_msk.omegaz = us_msk_caipi.omegaz(us_msk_caipi.kz) + blip_polarity * kz_off_rand(:).' * omegaz_single_blip; % Normal CAIPI kz encoding + random kz perturbation
        us_msk.kz = 1 : 1 : etl;
    end
    if p.mica_poisson
        max_ny = 480;
        maxetl = max_ny / p.inplane_R;
        fname = ['poisson_mask_arc' num2str(p.inplane_R) '_capfovshift' num2str(abs(p.cap_fov_shift)) '_maxetl' num2str(maxetl) '.raw'];
        us_msk = poisson_mask_read(fname, nky, p);
    end
else                                          % ~p.use_gzblips || p.mux = 1
    us_msk.ky = 1 : nky;
    us_msk.omegaz = 0;
    us_msk.kz = ones(1, nky);
    if inc_inp && (p.inplane_R > 1)
        us_msk = add_inplane_acc(us_msk, p);
    end
end

% FTz encoding phase added to each individual slice for each sample on the ky-omegaz plane
us_msk = encode_ftz_pha(us_msk, p.mux, 0);

return

function us_msk = add_inplane_acc(us_msk, p)
% Add inplane acceleration to undersample mask

us_msk.ky = bypass_lower_ileaves(us_msk.ky, p);
us_msk.kz = bypass_lower_ileaves(us_msk.kz, p);

return