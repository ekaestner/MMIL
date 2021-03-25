function [reconed, sz_out] = mux_recon_slice_grappa(dat_mux, dat_sing, us_msk, p)
%
% function [reconed, sz_out] = mux_recon_slice_grappa(dat_mux, dat_sing, us_msk, p)
%
% Reconstruct slice-multiplexed data using slice-GRAPPA.
% Reference: Kawin Setsompop, et al. MRM 2012;67(5):1210-24.
%
% Inputs
%   dat_mux  - Accelerated slice-multiplexed k-space data (no mux phase cycling time points).
%              Dim: [Kx(=nx), Ky(=ny_full), Echo(=nec), Slice(=nsl, with slice multiplexing), Coil(=nc), Time(=nt)].
%   dat_sing - Reference single-slice k-space data. Dim: [Kx(=nx), Ky(=ny_full), Echo(=nec), Slice(=nsl),
%              SimultaneousSlice Z(=nz, z indices -floor(nz/2):1:(ceil(nz/2)-1), Coil(=nc)]
%   us_msk   - Undersample mask with fields 'ky', 'kz', 'omegaz'.
%              See get_ky_omegaz_us_msk.m for details.
%   p        - Parameter structure. See mux_epi_params.m for details.
%              Fields used in this function: FE_DIM, PE_DIM, EC_DIM, SL_DIM,
%              C_DIM, T_DIM, KZ_DIM, mux, inplane_R, grappa_domain, mux_kersz_slice_grappa, mux_acssz_slice_grappa, debug.
%
% Outputs
%   reconed  - Reconctructed single-slice k-space data. Dim: [Kx(=nx), Ky(=ny_full), Echo(=nec),
%              Slice(nsl*nz, solved for slice multiplexing, z indices -floor(nz/2):1:(ceil(nz/2)-1), Coil(=nc), Time(=nt)].
%   sz_out   - Size of 'reconed'. A structure array with fields: x, y, ec, sl, c, t, kz.
%
% (c) Kangrong Zhu,     Stanford University     Aug 2014

datsz = get_dat_sz(dat_mux, p);
ny_full = datsz.y;
nz = p.mux;                                        % Number of simultaneous slices

% -- Check inputs
if size(dat_sing, 5) ~= nz
    error('Number of slices in the single-slice reference data doesn''t match number of simultaneous slices.');
end
if nz <= 1                                         % No slice multiplexing to solve
    reconed = dat_mux;
    sz_out  = datsz;
    return;
end

if p.debug
    if p.slice_grappa_odd_even_fit
        fprintf('   Fitting different kernels to odd and even ky lines...\n');
    else
        fprintf('   Fitting one kernel to all ky lines...\n');
    end
end

% -- Only use the acquired data if inplane acceleration was used
if p.inplane_R > 1
    dat_mux = dat_mux(:, us_msk.ky, :, :, :, :);   % Dim: [Kx, Ky(=ny_full/p.inplane_R=length(us_msk.ky)), Echo, Slice, Coil, Time]
    dat_sing = dat_sing(:, us_msk.ky, :, :, :, :); % Dim: [Kx, Ky(=ny_full/p.inplane_R=length(us_msk.ky)), Echo, Slice, SimultaneousSlice Z(=nz), Coil]
end
ny_acq = length(us_msk.ky);

% -- Set up single-slice data matrix
dat_sing = permute(dat_sing, [1,2,6,5,3,4]); % Dim: [Kx, Ky(=ny_acq), Coil, SimultaneousSlice Z(=nz), Echo, Slice]

% -- Interpolation kernels and slice-GRAPPA recon
[acspos, p.mux_acssz_slice_grappa] = grappa_acspos( struct('x',{datsz.x}, 'y',{ny_acq}), p.mux_acssz_slice_grappa, struct('x',{[1, datsz.x]}, 'y',{[1, ny_acq]}), p.debug);
reconed = zeros(datsz.x, ny_full, datsz.ec, datsz.sl, nz, datsz.c, datsz.t); % Dim: [Kx(=nx), Ky(=ny_full), Echo(=nec), Slice(=nsl), simultaneousSlice Z(=nz), Coil(=nc), Time(=nt)]
for echo = 1 : datsz.ec
    for slice = 1 : datsz.sl
        % Interpolation kernels
        dat_mux_for_kernel = reshape(dat_mux(acspos.x, acspos.y, echo, slice, :, 1), [p.mux_acssz_slice_grappa.x, p.mux_acssz_slice_grappa.y, datsz.c]); % Dim: [Kx(=p.mux_acssz_slice_grappa.x), Ky(=p.mux_acssz_slice_grappa.y), Coil(=nc)]
        dat_sing_for_kernel = dat_sing(acspos.x, acspos.y, :, :, echo, slice); % Dim: [Kx(=p.mux_acssz_slice_grappa.x), Ky(=p.mux_acssz_slice_grappa.y), Coil(=nc), SimultaneousSlice Z(=nz)]
        ker = slice_grappa_kernel(dat_mux_for_kernel, dat_sing_for_kernel, us_msk, [p.mux_kersz_slice_grappa.x, p.mux_kersz_slice_grappa.y], p.grappa_domain, [datsz.x, ny_acq], p.use_split_slice_grappa, p.slice_grappa_odd_even_fit);
        
        nker_per_slice = length(ker);
        if p.slice_grappa_odd_even_fit && (nker_per_slice ~= 2)
            error('Number of kernels per slice is not 2 for odd even kernel fitting.');
        end
        if ~p.slice_grappa_odd_even_fit && (nker_per_slice ~= 1)
            error('Number of kernels per slice is not 1 for fitting one kernel to all ky lines.');
        end
        
        % Slice-GRAPPA recon
        if datsz.t > 0
            for ker_idx = 1 : nker_per_slice
                for time = 1 : datsz.t
                    dat_mux_for_recon = reshape(dat_mux(:, :, echo, slice , :, time), [datsz.x, ny_acq, datsz.c]); % Dim: [Kx(=nx), Ky(=ny_acq), Coil(=nc)]
                    tmp = slice_grappa_recon(dat_mux_for_recon, ker{ker_idx}.ker, us_msk, p.grappa_domain); % Dim: [Kx(=nx), Ky(=ny_acq), Coil(=nc), SimultaneousSlice Z(=nz)]
                    if nker_per_slice == 2 % Odd even kernel fitting
                        switch mod(ker{ker_idx}.ky_start, 2)
                            case 0 % Even ky lines
                                ny_indices = 2 : 2 : ny_acq;
                            case 1
                                ny_indices = 1 : 2 : ny_acq;
                        end
                    else % Fit one kernel to all ky lines
                        ny_indices = 1 : 1 : ny_acq;
                    end
                    reconed(:, us_msk.ky(ny_indices), echo, slice, :, :, time) = permute(tmp(:, ny_indices, :, :), [1,2,5,6,4,3,7]);
                end
            end
        end
    end
end

% -- Reshape data
reconed = reshape(reconed, [datsz.x, ny_full, datsz.ec, datsz.sl*nz, datsz.c, datsz.t]); % Dim: [Kx(=nx), Ky(=ny_full), Echo(=nec), Slice(=nsl*nz), Coil(=nc), Time(=nt)]

sz_out = get_dat_sz(reconed, p);

% -- Residual ghost correction for odd-even kernel fitting
if p.slice_grappa_odd_even_fit
    reconed = slice_grappa_odd_even_fit_residual_ghost_correct(reconed, p.pha_coe, p);
end

return