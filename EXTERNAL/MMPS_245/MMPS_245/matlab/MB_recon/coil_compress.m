function [dat, p] = coil_compress(dat, p)
%
% function [dat, p] = coil_compress(dat, p)
%
% Conduct coil compression.
% Reference: Zhang T, et al. MRM 2013;69(2):571-82.
%
% Inputs
%   dat  - Input prewhitened data. Dim: [FE(=dsz.x), PE(=dsz.y), Echo(=dsz.ec), Slice(=dsz.sl),
%          Coil(=p.num_coils, original number of coils before compression), Time(=dsz.t)].
%   p    - Parameter structure. See mux_epi_params.m for details. The following
%          fields are used in this function: FE_DIM, PE_DIM, EC_DIM, SL_DIM,
%          C_DIM, T_DIM, debug, cc_method('single' or 'geometric'), gcc_slwin.
%
% Outputs
%   dat  - Coil compressed data. Dim: [FE(=dsz.x), PE(=dsz.y), Echo(=dsz.ec), Slice(=dsz.sl),
%          Coil(=p.num_vcoils, after compression), Time(=dsz.t)].
%   p    - Output parameter structure. The following fields might have been modified or added:
%          num_coils - Set to value of p.num_vcoils after coil compression.
%          ccmtx     - Coil compression matrices. Added if
%                      ((p.calc_snr_pmr || p.calc_sense_rsnr) && (~isfield(p, 'ccmtx') || isempty(p.ccmtx))).
%                      Dim: [Coil(=p.num_coils, before compression), Coil(=p.num_vcoils, after
%                      compression), X(=1 if p.cc_method is 'single'; =dsz.x if p.cc_method is 'geometric'),
%                      Echo(=dsz.ec), Slice(=dsz.sl)]
%          smap      - Sensitivity maps. If input sensitivity maps exist and are maps before coil compression, output maps are after coil compression.
%                      Dim: [X(p.nx_pres), Y(=p.ny_pres), Echo, Slice, SimultaneousSliceZ(=nz, z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))), Coil, SetOfSensitivityMaps(must = 1)]
%
% (c) Kangrong Zhu,     Stanford University     July 2013

if ~isfield(p, 'cc_method') || isempty(p.cc_method)
    p.cc_method = 'geometric';
end
if ~isfield(p, 'cc_dat') || isempty(p.cc_dat)
    p.cc_dat = 'muxcal';
end
input_smap = (isfield(p, 'smap') && ~isempty(p.smap));
if input_smap
    nc_input_smap = size(p.smap, 6);
end
if strcmp(p.cc_dat, 'smap') && ~input_smap
    error('No input sensitivity maps.');
end
if strcmp(p.cc_dat, 'smap') && strcmp(p.smap_type, 'espirit') && (p.espirit_nmap > 1)
    error('Number of ESPIRiT sensitivity maps larger than 1. Can''t be used as calibration for coil compression.');
end

% The fully sampled dimension in the accelerated data
fullsamp_dim = p.FE_DIM;

% Data size before coil compression
dsz = get_dat_sz(dat, p);

% Whether to use previously saved coil compression matrices
use_previous_ccmtx = (isfield(p, 'ccmtx') && ~isempty(p.ccmtx)); % True: Use previously saved coil compression matrices
if p.calc_snr_pmr || p.calc_sense_rsnr
    if ~use_previous_ccmtx % Will calculate coil compression matrices and save them into parameter structure
        switch p.cc_method
            case 'single'
                nx_ccmtx = 1; % All x positions use the same coil compression matrix
            case 'geometric'
                nx_ccmtx = dsz.x; % Each x position uses a unique coil compression matrix
        end
        p.ccmtx = zeros(p.num_coils, p.num_vcoils, nx_ccmtx, dsz.ec, dsz.sl); % Dim: [Coil(=p.num_coils, Before compression), Coil(=p.num_vcoils, After compression), X(=1 if p.cc_method is 'single'; =dsz.x if p.cc_method is 'geometric'), Echo(=dsz.ec), Slice(=dsz.sl)]
    end
end

% Conduct coil compression
dat = permute(dat, [1,2,6,5,3,4]); % Dim: [Kx, Ky, Time, Coil(=original number of coils), Echo, Slice]
if input_smap && (nc_input_smap > p.num_vcoils)
    p.smap = permute(p.smap, [1,2,5,6,3,4]); % Dim: [X(=p.nx_pres), Y(=p.ny_pres), SimultaneousSliceZ, Coil(=original number of coils), Echo, Slice]
end
for echo = 1 : dsz.ec
    for slice = 1 : dsz.sl

        % Data before coil compression for this echo and this slice
        d = dat(:, :, :, :, echo, slice); % Dim: [Kx, Ky, Time, Coil(=original number of coils)]
        
        % Sensitivity maps before coil compression for this echo and this slice
        if input_smap && (nc_input_smap > p.num_vcoils)
            s = p.smap(:, :, :, :, echo, slice); % Dim: [X(=p.nx_pres), Y(=p.ny_pres), SimultaneousSliceZ, Coil(=original number of coils)]
            if strcmp(p.smap_type, 'espirit') && p.espirit_catim
                ssz = size(s);
                s = reshape(s, [ssz(1), ssz(2)*ssz(3), 1, ssz(4)]);
                s = fft2c(s); % Dim: [Kx(=p.nx_pres), Ky(=p.ny_pres*p.mux), Kz(=1), Coil(=original number of coils)]
            else
                s = fft2c(mux_dftz(s, 3, p.cap_fov_shift_cal, p.mux, 'encode')); % Dim: [Kx, Ky, Kz, Coil(=original number of coils)]
            end
        end
        
        % Coil compression matrix
        if use_previous_ccmtx
            if p.debug
                fprintf('  Using previously saved coil compression matrices...\n');
            end
            
            ccmtx = p.ccmtx(:, :, :, echo, slice);
        else
            if p.debug
                fprintf('  Calculating %s coil compression matrices...\n', p.cc_method);
            end
            
            % Calibration data for coil compression
            switch p.cc_dat
                case 'muxcal'
                    cal = d(:, :, p.mux*(p.num_mux_cycle-1)+1 : p.mux*p.num_mux_cycle, :); % Dim: [Kx, Ky, Kz(i.e. Time), Coil(=original number of coils)]
                case 'smap'
                    cal = s; % Dim: [Kx, Ky(=p.ny_pres*p.mux if (strcmp(p.smap_type, 'espirit') && p.espirit_catim); =p.ny_pres otherwise), Kz(=1 if (strcmp(p.smap_type, 'espirit') && p.espirit_catim); =p.mux otherwise), Coil(=original number of coils)]
            end
            
            % Coil compression matrix
            if strcmp(p.cc_method, 'single')
                cal = reshape(cal, [dsz.x*dsz.y*p.mux, dsz.c]); % Dim: [Kx-by-Ky-by-Kz, Coil(=original number of coils)]
                [ccmtx, ed] = eig(cal' * cal); % Returned values are in ascending order of the eigenvalues
                ccmtx = fliplr(ccmtx); % Change into descending order
                ccmtx = ccmtx(:, 1:p.num_vcoils); % Dim: [Nc, p.num_vcoils]
                if p.debug
                    ed = flipud(diag(ed));
                    fprintf('  Echo=%d; Slice=%d; The energy in the retained virtual coils in descending order is\n', echo, slice); disp(round(ed(1:p.num_vcoils))');
                end
            else % p.cc_method == 'geometric'
                ccmtx = calcGCCMtx(cal, fullsamp_dim, p.gcc_slwin); % Not aligned. Dim: [Nc, Nc, Nx]
                ccmtx = alignCCMtx(ccmtx(:, 1:p.num_vcoils, :)); % Aligned. Dim: [Nc, p.num_vcoils, Nx]
            end
            
            % For SNR simulation: Save coil compression matrices into parameter structure
            if p.calc_snr_pmr || p.calc_sense_rsnr
                p.ccmtx(:, :, :, echo, slice) = ccmtx;
            end
        end

        % Coil compression on the whole time series
        if size(ccmtx, 3) == 1 % 'cc_method' is 'single'
            ccmtx = repmat(ccmtx, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, dsz.x, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ]);
        end
        if (size(ccmtx, 5) == 1) && (dsz.sl > 1) % All slices use the same ccmtx
            ccmtx = repmat(ccmtx, [p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, p.KEEP_ORIG_SZ, dsz.sl]);
        end
        d = CC(d, fullsamp_dim, ccmtx); % Treat the time dimension as the Kz dimension when calling function CC.
        dat(:, :, :, 1:p.num_vcoils, echo, slice) = d; % Dim: [Kx, Ky, Time, Coil(=original number of coils), Echo, Slice]. Save compressed data back to matrix 'dat' to save memory.

        % Coil compression on the sensitivity maps
        if input_smap && (nc_input_smap > p.num_vcoils)
            s = CC(s, fullsamp_dim, ccmtx);
            if strcmp(p.smap_type, 'espirit') && p.espirit_catim
                s = ifft2c(s);
                s = reshape(s, [ssz(1), ssz(2), ssz(3), p.num_vcoils]); % Dim: [X, Y, SimultaneousSliceZ, Coil(=number of virtual coils)]
            else
                s = ifft2c(mux_dftz(s, 3, p.cap_fov_shift_cal, p.mux, 'decode')); % Dim: [X, Y, SimultaneousSliceZ, Coil(=number of virtual coils)]
            end
            p.smap(:, :, :, 1:p.num_vcoils, echo, slice) = s; % Dim: [X, Y, SimultaneousSliceZ, Coil(=original number of coils), Echo, Slice]. Save compressed data back to matrix 'p.smap' to save memory.
        end
    end
end

dat = permute(dat(:, :, :, 1:p.num_vcoils, :, :), [1,2,5,6,4,3]); % Coil compressed data, Dim: [FE, PE, Echo, Slice, Coil(=number of virtual coils), Time]
if input_smap && (nc_input_smap > p.num_vcoils)
    p.smap = permute(p.smap(:, :, :, 1:p.num_vcoils, :, :), [1,2,5,6,3,4]); % Coil compressed sensitivity maps, Dim: [X, Y, Echo, Slice, SimultaneousSliceZ, Coil(=number of virtual coils)]
end
p.num_coils = p.num_vcoils;

return