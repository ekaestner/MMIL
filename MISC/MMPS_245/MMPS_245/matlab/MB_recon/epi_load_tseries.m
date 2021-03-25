function [dat, p, dat_ecc] = epi_load_tseries(pfile, ref, vrgf, refp, slices, slices_pfile, slices_ref, p)
%
% function [dat, p, dat_ecc] = epi_load_tseries(pfile, ref, vrgf, refp, slices, slices_pfile, slices_ref, p)
%
% Load and correct the data for an EPI time series.
%
% Inputs:
%   pfile        - Path to pfile.
%   ref          - Path to 'ref.dat'.
%   vrgf         - Path to 'vrgf.dat'.
%   refp         - Path to reference scan pfile.
%   slices       - Desired slices(Normally ordered indices) to load. Default: all slices.
%   slices_pfile - Slice indices in the p-file corresponding to 'slices'.
%                  'slices_pfile' and 'slices' are different for interleaved
%                  slice acquisition. Default: the same as 'slices'.
%   slices_ref   - Slice indices(Normally ordered indices, not pfile indices)
%                  in the reference scan corresponding to 'slices' in the actual scan.
%                  Only needed if the eddy current correction algorithm will use reference scan pfile.
%   p            - Parameter structure. See mux_epi_params.m for details.
%                  The following fields are used in this function: ychop,
%                  frames, echoes, coils, num_slices, num_coils, nx_pres,
%                  tpoints_to_load, pfile_header(used in function load_raw_tseries).
%
% Outputs:
%   dat          - K-space data. Dim: [FE, PE, Echo, Slice, Coil, Time].
%   p            - The output parameter structure, having the same fields as the input 'p'.
%                  The following fields might be changed inside this function if 
%                  the scan was aborted midway and the p-file was only
%                  partially collected: num_mux_cycle, num_passes, 
%                  nt_to_recon, tpoints_to_load.
%                  The following field will be added if p.md_ecc==true:
%                  pha_coe(EPI x-ky phase correction coefficients, will be
%                  passed to mux_epi_process_data_sense.m. Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Echo,
%                  Slice, SimultaneousSlice Z(indices [ifftshift(-floor(mux/2):1:ceil(mux/2)-1])),
%                  Coil(=1 if p.coil_compress==true; =nc if p.coil_compress==false)]
%                  The following fields will be changed if ((p.mux_encoded
%                  > p.mux_excited) && ~p.recon_cal_use_mux_extra):
%                  cap_fov_shift_cal, mux_encoded.
%   dat_ecc      - if p.cap_get_ecc==true, this is the reference k-space data for the
%                  eddy current correction which is collected by the mux sequence.
%                  Dim: [Kx, Ky, Echo, Slice*SimultaneousSlice, Coil(=nc), Time(=1)]
%                  if p.cap_get_ecc==false, this is [].
%
% (c) Kangrong Zhu,     Stanford University     Aug 2012

% -- Set defaults, check file existence
if ~exist(pfile, 'file')
    error('No p-file.');
end
if ~exist(ref, 'file')
    error('No ref file.');
end
if ~exist(vrgf, 'file') && p.vrgf
    error('No vrgf file.');
end

if ~exist('slices', 'var')
    slices =[];
end
if isempty(slices)
    slices = 1 : p.num_slices;
end
if ~exist('slices_pfile', 'var') || isempty(slices_pfile)
    slices_pfile = slices;
end

% -- Load k-space data
[dat, p] = load_raw_tseries(pfile, p, slices_pfile, p.tpoints_to_load); % Dim: [FE, PE, Echo, Slice, Coil, Time]

% -- Reference data for Nyquist ghost correction
datsz = get_dat_sz(dat, p);
if p.cap_get_ecc
    dat_ecc = dat(:, :, :, :, :, 1:p.mux_excited);
    if p.cap_get_ecc_z                                                  % Need to correct the FTz encoding phase
        us_msk = get_ky_omegaz_us_msk(p, datsz.y, 1, false);            % TODO: Deal with length(us_msk)~=1
        dat_ecc = permute(dat_ecc, [2,6,1,3,4,5]);                      % Dim: [Ky, SimultaneousSliceZ(=mux), Kx, Echo, Slice, Coil]
        for idx = 1 : datsz.x*datsz.ec*datsz.sl*datsz.c
            dat_ecc(:, :, idx) = dat_ecc(:, :, idx) .* conj(us_msk(1).ftz_pha);
        end
        dat_ecc = permute(dat_ecc, [3,1,4,5,2,6]);                      % Dim: [Kx, Ky, Echo, Slice, SimultaneousSliceZ(=mux), Coil]
    else
        dat_ecc = permute(dat_ecc, [1,2,3,4,6,5]);                      % Dim: [Kx, Ky, Echo, Slice, SimultaneousSliceZ(=mux), Coil]
    end
    dat_ecc = reshape(dat_ecc, [datsz.x, datsz.y, datsz.ec, datsz.sl*p.mux_excited, datsz.c]); % Dim: [Kx, Ky, Echo, Slice->SimultaneousSliceZ(num_muxed_slices*mux_excited), Coil]
    
    dat = dat(:, :, :, :, :, p.mux_excited+1:end);
else
    dat_ecc = [];
end

% -- Deal with the case where extra slices are encoded in calibration time points but won't be used in recon
if (p.mux_encoded > p.mux_excited) && ~p.recon_cal_use_mux_extra
    cal_dat_tpoints = 1 : p.mux_encoded*p.num_mux_cycle;
    dat_cal = dat(:, :, :, :, :, cal_dat_tpoints);
    dat = dat(:, :, :, :, :, p.mux_encoded*p.num_mux_cycle+1:end);
    dat_cal = reshape(dat_cal, [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, p.mux_encoded, p.num_mux_cycle]);
    dat_cal = mux_dftz(dat_cal, p.T_DIM, p.cap_fov_shift_cal, p.mux_encoded, 'decode');
    dat_cal = dat_cal(:, :, :, :, :, [1:ceil(p.mux_excited/2), p.mux_encoded+1-floor(p.mux_excited/2):p.mux_encoded], :);
    p.cap_fov_shift_cal = sign(p.cap_fov_shift_cal)*p.mux_excited;
    p.mux_encoded = p.mux_excited;
    dat_cal = mux_dftz(dat_cal, p.T_DIM, p.cap_fov_shift_cal, p.mux_encoded, 'encode');
    dat_cal = reshape(dat_cal, [datsz.x, datsz.y, datsz.ec, datsz.sl, datsz.c, p.mux_encoded*p.num_mux_cycle]);
    dat = cat(p.T_DIM, dat_cal, dat);
end

% -- Save the raw data for pseudo-multiple replica simulation
if p.calc_snr_pmr || p.calc_sense_rsnr
    p.pmr_raw_data = dat(:, :, :, :, :, p.mux_encoded * (p.num_mux_cycle-1) + 1 : p.mux_encoded * p.num_mux_cycle);
end

% -- Process the raw data
[dat, p, dat_ecc] = epi_process_rawdata(dat, p, dat_ecc, ref, vrgf, refp, slices, slices_ref);

return