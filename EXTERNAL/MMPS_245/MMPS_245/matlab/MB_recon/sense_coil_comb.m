function dat = sense_coil_comb(dat, p)
%
% function dat = sense_coil_comb(dat, p)
%
% Coil-combine single-coil single-slice data using SENSE.
%
% Inputs
%   dat - Single-coil single-slice k-space data. Dim: [Kx(=nx),
%         Ky(=ny_part, if p.partial_ky == true; =p.ny_pres, if p.partial_ky == false),
%         Echo(=nec), Slice(=nsl), Coil(=nc), Time(=p.num_mux_cycle+nt,
%         where nt is number of slice-accelerated time points)].
%   p   - Parameter structure. See mux_epi_params.m for details. Fields
%         used in this function: FE_DIM(=1), PE_DIM(=2), EC_DIM(=3), SL_DIM(=4),
%         C_DIM(=5), T_DIM(=6), KZ_DIM(=7), KEEP_ORIG_SZ(=1), num_mux_cycle,
%         ny_pres, sense_lambda, sense_lambda_default_ratio,
%         add_vcc, ONLY_ACTUAL_COILS, ACTUAL_AND_VIRTUAL_COILS, ONLY_VIRTUAL_COILS.
%
% Output
%   dat - Coil-combined single-slice images. Dim: [X(=nx), Y(=p.ny_pres),
%         Echo(=nec), Slice(=nsl), Coil(=1), Time(=p.num_mux_cycle+nt)].
%
% (c) Kangrong Zhu      Stanofrd University     Dec 2014

datsz = get_dat_sz(dat, p);

% MUX1 undersampling mask
msk.ky = 1 : datsz.y;
msk.kz = ones(1, datsz.y);
msk.omegaz = 0;

% Parameters
p.mux = 1;
p.cap_fov_shift_cal = 1;
p.cap_fov_shift = 1;

% Sensitivity maps
smap = sense_dat_cal_to_smap(dat(:, :, :, :, :, p.num_mux_cycle), p);

% Eddy current effects
eddy = []; % For single-slice data, this should always be empty.

% Combine coil images
dat = mux_recon_sense(dat, msk, p, smap, eddy); % Dim: [X(=nx), Y(=p.ny_pres), Echo(=nec), Slice(=nsl), Coil(=1), Time(=p.num_mux_cycle+nt), SetOfSensitivityMaps].

% Combine images from different sets of sensitivity maps, if multiple sensitivity maps were used in ESPIRiT
if strcmp(p.smap_type, 'espirit') && p.espirit_nmap > 1
    dat = sum(dat, 7);
end

return