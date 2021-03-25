function [smap, p] = sense_dat_cal_to_smap(dat_cal, p)
%
% [smap, p] = sense_dat_cal_to_smap(dat_cal, p)
%
% Calculates sensitivity maps using fully sampled single-slice k-space data.
%
% Inputs
%   dat_cal - Fully sampled single-slice k-space data. 
%             Dim: [Kx(=nx), Ky(=ny_part, if p.partial_ky == true; =p.ny_pres, if p.partial_ky == false),
%             Echo(=nec), Slice(=nsl), Coil(=nc), SimultaneousSlice Z(=nz)].
%   p       - Parameter structure. See mux_epi_params.m for details.
%             Fields used in this function: partial_ky, smap_type, smap_acssz, debug,
%               FE_DIM, PE_DIM, EC_DIM, SL_DIM, C_DIM, T_DIM, KZ_DIM.
%               If p.smap_type is 'espirit': espirit_nmap, espirit_ksize,
%               espirit_eigThresh_1, espirit_eigThresh_2, espirit_tmp_dir;
%               If p.smap_type is 'coil_over_sos': crop_smap.
%
% Outputs
%   smap    - Sensitivity maps. Dim: [X(=nx), Y(=p.ny_pres), Echo(=nec), Slice(=nsl), SimultaneousSlice Z(=nz), Coil(=nc), SetOfSensitivityMaps(=p.espirit_nmap if p.smap_type is 'espirit'; =1 if p.smap_type is 'coil_over_sos')]
%   p       - Parameter structure. With field 'smap_acssz' possibly having been modified.
%
% (c) Kangrong Zhu      Stanford University     June 2013

% Sizes
sz = get_dat_sz(dat_cal, p);
sz.nz = sz.t;
if p.partial_ky
    ny_part = sz.y;
else
    ny_part = [];
end
sz_smap = [sz.x, p.ny_pres];

if isempty(p.smap_acssz.x) || (p.smap_acssz.x > sz.x)
    p.smap_acssz.x = sz.x;
end
if ~(strcmp(p.smap_type, 'espirit') && p.espirit_catim) && (isempty(p.smap_acssz.y) || (p.smap_acssz.y > p.ny_pres))
    p.smap_acssz.y = p.ny_pres;
end
if strcmp(p.smap_type, 'espirit') && p.espirit_catim && (isempty(p.smap_acssz.y) || (p.smap_acssz.y > p.ny_pres*sz.nz))
    p.smap_acssz.y = p.ny_pres * sz.nz;
end
if strcmp(p.smap_type, 'espirit') && ~p.espirit_catim && p.partial_ky && (p.smap_acssz.y > ny_part)
    p.smap_acssz.y = ny_part;
end

% Pad zeros for partial ky
if p.partial_ky
    dat_cal = cat(p.PE_DIM, dat_cal, zeros(sz.x, p.ny_pres-sz.y, sz.ec, sz.sl, sz.c, sz.nz)); % Dim: [Kx, Ky, Echo, Slice, Coil, SimultaneousSlice Z]
    sz.y = p.ny_pres;
end

% Sensitivity maps
if strcmp(p.smap_type, 'espirit') % ESPIRiT sensitivity maps
    if p.debug
        fprintf('   ESPIRiT maps...\n');
    end
    ncalib = [p.smap_acssz.x, p.smap_acssz.y];
    if ~p.espirit_catim % Calculate sensitivity maps slice by slice
        dat_cal = permute(dat_cal, [1,2,5,3,4,6]); % Dim: [Kx, Ky, Coil, Echo, Slice, SimultaneousSlice Z]
        smap = zeros(sz_smap(1), sz_smap(2), sz.c, p.espirit_nmap, sz.ec, sz.sl, sz.nz); % Dim: [X, Y, Coil, SetOfSensitivityMaps, Echo, Slice, SimultaneousSlice Z]
        for im_idx = 1 : sz.ec * sz.sl * sz.nz
            smap(:, :, :, :, im_idx) = sense_map_espirit(dat_cal(:, :, :, im_idx), p.espirit_nmap, ncalib, p.espirit_ksize, p.espirit_eigThresh_1, p.espirit_eigThresh_2, p.partial_ky, ny_part, false, p.espirit_tmp_dir);
        end
        smap = permute(smap, [1,2,5,6,7,3,4]); % Dim: [X, Y, Echo, Slice, SimultaneousSlice Z, Coil, SetOfSensitivityMaps]
    else % Concatenate simultaneous slices into one image
        dat_cal = permute(dat_cal, [1,2,6,5,3,4]); % Dim: [Kx, Ky, SimultaneousSliceZ, Coil, Echo, Slice]
        dat_cal = ifft2c(dat_cal);
        dat_cal = reshape(dat_cal, [sz.x, sz.y*sz.nz, sz.c, sz.ec, sz.sl]); % Dim: [X, Y->SimultaneousSliceZ, Coil, Echo, Slice]
        dat_cal = fft2c(dat_cal); % Dim: [Kx, Ky, Coil, Echo, Slice]
        smap = zeros(sz_smap(1), sz_smap(2)*sz.nz, sz.c, p.espirit_nmap, sz.ec, sz.sl); % Dim: [X, Y->SimultaneousSliceZ, Coil, SetOfSensitivityMaps, Echo, Slice]
        for im_idx = 1 : sz.ec * sz.sl
            smap(:, :, :, :, im_idx) = sense_map_espirit(dat_cal(:, :, :, im_idx), p.espirit_nmap, ncalib, p.espirit_ksize, p.espirit_eigThresh_1, p.espirit_eigThresh_2, false, [], false, p.espirit_tmp_dir);
        end
        smap = reshape(smap, [sz_smap(1), sz_smap(2), sz.nz, sz.c, p.espirit_nmap, sz.ec, sz.sl]);
        smap = permute(smap, [1,2,6,7,3,4,5]); % Dim: [X, Y, Echo, Slice, SimultaneousSliceZ, Coil, SetOfSensitivityMaps]
    end
else % Coil images divided by SOS image
    if p.debug
        fprintf('   Coil over SOS maps...\n');
    end
    dat_cal = permute(dat_cal, [1,2,3,4,6,5]); % Dim: [Kx, Ky, Echo, Slice, SimultaneousSlice Z, Coil]
    acspos = grappa_acspos(struct('x', {sz.x}, 'y', {sz.y}), p.smap_acssz, struct('x',{[1, sz.x]}, 'y',{[1, sz.y]}), p.debug); % Calculate ACS positions using full matrix size
    smap = sense_map_coil_over_sos(dat_cal(acspos.x, acspos.y, :, :, :, :), sz_smap, 6, p.crop_smap); % Dim: [X, Y, Echo, Slice, SimultaneousSlice Z, Coil, SetOfSensitivityMaps(=1)]
end

return