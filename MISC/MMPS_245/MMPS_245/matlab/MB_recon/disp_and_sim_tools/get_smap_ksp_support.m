function ksup = get_smap_ksp_support(smap, fov, slthick, slspacing)
%
% function ksup = get_smap_ksp_support(smap, fov, slthick, slspacing)
%
% Gets k-space kernel support of the input sensitivity maps in ky and kz.
%
% Inputs
%   smap      - Sensitivity map of a 3D volume. Dim: [X(=nx), Y(=ny), Echo(=1), Slice(=nz), Coil(=nc)].
%   fov       - In-plane field of view for the sensitivity maps, in mm.
%   slthick   - Slice thickness for the sensitivity maps, in mm.
%   slspacing - Slice spacing for the sensitivity maps, in mm.
%
% Output
%   ksup      - A structure for k-space kernel support.
%               Fields 'dky_cycles_per_mm' and 'dkz_cycles_per_mm': k-space kernel support in ky and kz. In cycles/mm.
%               Fields 'dky_nlines' and 'dkz_nlines': number of ky and kz lines the k-space kernel support corresponds to, if the acquisition uses the parameters used for measuring the sensitivity maps.
%
% (c) Kangrong Zhu  Stanford University     June 2014

ONE = 1;
PE_DIM = 2;
SL_DIM = 4;

k_thresh1 = 0.5;                                             % Threshold to exclude x positions with low signals for a particular coil
k_thresh2 = 0.3;                                             % Threshold to determine ky-kz space kernel support

[nx, ny, nec, nz, nc] = size(smap);

% Estimate kernel support in ky and kz for each coil at multiple x positions
k = fftc(fftc(smap, PE_DIM), SL_DIM);                        % Dim: [X(=nx), Ky(=ny), Echo(=1), Kz(=nz), Coil(=nc)]
maxk_in_coils = max(reshape(abs(k), [nx*ny*nz, nc]), [], 1); % Dim: [1, nc]
all_nky_nkz(nc) = struct('dy', {[]}, 'dz', {[]});
for coil = 1 : nc
    maxk_this_coil = maxk_in_coils(coil);
    nsamp = 0;
    for x = 1 : nx
        k_this_coil_this_x = reshape(k(x, :, ONE, :, coil), [ny, nz]);   % Dim: [Ky(=ny), Kz(=nz)]
        if max(abs(k_this_coil_this_x(:))) >= k_thresh1 * maxk_this_coil % Signal larger than k_thresh1 times the maximum signal in this coil
            nsamp = nsamp + 1;
            [all_nky_nkz(coil).nky(nsamp), all_nky_nkz(coil).nkz(nsamp)] = get_ky_kz_support(k_this_coil_this_x, k_thresh2);
        end
    end
end

% Average across x positions for each coil
coil_nky = zeros(nc, 1);
coil_nkz = zeros(nc, 1);
for coil = 1 : nc
    coil_nky(coil) = mean(all_nky_nkz(coil).nky);
    coil_nkz(coil) = mean(all_nky_nkz(coil).nkz);
end

% Average across coils
nky = mean(coil_nky);
nkz = mean(coil_nkz);

% Output
delta_ky = 1 / fov;                      % Ky resolution, in cycles/mm
delta_kz = 1 / (nz*(slthick+slspacing)); % Kz resolution, in cycles/mm
ksup.dky_cycles_per_mm = nky * delta_ky;
ksup.dkz_cycles_per_mm = nkz * delta_kz;
ksup.dky_nlines = nky;
ksup.dkz_nlines = nkz;

return

function [dky, dkz] = get_ky_kz_support(dat_kykz, thresh)
%
% Finds the largest span in ky and kz for the signal region whose amplitude is equal to or larger than thresh*maxKyKzSignalAmplitude.
%
% Inputs
%   dat_kykz - Data in ky-kz space. Dim: [Ky(=nky), Kz(=nkz)]
%   thresh   - Threshold on the magnitude.
%
% Outputs
%   dky, dkz - Largest span in ky and kz for the signal region whose amplitude is equal to or larger than thresh*maxKyKzSignalAmplitude.
%

nky = size(dat_kykz, 1);
maxk = max(abs(dat_kykz(:)));
idx = find(abs(dat_kykz) >= thresh*maxk);
ky_idx = mod(idx-1, nky) + 1;
kz_idx = ceil(idx/nky);
dky = max(ky_idx) - min(ky_idx);
dkz = max(kz_idx) - min(kz_idx);

return