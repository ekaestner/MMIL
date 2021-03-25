function us_msk = get_caipi_us_msk(p, nky, inc_inp)
%
% Calculate CAIPI undersample mask.
%
% Inputs
%   p       - Parameter structure. See mux_epi_params.m for details.
%             The following fields are used in this function:
%             cap_blip_start, inplane_R, kydir, BOTTOM_UP, cap_fov_shift.
%   nky     - Size of the acquired k-space data (raw k-space data loaded
%             from p-file. Inplane acceleration interleaved with zeros,
%             partial ky not zero-padded) in the PE dimension.
%   inc_inp - True: Include inplane acceleration in the undersample mask.
%             False: Don't include inplane acceleration.
%
% Output
%   us_msk  - CAIPI undersampling mask with fields 'ky', 'kz', 'omegaz'. See get_ky_omegaz_us_msk.m for details.
%
% (c) Kangrong Zhu      Stanford University     July 2014

% Ky indices
us_msk.ky = 1 : 1 : nky;                      % With no inplane acceleration.

% Sampled frequencies, omegaz
[us_msk.omegaz, cal_blips] = encode_dftz_omegaz(p.cap_fov_shift);

% Sampling order of omegaz
fov_shift = abs(p.cap_fov_shift);
acc_blips = mod( (p.cap_blip_start : 1 : p.cap_blip_start+(nky-1)), fov_shift); % Blip indices for the accelerated data when there were no inplane acceleration. Blip indices 0~(abs(p.cap_fov_shift)-1) correspond to -kzmax~kzmax. This is [1,2,0,1,2,0,1,...] for p.cap_blip_start=1(start at kz~=0), fov_shift=3, [2,3,0,1,2,3,0,1...] for p.cap_blip_start=2(start at kz~=0), fov_shift=4.
us_msk.kz = zeros(1, nky);
for idx = p.inplane_R : 1 : nky               % Need to account for the fact that the lower interleaves are by-passed when p.inplane_R > 1
    us_msk.kz(idx) = find(cal_blips == acc_blips(idx-p.inplane_R+1));
end
for idx = p.inplane_R-1 : -1 : 1              % Fill in the lower interleaves
    us_msk.kz(idx) = find(cal_blips == (mod(p.cap_blip_start-(p.inplane_R-idx), fov_shift)));
end

if p.kydir == p.BOTTOM_UP
    us_msk.kz = us_msk.kz(end : -1 : 1);
end

% Inplane acceleration
if inc_inp && (p.inplane_R > 1)
    us_msk = add_inplane_acc(us_msk, p);
end

return

function us_msk = add_inplane_acc(us_msk, p)
% Add inplane acceleration to undersample mask

us_msk.ky = bypass_lower_ileaves(us_msk.ky, p);
us_msk.kz = bypass_lower_ileaves(us_msk.kz, p);

return