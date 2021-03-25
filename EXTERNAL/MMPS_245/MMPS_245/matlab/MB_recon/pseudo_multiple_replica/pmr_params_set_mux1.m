function p = pmr_params_set_mux1(p)
%
% Pseudo-multiple replica: Set up the parameter structure for MUX1 scan.
%
% Input
%   p - Original parameter structure. See mux_epi_params.m for details.
%
% Output
%   p - Modified parameter structure for a corresponding mux1 scan. Inplane
%       acceleration factor has not been changed.
%
% (c) Kangrong Zhu      Stanford University     Jan 2015

p.mux_excited = 1;
p.mux_encoded = 1;
p.mux = 1;
p.cap_fov_shift = 1;
p.cap_fov_shift_cal = 1;
p.caipi = true; % Use a CAIPI-type acquisition to simulate MUX1 scan.
p.mica_br = false;
p.mica_rand = false;
p.mica_perturbed_caipi = false;
p.mica_poisson = false;
p.cap_blip_start = floor(abs(p.cap_fov_shift)/2);
p.cap_blip_inc = p.inplane_R;
p.cap_fov_shift = (-1)^strcmp(p.dftz_type, 'inverse') * abs(p.cap_fov_shift);

% Clear previously saved parameters that correspond to slice-multiplexed scans. MUX1 data shouldn't use these.
if isfield(p, 'ccmtx')
    p = rmfield(p, 'ccmtx'); % coil compresion matrices
end
if isfield(p, 'smap')
    p = rmfield(p, 'smap'); % coil sensitivity maps
end

return