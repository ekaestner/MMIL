function p = mux_epi_params_set_tpoints(p, nt_to_recon)
% function p = mux_epi_params_set_tpoints(p, nt_to_recon)
%
% Set up time-point-related fields in the parameter structure.
%
% Inputs
%   p           - Parameter structure. See mux_epi_params.m for details.
%                 Fields used in this function: num_passes, mux, num_mux_cycle.
%   nt_to_recon - Number of time points to reconstruct, excluding the first
%                 few mux phase cycling time points. Default: All time points.
%
% Output
%   p           - Parameter structure, with the following fields updated according to the input value of nt_to_recon:
%                 nt_to_recon - Number of accelerated time points(excluding the first few slice phase cycling time points) to reconstruct.
%                 tpoints_to_load - Time points to be loaded from the mux epi p-file.
%
% (c) Kangrong Zhu      Stanford University     Jan 2014

if ~exist('nt_to_recon', 'var') || isempty(nt_to_recon)
    nt_to_recon = p.num_passes - p.mux_excited*p.cap_get_ecc - p.mux_encoded*p.num_mux_cycle;
end
if numel(nt_to_recon) ~= 1
    error('''nt_to_recon'' must be a scalar.');
end

p.nt_to_recon = nt_to_recon;

end_time_indx = p.mux_excited*p.cap_get_ecc + p.mux_encoded*p.num_mux_cycle + p.nt_to_recon;
if end_time_indx > p.num_passes
    error('Number of time points to reconstruct exceeds total number of time points in data.');
end
if end_time_indx <= 0
    error('Number of time points to reconstruct requires data collectd at negative time.');
end
p.tpoints_to_load = 1 : end_time_indx;

return