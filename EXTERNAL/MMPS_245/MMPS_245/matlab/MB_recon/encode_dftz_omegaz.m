function [omegaz, cal_blips] = encode_dftz_omegaz(fov_shift)
%
% function [omegaz, cal_blips] = encode_dftz_omegaz(fov_shift)
%
% Calculate Discrete Fourier Transform (DFT) encoding frequencies for the slice dimension.
%
% Input
%   fov_shift - CAIPI FOV shift. abs(fov_shift) is the length of the DFTz encoding. Negative for inverse DFTz encoding, positive for forward DFTz encoding.
%
% Outputs
%   omegaz    - DFTz encoding frequencies, in range [-pi, pi]. The order of omegaz corresponds to the order used in the mux cycling calibration time points. Dim: [1, abs(fov_shift)].
%   cal_blips - Blip indices for the calibration time points if the calibration were acquired using a CAIPI FOV shift of 'fov_shift'. Blip indices 0~(abs(fov_shift)-1) correspond to -kzmax~kzmax. cal_blips is [1,2,0] for abs(fov_shift)=3, [2,3,0,1] for abs(fov_shift)=4.
%
% (c) Kangrong Zhu  Stanford University     Feb 2014

blip_polarity = sign(fov_shift);                        % -1 for inverse DFTz encoding, +1 for forward DFTz encoding.
fov_shift = abs(fov_shift);
cap_blip_start_cal = floor(fov_shift/2);                % Staring blip index of the calibration data when the calibration were acquired using a CAIPI FOV shift of 'fov_shift'. Indices 0~(abs(fov_shift)-1) correspond to -kzmax~kzmax.
cal_blips = mod( (cap_blip_start_cal : 1 : cap_blip_start_cal+fov_shift-1), fov_shift);
cal_kz = blip_polarity * (cal_blips - (fov_shift-1)/2); % kz values. When blip_polarity = -1, this is [0, -1, 1] for abs(fov_shift) = 3, and [-0.5, -1.5, 1.5, 0.5] for abs(fov_shift) = 4.
omegaz = - 2*pi * cal_kz / fov_shift;

return