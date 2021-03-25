function dftz_mtx = encode_dftz_mtx(fov_shift, nz)
%
% function dftz_mtx = encode_dftz_mtx(fov_shift, nz)
%
% Calculate Discrete Fourier Transform (DFT) encoding matrix for the slice dimension.
%
% Inputs
%   fov_shift - CAIPI FOV shift. abs(fov_shift) = nkz is the length of the DFTz encoding. Negative for inverse DFTz encoding, positive for forward DFTz encoding.
%   nz        - Number of simultaneous slices.
%
% Output
%   dftz_mtx  - DFTz encoding matrix. Dim: [nkz, nz(z indices ifftshift(-floor(nz/2):1:(ceil(nz/2)-1)))].
%
% (c) Kangrong Zhu      Stanford University     Mar 2014

nkz = abs(fov_shift);
msk.omegaz = encode_dftz_omegaz(fov_shift); % Dim: [1, nkz]
msk.kz = 1 : nkz;
dftz_mtx = encode_ftyz_mtx(msk, 1, nz, 0);  % Dim: [nkz, nz]

return