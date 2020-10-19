function msk = encode_ftz_pha(msk, nz, neg)
%
% function msk = encode_ftz_pha(msk, nz, neg)
%
% Calculate the encoding phase added to each individual slice for each sample on the ky-omegaz plane.
%
% Inputs
%   msk - Structure for sampling mask on the ky-omegaz plane. Fields 'omegaz' and 'kz' are used in this function. See get_ky_omegaz_us_msk.m for details. 
%   nz  - Number of simultaneous slices.
%   neg - True: use negative value of specified frequencies for slice encoding. False: use specified frequencies.
%
% Output
%   msk - The input 'msk' with the field 'ftz_pha' added.
%         ftz_pha: The FTz encoding phase added to each individual slice for each sample on the ky-omegaz plane. Dim: [Sample(=nsamp), SimultaneousSliceZ(=nz, z indices (-floor(nz/2):1:ceil(nz/2)-1)].
%
% (c) Kangrong Zhu      Stanford University     Feb 2014

if ~exist('neg', 'var') || isempty(neg)
    neg = 0;
end
if neg
    neg = 1;                                   % Ensure a known number is used for power of -1.
else
    neg = 0;
end

nmsk = length(msk);
nsamp = length(msk(1).kz);
z_indices = -floor(nz/2) : 1 : (ceil(nz/2)-1); % Slice indices must center around index 0. e.g. [-1,0,1] for nz=3, [-2,-1,0,1] for nz=4. ([-1,0,1] and [2,0,1] are the same for CAIPI(DFTz encoding) when abs(cap_fov_shift)==nz, but different for MICA(Randomly sampled DTFTz encoding) or for CAIPI when abs(cap_fov_shift)~=nz)
for msk_idx = 1 : nmsk
    msk(msk_idx).ftz_pha = zeros(nsamp, nz);
    for slice = 1 : nz
        z = z_indices(slice);                  % Slice index
        msk(msk_idx).ftz_pha(:, slice) = exp(1i * ((-1)^neg) * msk(msk_idx).omegaz(msk(msk_idx).kz(:)) * z);
    end
end

return