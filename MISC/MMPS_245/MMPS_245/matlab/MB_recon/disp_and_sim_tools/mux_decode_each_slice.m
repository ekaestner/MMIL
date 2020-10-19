function dat = mux_decode_each_slice(dat, nz, us_msk)
% function dat = mux_decode_each_slice(dat, nz, us_msk)
%
% Decode the encoding phase added to each of the simultaneous slices.
%
% Inputs
%   dat    - K-space data. Dim: [Kx, Ky(=ne, number of encodings), Echo, Slice, Coil, Time]
%   nz     - Number of simultaneous slices.
%   us_msk - Undersample mask on the ky-omegaz plane. See get_ky_omegaz_us_msk.m for details.
%
% Output
%   dat    - Decoded data. Dim: [Kx, Ky(=ne), Echo, Slice, Coil, Time,
%            IndividualSlice(e.g. index 1 means the encoding added to the 
%            1st slice has been removed)].
%
% (c) Kangrong Zhu      Stanford University     Oct 2013

PE_DIM = 2;
KEEP_ORIG_SZ = 1;

nsamp = length(us_msk(1).kz);
if nsamp ~= size(dat, PE_DIM)
    error('Number of samples mismatch between undersample mask and data.');
end
nmsk = length(us_msk);
[nx, ny, nec, nsl, nc, nt] = size(dat);

dat = reshape(dat, [nx, ny, nec*nsl*nc, nt]);
res = zeros([nx, ny, nec*nsl*nc, nt, nz]);
for time = 1 : nt
    if nmsk == 1
        us_msk_this_time = us_msk(1);
    else
        us_msk_this_time = us_msk(time);
    end
    
    for decode_slice = 1 : nz
        decode_pha = conj(us_msk_this_time.ftz_pha(:, decode_slice)); % Conjugate of the phase added to each slice for the FTz encoding. Dim: [sample(=nsamp), nz]
        
        for im = 1 : nec*nsl*nc
            res(:, :, im, time, decode_slice) = dat(:, :, im, time) .* repmat(decode_pha(:).', [nx, KEEP_ORIG_SZ]);
        end
    end
end
dat = reshape(res, [nx, ny, nec, nsl, nc, nt, nz]);

return
