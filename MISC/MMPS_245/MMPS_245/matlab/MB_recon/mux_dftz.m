function dat = mux_dftz(dat, dim, fov_shift, nz, enc_or_dec)
%
% function dat = mux_dftz(dat, [dim=ndims(dat)], [fov_shift=-size(dat,dim)], [nz=size(dat,dim)], [enc_or_dec='decode'])
%
% Encode or decode the DFTz encoding on simultaneous slices.
%
% Inputs
%   dat        - Original data. Slice indices are ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))).
%   dim        - The dimension for the simultaneous slices, along which the encoding or decoding will be conducted.
%   fov_shift  - CAIPI FOV shift. abs(fov_shift) is the length of the DFTz encoding. Negative for inverse DFTz encoding, positive for forward DFTz encoding.
%   nz         - Number of simultaneous slices.
%   enc_or_dec - 'encode' or 'e': Encode the input 'dat'.
%                'decode' or 'd': Decode the input 'dat'.
%
% Output
%   dat        - Input 'dat' with forward or inverse DFTz encoded or decoded along 'dim'. Slice indices are ifftshift(-floor(nz/2):1:(ceil(nz/2)-1))).
%
% (c) Kangrong Zhu,     Stanford University     June 2013

if ~exist('dim', 'var') || isempty(dim)
    dim = ndims(dat);
end

if ~exist('fov_shift', 'var') || isempty(fov_shift)
    fov_shift = - size(dat, dim);
end

if ~exist('nz', 'var') || isempty(nz)
    nz = size(dat, dim);
end

if ~exist('enc_or_dec', 'var') || isempty(enc_or_dec)
    enc_or_dec = 'decode';
end

sz_orig = size(dat);
dims_all = 1 : length(sz_orig);
dims_non_slice = dims_all(dims_all ~= dim);                       % Dimensions except the slice dimension
dims_to_permute = [dim, dims_non_slice];                          % Matrix permutation order to make the slice dimension to be the 1st dimension
dims_to_permute_back = [2:dim, 1, dim+1:length(dims_to_permute)]; % Matrix permutation order to permute the data back to its orginal order

% Prepare data
dat = permute(dat, dims_to_permute);                              % Make the slice dimension to be the 1st dimension
sz_permuted = size(dat);
dat = reshape(dat, [sz_permuted(1), prod(sz_permuted(2:end))]);

% Calculate encoding matrix
encode_mtx = encode_dftz_mtx(fov_shift, nz);                         % Dim: [abs(fov_shift)(i.e. nkz), nz]

% Encode or decode
switch lower(enc_or_dec(1))
    case 'e'                                                      % Encode
        dat = encode_mtx * dat;
    case 'd'                                                      % Decode
        dat = encode_mtx \ dat;                                   % i.e. inv(encode_mtx)*dat
end
dat = reshape(dat, [size(dat,1), sz_permuted(2:end)]);
dat = permute(dat, dims_to_permute_back);

return