function dat = mux_shift(dat, shift_dim, shift_amt, dat_type)
%
% function dat = mux_shift(dat, shift_dim, shift_amt [, dat_type])
%
% Shift each slice in the specified dimension by a certain amount in the PE direction.
%
% Inputs : 
%   dat       - Input k-space or image-space data. Dim: [FE, PE, ...].
%   shift_dim - The dimension for the slice shifting.
%   shift_amt - A vector specifying the amount of shift which will be
%               applied to each slice. The m-th (m=0,1,2...) slice in the 
%               specified 'shift_dim' dimension will be shifted by 
%               shift_amt(m)*FOVy in the image space in the PE direction.
%   dat_type  - 'ksp'(Default): 'dat' is k-space data.
%               'im'          : 'dat' is image-space data.
%
% Output: 
%   dat       - K-space or image-space data after slice shifting, having
%                the same size and data type as the input 'dat'.
%
% (c) Kangrong Zhu,     Stanford University,    Oct 2012

if ~exist('dat_type', 'var') || isempty(dat_type)
    dat_type = 'ksp';
end

if strcmp(dat_type, 'im')
    dat = fft2c(dat);
else if ~strcmp(dat_type, 'ksp')
        error('The input ''dat_type'' must be either ''ksp'' or ''im''.');
    end
end

PE_DIM = 2;                   % The dimension for phase encoding in the data matrix
KEEP_ORIG_SZ = 1;             % Keep original matrix size when using 'repmat'

sz = size(dat);
ny = sz(PE_DIM);              % Number of ky lines.
nim = sz(shift_dim);          % Number of images to shift.
ndim = length(sz);            % Number of dimensions in 'dat'.

shift_amt = (shift_amt(:)).'; % Make shift_amt a row vector

% Basic phases to apply.
pha = 2*pi * (1:ny).';        % Applying a phase of -n*pha in the PE direction in k-space will shift the image by +n*FOVy in the image space.
pha = repmat(pha, [KEEP_ORIG_SZ, nim]) .* repmat(shift_amt, [ny, KEEP_ORIG_SZ]);

% Reshape and repeat to the same size as the input 'dat'.
reshape_sz = ones(1, ndim);
reshape_sz(PE_DIM) = ny;
reshape_sz(shift_dim) = nim;

rep_sz = sz;
rep_sz(PE_DIM) = KEEP_ORIG_SZ;
rep_sz(shift_dim) = KEEP_ORIG_SZ;

pha = repmat( reshape(pha, reshape_sz), rep_sz);

% Shift the images
dat = dat .* exp(1i * pha);

if strcmp(dat_type, 'im')
    dat = ifft2c(dat);
end

return;