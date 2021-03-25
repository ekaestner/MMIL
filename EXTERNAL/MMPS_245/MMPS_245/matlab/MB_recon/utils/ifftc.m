function res = ifftc(x, dim)
%
% res = ifftc(x [, dim])
%
% Conduct orthonormal inverse 1D FFT on the specified dimension of the input array.
%
% Inputs:
%   x   - The array to conduct orthonormal iFFT on.
%   dim - The dimension to conduct orthonormal iFFT on. Default: the 
%         last dimension of 'x'.
% 
% Output:
%   res - The resulting array after the orthonormal iFFT.
%
% (c) Kangrong Zhu,     Stanford University     2011

if ~exist('dim','var')
    if iscolumn(x)
        dim = 1;
    else
        dim = length( size(x) );
    end
end

res = sqrt(size(x,dim)) * fftshift(ifft(ifftshift(x,dim),[],dim),dim);

