function res = fftc(x, dim)
%
% res = fftc(x [, dim])
%
% Conduct orthonormal forward 1D FFT on the specified dimension of the input array.
% 
% Inputs:
%   x   - The array to conduct the orthonormal FFT on.
%   dim - The dimension to conduct the orthonormal FFT on. Default: the
%         last dimension of 'x'.
% 
% Output:
%   res - The resulting array after the orthonormal FFT.
%
% (c) Kangrong Zhu,     Stanford University     2011

if ~exist('dim','var')
    if iscolumn(x)
        dim = 1;
    else
        dim = length( size(x) );
    end
end

res = 1/sqrt(size(x,dim)) * fftshift(fft(ifftshift(x,dim),[],dim),dim);

