function im = sos(ic, dim)
%
% im = sos(ic [, dim])
%
% Square root of sum of squares(SOS) reconstruction of multicoil images.
% 
% Inputs:
%   ic  - Single coil images. 
%   dim - The dimension for coil in matrix 'ic'. Default: the last dimension of 'ic'.
%
% Output:
%   im  - The SOS combined images.
%
% (c) Kangrong Zhu,     Stanford University     2011

if ~exist('dim','var')
    dim = length( size(ic) );
end

im = sqrt( sum( ic.*conj(ic), dim) );