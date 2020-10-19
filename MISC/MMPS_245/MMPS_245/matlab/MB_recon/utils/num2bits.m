function bits = num2bits(n, nbits)
%
% Converts a number to binary bits. Keeps the least significant 'nbits' bits.
%
% Inputs
%   n     - A number.
%   nbits - Number of bits to convert to. Default: 32.
%
% Output
%   bits  - Binary bits in a vector of length nbits. bits(1) is the most significant bit... bits(nbits) is the least significant bit.
%
% (c) Kangrong Zhu,     Stanford University     April 2014

if ~exist('nbits', 'var') || isempty(nbits)
    nbits = 32;
end

bits = zeros(nbits, 1);
for bit_idx = nbits-1 : -1 : 0
    this_bit = mod(floor(n/2^bit_idx), 2); % Use modulus 2 to account for floor(n/2^bit_idx)>1
    bits(nbits - bit_idx) = this_bit;
end