function n = bits2num(bits)
%
% Converts binary bits to a number.
%
% Input
%   bits - Binary bits in a vector. bits(1) is the most significant bit... bits(end) is the least significant bit.
%
% Output
%   n    - The converted number.
%
% (c) Kangrong Zhu,     Stanford University     April 2014

bitvals = unique(bits);
if (length(bitvals) > 2) || (bitvals(1)~=0) || (bitvals(2)~=1)
    error('Bit values must be either 0 or 1.');
end

nbits = length(bits);
n = 0;
for bit_idx = nbits-1 : -1 : 0
    this_bit = bits(nbits - bit_idx);
    n = n + this_bit * 2^bit_idx;
end