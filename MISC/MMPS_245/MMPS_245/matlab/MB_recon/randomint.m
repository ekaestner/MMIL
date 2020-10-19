function  s = randomint(n, seed, RAND_MOD)
%
% Generate n random integers, each in the range 0,1,2...RAND_MOD-1.
% Used for MICA whose kz encoding is CAIPI kz encoding with random kz perturbations added.
% Corresponds to the randomint function in the muxarcepi sequence.
%
% Inputs
%   n        - Number of random integers to generate.
%   seed     - Seed number. 0,1,2...
%   RAND_MOD - Each generated random integer will be be in the range 0,1,2,...RAND_MOD-1.
%
% Output
%   s        - n generated random integers in a column.
%
% (c) Kangrong Zhu,     Stanford University     April 2014

if ~exist('RAND_MOD', 'var') || isempty(RAND_MOD)
    RAND_MOD = hex2dec('7fffffff');
end

s = zeros(n, 1);
rand_num = seed;
for idx = 0 : n-1
    rand_num = rand_lcg_gcc(rand_num);
    s(idx + 1) = mod(rand_num, RAND_MOD);
end

return