function rand_num = rand_lcg_gcc(seed)
%
% Linear congruential generator (LCG) which yields a sequence of randomized
% numbers calculated with a linear equation.
% Used for generating random numbers for MICA.
% Corresponds to the rand_lcg_gcc function in the muxarcepi sequence.
%
% Input
%   seed     - Seed number. 0,1,2...
%
% Output
%   rand_num - Generated random number. 0,1,2,...,RAND_MOD(i.e. hex2dec('7fffffff')).
%

% The following code mimics this C code:
%
% #define RAND_MOD (0x7fffffff)
% #define RAND_A (1103515245)
% define RAND_C (12345)
% int rand_lcg_gcc (int seed){
%   long int rand_num;
%   rand_num = (seed * RAND_A + RAND_C) & RAND_MOD;
%  return (rand_num);
% }


RAND_MOD = int32(hex2dec('7fffffff'));
RAND_A = int64(1103515245);
RAND_C = int64(12345);

rand_num = bitand(int32(mod((int64(seed) * RAND_A + RAND_C), 2^31)), RAND_MOD);

return


% OLD version, using the FP toolbox.
%
% To compare: s=zeros(2,1e5); for(i=0:1e5), s(:,i+1) = [rand_lcg_gcc(i), rand_lcg_gcc_ORIG(i))]; end
% Use fixed-point object with embedded.numerictype T and embedded.fimath F. To mimic the c language.
%T = numerictype('WordLength', 32, 'FractionLength', 0, 'Signed', 1);
%F = fimath('MaxProductWordLength', 32, 'OverflowMode', 'Wrap', 'ProductMode', 'KeepLSB', 'ProductWordLength', 32, 'SumMode', 'KeepLSB', 'SumWordLength', 32);
%
%seed   = fi(seed, T, F);
%RAND_A = fi(1103515245, T, F);
%RAND_C = fi(12345, T, F);
%
%% Generate random number
%rand_num = seed * RAND_A + RAND_C; % rand_num is a fixed-point object
%rand_num = double(rand_num);       % rand_num is double
%bits = num2bits(rand_num, 32);     % Get least significant 32 bits
%bits = bits(2 : 32);               % Keep the least significant 31 bits. Equivalent to (seed * RAND_A + RAND_C) & RAND_MOD in EPIC. Where RAND_MOD = hex2dec('7fffffff');
%rand_num = bits2num(bits);
%
%return
