function  s = shuffle(n, seed)
%
% Shuffle index series 0,1,...n-1 into a random order. Used for random MICA.
% Corresponds to the shuffle function in the muxarcepi sequence.
%
% Inputs
%   n    - The index series to shuffle is 0,1,...n-1.
%   seed - Seed number for shuffling. 0,1,2...
%
% Output
%   s    - Shuffled index series.
%
% (c) Kangrong Zhu,     Stanford University     April 2014

s = (0 : 1 : (n-1)).';
rand_num = seed;
for shuffle_idx = 0 : n-2
    rand_num = rand_lcg_gcc(rand_num);
    idx = mod(rand_num, n-shuffle_idx);
    
    % Save value
    tmp = s(idx+1);
    
    % Move all unchosen values up in array
    s(idx+1 : n-shuffle_idx-1) = s(idx+2 : n-shuffle_idx);
    
    % Store chosen value at end of unchosen values
    s(n-shuffle_idx) = tmp;
end

return