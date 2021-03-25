function [arrout, numelout] = checkrange(arr, amin, amax, atype, debug)
%
% Check values in array are within range, and remove ones that are not.
% 
% Inputs:
%   arr      - The input array.
%   amin     - Minimum for the range(inclusive).
%   amax     - Maximum for the range(inclusive).
%   atype    - A string specifying the type of the array. Used for
%              displaying the out-of-range elements in the array.
%              Only needed when debug==true.
%   debug    - True: Print out some messages for debugging;
%              False(Default): No messages.
%
% Output:
%   arrout   - The output array with out-of-range elements removed.
%   numelout - The number of elements in the output array 'arrout';
%
% (c) Brian Hargreaves,         Stanford University     2008
% Modified by Kangrong Zhu,     Stanford Univeristy     2012

if ~exist('debug', 'var') || isempty('debug')
    debug = false;
end

arrout = arr(arr >= amin);
arrout = arrout(arrout <= amax);

numelout = numel(arrout);

if debug && (numelout < numel(arr))
    arr = arr(:);
    f1 = find(arr < amin);
    f2 = find(arr > amax);
    fprintf('   %d %s out of range: %s \n', ...
        length(f1)+length(f2), atype, num2str([arr(f1); arr(f2)].'));
    fprintf('   Removed those %s from the array.\n', atype);
end;

return;
