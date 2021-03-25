function res = zpad(x, s1, s2, s3, s4, s5, s6)
%
% res = zpad(x, s1, s2)
% Zero pads a 2D matrix around its center to a size of [s1, s2].
%
% res = zpad(x, [s1, s2])
% Same as the previous example.
%
% res = zpad(x, s1, s2, s3)
% Zero pads a 3D matrix around its center.
%
% ...
%
% res = zpad(x, s1, s2, s3, s4, s5, s6)
% Zero pads a 6D matrix around its center.
%
% (c) Michael Lustig,           Stanford University
% Modified by Kangrong Zhu,     Stanford University      2012

if nargin < 2
    error('must have a target size');
end

s = [];
for sin = 1 : (nargin-1)
    s = [ s, eval(sprintf('s%d', sin)) ];
end

m = size(x);
if length(m) < length(s)
    m = [ m, ones(1, length(s)-length(m))];
end

if sum(m==s) == length(m)
    res = x;
    return;
end

res = zeros(s);

for n = 1:length(s)
    idx{n} = floor(s(n)/2)+1+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2);
end

% this is a dirty ugly trick
cmd = 'res(idx{1}';
for n = 2:length(s)
    cmd = sprintf('%s, idx{%d}', cmd, n);
end
cmd = sprintf('%s) = x;', cmd);
eval(cmd);



