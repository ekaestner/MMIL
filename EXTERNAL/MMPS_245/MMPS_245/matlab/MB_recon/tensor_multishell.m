function [grad_vec ]=tensor_multishell_NEW(ndirs)

% The location of your tensor file needs to be changed
A = importdata('/home/kuperman/matlab/MB_recon/tensor_ABCD.dat');

sz = size(find((A.data(:,1)==ndirs)));

if sz(1,1) ~= 1
    error('There are multiple entries within tensor.dat for %d directions',ndirs);
end

grad_vec = zeros(ndirs,3);
start = find((A.data(:,1)==ndirs)) + 1;

c = 0;
for i = 1:ndirs
    grad_vec(i,:) = A.data((3*c)+start:(3*c)+start+2,1)';
    c = c + 1;
end

end

