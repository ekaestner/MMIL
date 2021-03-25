function [sub1,sub2,sub3,ind] = FindMax(M)
[val,ind] = max(M(:));
[sub1,sub2,sub3] = ind2sub(size(M),ind);
