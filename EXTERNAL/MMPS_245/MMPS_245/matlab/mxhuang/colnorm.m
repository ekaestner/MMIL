function [nrm,Anrm] = colnorm(A);
%COLNORM - Calculate L2 norm for each column of a matrix and normalize
% function [nrm,Anrm] = colnorm(A);
% calculate the Euclidean norm of each COLUMN in A, return as
% a row vector with same number of cols as A.
% Optionally, return A with each column now normalized


[m,n] = size(A);

if(m>1),			% multiple rows
  nrm = sqrt(sum([A.*conj(A)]));
else				% A is row vector
  nrm = abs(A);			% just return mag of each column
end

if(nargout > 1),
  ndx = find(nrm>0);		% any zero norm?
  Anrm = zeros(size(A));
  % normalize any non-zero columns
  Anrm(:,ndx) = A(:,ndx) ./ nrm(ones(m,1),ndx);
end

return
