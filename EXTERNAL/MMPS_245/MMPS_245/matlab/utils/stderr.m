function y = stderr(x,w,dim)
%function y = stderr(x,w,dim)

if nargin < 2 || isempty(w), w = 0; end

if nargin < 3
    % The output size for [] is a special case when DIM is not given.
    if isequal(x,[]), y = NaN(class(x)); return; end

    % Figure out which dimension sum will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end
n = size(x,dim);

y = sqrt(var(x,w,dim)/n);

