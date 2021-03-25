function Mm = c2matlab(Mc);

% function Mm = c2matlab(Mc);

Mm = Mc*[1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1]; % Convert from 0-based to 1-based indexing

return