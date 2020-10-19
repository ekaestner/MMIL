function dim = mmil_nonsingleton_dim(x)
%function dim = mmil_nonsingleton_dim(x)
%
% Purpose: determine the first non-singleton dimension of a matrix
%
% Created:  10/22/12 by Don Hagler
% Last Mod: 10/22/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

dims = size(x);
ind = find(dims~=1);
if isempty(ind)
  dim = length(dims);
else
  dim = min(ind);
end;

