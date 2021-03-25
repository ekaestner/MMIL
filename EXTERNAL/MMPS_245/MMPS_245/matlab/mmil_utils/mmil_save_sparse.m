function mmil_save_sparse(vol,fname,M)
%function mmil_save_sparse(vol,fname,M)
% 
% Purpose:
%   Convert from mgh volume to sparse matrix (stores values and indices)
%
% Required Parameters:
%  vol: volume in mgh format
%  M: M matrix of the input volume
%  fname: Output filename for the sparse matrix (.mat)
%
% Created:  04/08/11 by Vijay Venkatraman
% Last Mod: 05/21/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

if isempty(vol) | isempty(M) | isempty(fname)
  error('vol or M or fname is empty');
end;

sparse_vol = [];
sparse_vol.indices = find(vol~=0);
sparse_vol.values = vol(sparse_vol.indices);
sparse_vol.volsz = size(vol);
sparse_vol.M = M;
save(fname,'sparse_vol');

return;
