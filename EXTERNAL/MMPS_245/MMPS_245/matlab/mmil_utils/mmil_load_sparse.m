function [vol,M,volsz] = mmil_load_sparse(fname)
%function [vol,M,volsz] = mmil_load_sparse(fname)
% 
% Purpose:
%   Loads sparse matrix and converts it to output volume 
%
% Required Parameters:
%   fname: Sparse input filename
%
% Output: 
%  vol: output volume 
%  M: M matrix of the volume
%  volsz: vector of number of elements for each dimension
%
% Created:  04/08/11 by Vijay Venkatraman
% Last Mod: 05/21/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
sparse_vol = [];

if ~exist(fname,'file'), error('file %s not found',fname); end;
load(fname);
if isempty(sparse_vol), error('file %s is missing sparse_vol',fname); end;

volsz = sparse_vol.volsz;
vol = zeros(volsz);
vol(sparse_vol.indices) = sparse_vol.values;
M = sparse_vol.M;

return;
