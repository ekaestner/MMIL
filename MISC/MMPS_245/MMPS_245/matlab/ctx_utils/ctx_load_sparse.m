function [ctx_vol,volsz] = ctx_load_sparse(fname)
%function [ctx_vol,volsz] = ctx_load_sparse(fname)
%
% created:   04/09/12 by Don Hagler
% last mod:  04/09/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

[vol,M,volsz] = mmil_load_sparse(fname);
ctx_vol = ctx_mgh2ctx(vol,M);

