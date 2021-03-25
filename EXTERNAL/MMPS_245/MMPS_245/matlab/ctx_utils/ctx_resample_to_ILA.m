function ctx_vol_res = ctx_resample_to_ILA(ctx_vol,M_reg,interpm,bclamp)
%function ctx_vol_res = ctx_resample_to_ILA(ctx_vol,M_reg,[interpm],[bclamp])
%
% like resample_to_LIA, but resample to ILA instead of LIA
%
%   interpm : 0: : Nearest 1:Linear 2:cubic (default)
%             3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%   bclamp : set negative value to zero default = true
%
% Created:  03/06/13 by Don Hagler
% Last Mod: 03/06/13 by Don Hagler
%

ctx_vol_res = [];
if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('interpm','var') | isempty(interpm), interpm = 2; end
if ~exist('bclamp','var') | isempty(bclamp), bclamp = 1; end

size_out = [256 256 256];
M_cor_ideal =  [ 0    -1     0   129;...
                 0     0     1  -129;...
                -1     0     0   129;...
                 0     0     0     1];

[vol,M] = ctx_ctx2mgh(ctx_vol);
[vol_res,M_res] = mmil_resample_vol(vol,M,...
  'M_ref',M_cor_ideal,...
  'nvox_ref',size_out,...
  'M_reg',M_reg,...
  'interpm',interpm,...
  'bclamp',bclamp);
ctx_vol_res = ctx_mgh2ctx(vol_res,M_res);

