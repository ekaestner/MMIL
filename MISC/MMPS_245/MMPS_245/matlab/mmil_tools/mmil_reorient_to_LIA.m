function mmil_reorient_to_LIA(fname_in,fname_out,trans_vec)
%function mmil_reorient_to_LIA(fname_in,fname_out,[trans_vec])
%
% Purpose: reorients a volume to LIA without resampling
%
% Required:
%   fname_in: input file name (mgh format)
%   fname_out: output file name (mgh format)
%
% Optional:
%   trans_vec: vector of [x,y,z] translation in vox2ras matrix
%     if empty, calculate from input vox2ras
%     {default = []}
%
% Early Mod: 07/02/08 by Don Hagler
% Last Mod:  04/23/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,2)) return; end;

if ~exist('trans_vec','var'), trans_vec=[]; end;

if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;

vol = ctx_load_mgh(fname_in);
vol = ctx_reorient_to_LIA(vol,trans_vec);
ctx_save_mgh(vol,fname_out);

