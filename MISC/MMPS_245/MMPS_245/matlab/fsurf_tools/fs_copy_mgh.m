function fs_copy_mgh(fname_in,fname_out)
%function fs_copy_mgh(fname_in,fname_out)
%
% Purpose: copy an mgh file (possibly uncompressing)
%
% Required Input:
%   fname_in: full path name of input mgh/mgz file
%   fname_out: full path name of output mgh/mgz file
%
% Created:  02/14/10 by Don Hagler
% Last Mod: 03/10/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
keepsingle = 1;

[vol,M,mrparms] = fs_load_mgh(fname_in,[],[],[],keepsingle);
fs_save_mgh(vol,fname_out,M,mrparms);

