function cond_info = rc_get_cond_info_from_prefix(prefix)
%function cond_info = rc_get_cond_info_from_prefix(prefix)
%
% Created:  06/11/09 by Don Hagler
% Last Mod: 11/15/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

cond_info = [];
matfile = sprintf('matfiles/%s_ret_forward.mat',prefix);
if ~exist(matfile)
  error('file %s not found\n',matfile);
end;

load(matfile)
cond_info = mmil_getfield(retforward.retmap,'cond_info',[]);

