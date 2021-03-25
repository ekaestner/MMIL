function [arg_groups,extra_tags,strip_flag] = mmil_args_DTI_MMIL_Average_Atlas_Fibers()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_DTI_MMIL_Average_Atlas_Fibers()
%
% Purpose: for use with MMIL_Args in selecting context relevant parameters
%   from parameter struct arrays (e.g. parms, ProjInfo)
%
% Output:
%   arg_groups: cell array of fieldname prefixes (e.g. 'BOLD', 'PROC', etc.)
%     specifying group of parameter names
%   extra_tags: cell array of parameter names
%   strip_flag: [0|1] whether arg_group prefixes should be stripped from
%     final parameter names
%
% Created:  07/20/11 by Vijay Venkatraman
% Last Mod: 07/20/11 by Vijay Venkatraman
%

arg_groups = {};
extra_tags = {'ATLDIR',...
'tensor_smooth_sigma',... 
'countflag',...
'first_only_flag',...
'min_tensor_count',...
'smf',...
'forceflag'
};
strip_flag = 0;
