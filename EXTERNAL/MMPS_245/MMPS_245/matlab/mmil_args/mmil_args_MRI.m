function [arg_groups,extra_tags,strip_flag] = mmil_args_MRI()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MRI()
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
% Created:  11/07/12 by Don Hagler
% Last Mod: 11/07/12 by Don Hagler
%

arg_groups = {'MRI'};
extra_tags = {};
strip_flag = 1;
