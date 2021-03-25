function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Compile_Long_Analysis()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Compile_Long_Analysis()
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
% Created:  04/11/14 by Don Hagler
% Last Mod: 04/11/14 by Don Hagler
%

arg_groups = {};
extra_tags = {'StudyInfo' 'RootDirs' 'batchname' 'qcflag',...
              'outdir','baseflag','regtype','aseg_flag',...
              'subhippo_flag','aseg_roigroups_flag','aparc_flag'};
strip_flag = 0;


