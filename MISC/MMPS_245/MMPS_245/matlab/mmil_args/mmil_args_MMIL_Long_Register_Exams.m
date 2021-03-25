function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Long_Register_Exams()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Long_Register_Exams()
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
% Created:  06/30/11 by Vijay Venkatraman
% Last Mod: 04/11/14 by Don Hagler
%

arg_groups = {};
extra_tags = {'StudyInfo' 'RootDirs' 'batchname' 'qcflag',...
  'bindir','parmdir','dirprefix','imgtypes','ext','nobiasflag',...
  'baseflag','LongDTI_flag','forceflag'...
};
strip_flag = 0;
