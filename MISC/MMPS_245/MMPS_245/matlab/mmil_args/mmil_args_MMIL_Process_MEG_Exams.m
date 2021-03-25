function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Process_MEG_Exams()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Process_MEG_Exams()
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
% Created:  03/23/11 by Don Hagler
% Last Mod: 03/23/11 by Don Hagler
%

arg_groups = {'RAW' 'PROC'};
extra_tags = {'StudyInfo' 'RootDirs' 'batchname' 'qcflag'...
  'procstep' 'newflag' 'forceflag'};
strip_flag = 0;
