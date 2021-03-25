function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_IcoResamp_FSRecon_Exams()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_IcoResamp_FSRecon_Exams()
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
% Created:  08/08/11 by Don Hagler
% Last Mod: 08/08/11 by Don Hagler
%

arg_groups = {};
extra_tags = {'StudyInfo' 'RootDirs' 'FS_version' 'ico' 'batchname'...
  'ProjID' 'logflag' 'logfile' 'check_complete_flag' 'forceflag'...
  'hemilist' 'required_rootdirs' 'source_typestr','dest_typestr'};
strip_flag = 0;

