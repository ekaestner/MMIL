function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Process_Exams()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Process_Exams()
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
% Created:  05/19/11 by Don Hagler
% Last Mod: 12/18/12 by Don Hagler
%

arg_groups = {'DCM' 'STRUCT' 'DTI' 'BOLD'};
extra_tags = {'StudyInfo' 'RootDirs' 'batchname' 'qcflag'...
  'procstep' 'newflag' 'forceflag' 'kernelWidthMax' 'lambda2'...
  'kernelWidthMax_vec' 'lambda2_vec' 'multi_opt_flag',...
  'STRUCTflag','DTIflag','BOLDflag'};
strip_flag = 0;
