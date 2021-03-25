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
% Created:  06/30/11 by Vijay Venkatraman
% Last Mod: 02/04/12 by Don Hagler
%

arg_groups = {};
extra_tags = {'StudyInfo' 'RootDirs' 'batchname' 'qcflag',...
  'required_containers', 'required_rootdirs', 'modality',...
  'LongDTI_flag','xcg_flag','masksf_flag','fibers','dirprefix','T1type',...
  'outext','fiberdir_resT1','fiber_infix','forceflag','DTI_snums_flag',...
  'DTI_nob0_flag','DTI_min_ndirs','DTI_min_bval',...
  'DTI_flex_flag','DTI_min_nb0','DTI_revflag',...
  'DTI_infix', 'DTI_DT_regT1flag','forceflag'};
strip_flag = 0;
