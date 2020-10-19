function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Freesurfer_Recon_Exams()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Freesurfer_Recon_Exams()
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
% Created:  03/09/12 by Don Hagler
% Last Mod: 12/23/14 by Don Hagler
%

arg_groups = {'FS'};
extra_tags = {'StudyInfo' 'RootDirs'...
  'STRUCT_T1type' 'STRUCT_BEMflag'...
  'old_rootdir' 'BEMtype' 'followup_flag' 'full_recon_flag'...
  'talairach_flag' 'nu_flag' 'nu3T_flag' 'labelV1_flag' 'surfsegedit_flag'...
  'batchname' 'logflag' 'logfile' 'touchonly_flag' 'forceflag'};
strip_flag = 0;
