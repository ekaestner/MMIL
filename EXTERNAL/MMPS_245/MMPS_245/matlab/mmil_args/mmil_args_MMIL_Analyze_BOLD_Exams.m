function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Analyze_DTI_Exams()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Analyze_DTI_Exams()
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
% Last Mod: 08/30/12 by Don Hagler
%

arg_groups = {'FP' 'FE' 'GLM'};
extra_tags = {'StudyInfo' 'RootDirs' 'batchname' 'qcflag' 'forceflag'...
  'infix' 'regFS_flag' 'mc_flag' 'mc_inter_flag' 'fstats_type' 'cxfstatsflag'...
  'paint_flag' 'paint_dist' 'paint_frac' 'paint_frac_flag' 'paint_dist_avg'...
  'paint_fract_avg' 'resamp_flag' 'force_repaint_flag' 'sphere_flag'};
strip_flag = 0;


