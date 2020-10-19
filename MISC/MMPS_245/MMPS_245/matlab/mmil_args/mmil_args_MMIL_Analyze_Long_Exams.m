function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Analyze_Long_Exams()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Analyze_Long_Exams()
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
  'baseflag','sphsmoothsteps','mask_midbrain_flag','outdir',...
  'aseg_flag','aseg_roigroups_flag','subhippo_flag','aparc_flag',...
  'nobiasflag','forceflag'};
strip_flag = 0;

