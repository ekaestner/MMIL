function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Get_StudyInfo()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Get_StudyInfo()
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
% Last Mod: 04/03/16 by Don Hagler
%

arg_groups = [];
extra_tags = {'VisitIDs' 'SubjIDs' ...
 'ignore_VisitInfo_flag' 'checkflag' 'qcflag' 'modality' ...
 'required_containers' 'required_rootdirs' ...
 'ico' 'numvec_tags' 'user' ...
 'QC_raw' 'QC_recon' 'QC_dv' 'QC_DTI' 'QC_BOLD' 'QC_PET' 'metaflag'};
strip_flag = 0;

