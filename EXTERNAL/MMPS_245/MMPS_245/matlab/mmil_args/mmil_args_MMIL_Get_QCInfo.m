function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Get_QCInfo()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_MMIL_Get_QCInfo()
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
% Last Mod: 12/31/12 by Don Hagler
%

arg_groups = [];
extra_tags = {'QC_raw' 'QC_recon' 'QC_dv' 'QC_DTI' 'QC_BOLD'...
              'QC_PET' 'metaflag'};
strip_flag = 0;

