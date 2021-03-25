function args = MMIL_Args(parms,context)
%function args = MMIL_Args(parms,[context])
%
% Purpose: return arguments (cell array of 'key',value pairs)
%   from parameter structure, depending on context
%
% Required Input:
%   parms: matlab struct (e.g. parms, ProjInfo, StudyInfo)
%
% Optional Input:
%   context: name of function needing arguments from parms
%     must be a matching function called mmil_args_{context}
%     Allowed values:
%      'MMIL_Check_StudyInfo'
%      'MMIL_Get_StudyInfo'
%      'MMIL_Check_ProjID'
%      'MMIL_Get_QCInfo'
%      'MMIL_Process_Exams'
%      'MMIL_PreProcess_Exams'
%      'MMIL_Classify_Exams'
%      'MMIL_Freesurfer_Recon_Exams'
%      'MMIL_IcoResamp_FSRecon_Exams'
%      'MMIL_Create_BEM_Exams'
%      'MMIL_Analyze_BOLD_Exams'
%      'MMIL_RetFit_BOLD_Exams'
%      'BOLD_MMIL_Ready_RetFit'
%      'BOLD_MMIL_GLM_ROI_GroupAvg'
%      'DTI_MMIL_Import_DTIStudio_Fibers_Exams' 
%      'DTI_MMIL_WarpToAtlas_Exams'
%      'DTI_MMIL_Average_Atlas_Fibers'
%      'MMIL_Long_Setup_Exams'
%      'MMIL_Long_Register_Exams'
%      'MMIL_Process_MEG_Exams'
%      'MMIL_Import_MEG'
%      'MMIL_Process_MEG'
%      'MMIL_Analyze_MEG_Exams'
%      'MEG_MMIL_dSPM'
%      'MEG_MMIL_dSPM_ROI'
%      'MEG_MMIL_RCSE'
%      'MEG_MMIL_dSPM_RCSE'
%      'MEG_MMIL_RCSE_GroupAvg'
%      'MEG_MMIL_GroupRCSE'
%      'MEG_MMIL_dSPM_ROI_GroupAvg'
%      'rc_fit_wforms'
%      'rc_analyze_wforms'
%      'rc_RCSE'
%      'ts_process_fif_data'
%     If empty, return all arguments
%     {default = []}
%
% Created:  02/17/11 by Don Hagler
% Last Mod: 05/23/12 by Vijay Venkatraman
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('context','var'), context = []; end;
args = [];

if isempty(parms)
  return;
end;
if isempty(context)
  args = mmil_parms2args(parms);
  return;
end;

try
  [arg_groups,extra_tags,strip_flag] = eval(['mmil_args_' context]);
catch
  error('invalid context %s, missing function mmil_args_%s',...
    context,context);
end;

inds = [];
all_tags = fieldnames(parms);
for g=1:length(arg_groups)
  % find tags that start with arg_group
  tmp = sprintf('^%s_',arg_groups{g});  
  ind = find(~cellfun(@isempty,regexp(all_tags,tmp)));
  if ~isempty(ind), inds = union(inds,ind); end;
end;
tags = all_tags(inds);
tags = cat(1,tags,extra_tags');

if isempty(tags), return; end;

args = mmil_parms2args(parms,tags);

if strip_flag
  nparms = length(args)/2;
  for n=1:nparms
    tmp_name = args{(2*(n-1)) + 1};
    ind = find(~cellfun(@isempty,regexp(tmp_name,arg_groups)));
    if ~isempty(ind)
      tmp = sprintf('^%s_',arg_groups{ind});  
      tmp_name = regexprep(tmp_name,tmp,'');
      args{(2*(n-1)) + 1} = tmp_name;
    end;
  end;
end;

