function StudyInfo = BOLD_MMIL_Check_StudyInfo(StudyInfo);
%function StudyInfo = BOLD_MMIL_Check_StudyInfo(StudyInfo);
%
% Required Parameters:
%   StudyInfo: StudyInfo struct from MMIL_Read_StudyInfo
%
% Created:  06/04/09 by Don Hagler (as BOLD_MMIL_Read_StudyInfo)
% Last Mod: 03/23/11 by Don Hagler
%

nsubs = length(StudyInfo);

% convert GLM_stim_fnames to nested cell array
%   with outer cell having one element for each snum
%   inner cell arrays containing variable number of stim files,
%     one for each non-baseline condition
if isfield(StudyInfo,'GLM_stim_fnames')
  for s=1:nsubs
    if ~isempty(StudyInfo(s).GLM_stim_fnames) &&...
       ischar(StudyInfo(s).GLM_stim_fnames)
      % split input string with ';' delimiter
      tmp_fnames = mmil_splitstr(deblank(StudyInfo(s).GLM_stim_fnames),';')';
      % split each of these with ' ' delimiter
      for j=1:length(tmp_fnames)
        tmp_fnames{j} = mmil_splitstr(tmp_fnames{j},' ')';
      end;
    else
      tmp_fnames = [];
    end;

    % reshape into linear cell array to avoid problems with mmil_args2parms
    %   in BOLD_MMIL_Analyze_Exam
    nstims = length(tmp_fnames);
    if nstims
      StudyInfo(s).GLM_num_conds = zeros(nstims,1,'uint16');
      StudyInfo(s).GLM_stim_fnames = {};
    else
      StudyInfo(s).GLM_num_conds = [];
      StudyInfo(s).GLM_stim_fnames = [];
    end;
    f = 1;
    for i=1:nstims
      if iscell(tmp_fnames{i})
        nconds = length(tmp_fnames{i});
        for j=1:nconds
          StudyInfo(s).GLM_stim_fnames{f} = tmp_fnames{i}{j};
          f=f+1;
        end;
      else
        nconds = 1;
        f=f+1;
      end;
      StudyInfo(s).GLM_num_conds(i) = nconds;
    end;
  end;
end;

