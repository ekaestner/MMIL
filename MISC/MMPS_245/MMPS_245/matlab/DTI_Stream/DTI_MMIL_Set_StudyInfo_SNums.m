function StudyInfo = DTI_MMIL_Set_StudyInfo_SNums(StudyInfo,snums_flag)
%function StudyInfo = DTI_MMIL_Set_StudyInfo_SNums(StudyInfo,[snums_flag])
%
% Usage:
%  StudyInfo = DTI_MMIL_Set_StudyInfo_SNums(StudyInfo,snums_flag)
%
% Required Parameters:
%  StudyInfo: struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%      may contain these fields
%        DTIScanNums, DTIScanNums2, DTIScanNums3
%
% Optional Parameters:
%  snums_flag: [0|1|2|3] which set of "snums" to use
%     0: use all available scans
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%   NOTE: subjects with values of -1 are excluded
%     subject with empty ScanNums will use all available
%    {default = 0}
%
% Output:
%   StudyInfo: struct array of study information
%     DTIScanNums field set according to snums_flag
%
% Created:  11/18/09 by Don Hagler
% Last Mod: 08/30/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,1)) return; end;
if isempty(StudyInfo), return; end;
if ~exist('snums_flag','var') || isempty(snums_flag), snums_flag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(StudyInfo)
  StudyInfo(i).DTI = 1;
  switch snums_flag
    case 0
      snums = [];
    case 1
      if ~isfield(StudyInfo(i),'DTIScanNums')
        snums = [];
      else
        snums = StudyInfo(i).DTIScanNums;
        % if -1, exclude this subject
        if ismember(-1,snums), StudyInfo(i).DTI = 0; end;
      end;
    case 2
      if ~isfield(StudyInfo(i),'DTIScanNums2')
        snums = [];
      else
        snums = StudyInfo(i).DTIScanNums2;
        % if -1, exclude this subject
        if ismember(-1,snums), StudyInfo(i).DTI = 0; end;
      end;
    case 3
      if ~isfield(StudyInfo(i),'DTIScanNums3')
        snums = [];
      else
        snums = StudyInfo(i).DTIScanNums3;
        % if -1, exclude this subject
        if ismember(-1,snums), StudyInfo(i).DTI = 0; end;
      end;
  end;
  StudyInfo(i).DTIScanNums = snums;
end;


