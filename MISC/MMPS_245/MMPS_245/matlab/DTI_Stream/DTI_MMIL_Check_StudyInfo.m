function [StudyInfo,flags] = DTI_MMIL_Check_StudyInfo(StudyInfo,RootDirs,varargin)
%function StudyInfo = DTI_MMIL_Check_StudyInfo(StudyInfo,RootDirs,[options])
%
% Usage:
%  [StudyInfo,flags]=DTI_MMIL_Check_StudyInfo(StudyInfo,RootDirs,'key1', value1,...);
%
% Required Parameters:
%  StudyInfo: struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%      must contain these fields:
%        SubjID
%        StudyDate
%      may contain these fields
%        proc_dti
%        fsurf
%      if proc_dti or fsurf are unspecified, will look for Containers
%       with SubjID and StudyDate (will choose first one if more than one)
%      if empty, use all subjects found in RootDirs.proc_dti
%  RootDirs: struct that must contain the following fields:
%       proc_dti, fsurf
%    and may contain the following fields:
%       home, batch, orig, raw, proc_dti, fsurf, fsico, long
%    these specify the full paths of root directories containing data containers
%
% Optional Parameters:
%  'snums_flag': [0|1|2|3] which set of "snums" to use
%     0: use all available scans
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%    {default = 0}
%   'qcflag': [0|1] whether to exclude subjects with QC=0
%     {default = 0}
%   'DTI_flag': [0|1[ whether to exclude subjects without DTI data
%     {default = 1}
%
% Output:
%   StudyInfo: struct array of study information
%     created from Container names (if not supplied as input)
%     or filtered by checking for valid scans and rejected subjects (QC=0)
%   flags: struct array containing flags indicating whether certain fields
%     are included in StudyInfo ('group', 'visitcode','age', and 'statsflag')
%
%
% Created:  08/11/08 by Don Hagler
% Rcnt Mod: 03/21/10 by Don Hagler
% Last Mod: 10/05/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums_flag',0,[0:3],...
  'qcflag',false,[false true],...
  'DTI_flag',true,[false true],...
...
  'required_containers',{'proc_dti'},[],...
});

flags = [];
flags.group = 0;
flags.visitcode = 0;
flags.age = 0;
flags.stats = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_StudyInfo');
[StudyInfo,RootDirs] = MMIL_Check_StudyInfo(StudyInfo,RootDirs,args{:});

if isempty(StudyInfo), return; end;
if isfield(StudyInfo(1),'Group'), flags.group = 1; end;
if isfield(StudyInfo(1),'VisitCode'), flags.visitcode = 1; end;
if isfield(StudyInfo(1),'Age'), flags.age = 1; end;
if isfield(StudyInfo(1),'StatsFlag'), flags.stats = 1; end;

% choose snums based on snums_flag
StudyInfo = DTI_MMIL_Set_StudyInfo_SNums(StudyInfo,parms.snums_flag);

% exclude subjects without DTI data
if parms.DTI_flag
  StudyInfo = StudyInfo(find(cell2mat({StudyInfo.DTI})));
end;


