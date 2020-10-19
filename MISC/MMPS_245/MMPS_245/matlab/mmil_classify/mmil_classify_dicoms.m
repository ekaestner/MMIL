function SeriesInfo = mmil_classify_dicoms(SeriesInfo,rules)
%function SeriesInfo = mmil_classify_dicoms(SeriesInfo,rules)
%
% Purpose: gather information from dicom header
%   and determine scan type
%
% Required Input:
%   SeriesInfo: struct array containing information about
%     multiple dicom image series returned from mmil_unpack_dicoms
%
% Optional Input:
%   rules: struct array containing rules for each series in fname_rules
%     organized by operation
%
%   Operations:
%     'match' : use regexp to match a part of the string
%     'exact' : use strcmp to match the entire string or == to match a value
%     'expr'  : use eval to evaluate a logical expresssion (e.g. >0, <100, >=5)
%               NOTE: value from SeriesInfo always on left side of expression
%     'set'   : set SeriesInfo field with this value
%
%   e.g. rules(1).SeriesType = 'MPR'
%        rules(1).match.SeriesDescription = 'MPRAGE'
%        rules(1).exact.SequenceName = 'EFGRE3D'
%        rules(1).expr.FlipAngle = '<10'
%        rules(1).set.rfrxcoiltype = 'BODY'
%
% Output:
%   SeriesInfo: will contain SeriesType field plus
%     additional information extracted from dicom headers
%     if SeriesInfo(s).ignore = 1, not valid dicom file series
%
%   NOTE: output SeriesInfo will be sorted by
%     StudyDate, StudyTime, SeriesDate, and SeriesTime
%
% Created:  05/29/11 by Don Hagler
% Last Mod: 11/16/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('rules','var'), rules = []; end;

% get info from dicom headers
SeriesInfo = mmil_init_serinfo(SeriesInfo);

% modify info using additional info and private tags
SeriesInfo = mmil_modify_serinfo(SeriesInfo);

% determine scan type based on codified rules
SeriesInfo = mmil_classify_by_code(SeriesInfo);

% classify according to user-supplied rules
if ~isempty(rules)
  SeriesInfo = mmil_classify_by_rules(SeriesInfo,rules);
end;

% for BOLD scans, determine slice timing patterns
SeriesInfo = mmil_slicetiming_BOLD(SeriesInfo);

% for DTI scans, get diffusion directions
SeriesInfo = mmil_qmat_DTI(SeriesInfo);

% sort series by StudyDate, StudyTime, SeriesDate, and SeriesTime
SeriesInfo = mmil_sort_serinfo(SeriesInfo);

