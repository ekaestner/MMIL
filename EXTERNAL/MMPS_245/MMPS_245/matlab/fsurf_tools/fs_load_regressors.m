function output = fs_load_regressors(fname)
%function output = fs_load_regressors(fname)
%
% Purpose: load regressors from csv file for fs_groupavg_glm
%
% Usage:
%  output = fs_load_regressors(fname, 'key1', value1,...);
%
% Required Input:
%  fname: full path name of csv file containing subject IDs
%     and regressor variables in additional columns
%     First column should have subject IDs, additional columns should contain
%       regressor values (e.g. group, test scores, etc.)
%     First row should be column headers, number of additional rows
%       is number of subjects included in analysis
%     Group columns should have binary (0 or 1) values to indicate membership
%
% Output:
%   output: structure containing SubjIDs, regressors, regnames
%
% Created:  10/19/11 by Don Hagler
% Last Mod: 10/21/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = [];
if ~mmil_check_nargs(nargin, 1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read regressor file

if ~exist(fname,'file'), error('file %s not found',fname); end;
try
  input_info = mmil_readtext(fname);
  SubjIDs = {input_info{2:end,1}};
catch
  error('failed to get info from regressors file %s -- check format',fname);
end;
regnames = input_info(1,2:end);
try
  regressors = cell2mat(input_info(2:end,2:end));
catch
  error('list of regressors in file %s must all be integers',fname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% package output

output.SubjIDs = SubjIDs;
output.regressors = regressors;
output.regnames = regnames;

