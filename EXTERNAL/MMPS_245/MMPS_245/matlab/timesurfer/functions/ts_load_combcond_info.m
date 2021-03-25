function [neweventcodes,combinations,info] = ts_load_combcond_info(fname)
%function [neweventcodes,combinations,info] = ts_load_combcond_info(fname)
%
% Purpose: load info from csv file for ts_combine_conditions
%
% Required Input:
%  fname: full path name of csv file
%   two columns headers required: 'neweventcodes' and 'combinations'
%   additional columns allowed for additional information
%
% Output:
%   neweventcodes: vector of event codes for new combination conditions
%   combinations: cell array of condition combination strings
%     e.g. '3-5' or '1+2+3'
%     use single quotes to prevent excel or open office from converting to dates
%   info: struct array with all columns from fname as fields
%
% See also: ts_combine_conditions
%
% Created:  10/21/11 by Don Hagler
% Last Mod: 10/24/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neweventcodes = [];
combinations = [];
info = [];
if ~mmil_check_nargs(nargin, 1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read file

info = mmil_csv2struct(fname);

if ~isfield(info,'neweventcodes')
  error('missing "neweventcodes" column in %s',fname);
end;
if ~isfield(info,'combinations')
  error('missing field "combinations" in %s',fname);
end;

neweventcodes = {info.neweventcodes};
combinations = {info.combinations};

if any(cellfun(@isempty,neweventcodes))
  error('%s has rows with empty cell for neweventcodes');
end;
if any(cellfun(@isempty,combinations))
  error('%s has rows with empty cell for combinations');
end;

neweventcodes = cell2mat(neweventcodes);

combinations = regexprep(combinations,'''','');

