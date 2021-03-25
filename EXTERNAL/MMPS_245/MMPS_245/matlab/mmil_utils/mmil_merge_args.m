function args = mmil_merge_args(arglist1,arglist2)
%function args = mmil_merge_args(arglist1,arglist2)
%
% Purpose:
%   Merge two named_argument-value pair list into one list.
%   If a named_argument is present in both input lists, the
%   associated value in arglist1 will be saved in args.
%
% Required Parameters:
%   arglist1 - a named_argument-value pair list (cell array).
%   arglist2 - a named_argument-value pair list (cell array).
%
% Optional Parameters:
%   args - a named_argument-value pair list (1-D cell array).
%   
%
% Created:  02/02/11 by C Roddey
% Prev Mod: 09/14/12 by Don Hagler
% Last Mod: 08/29/16 by Don Hagler
%

args = cell(0);

if ~mmil_check_nargs(nargin,2), return; end;

if rem(length(arglist1),2),
  error('length of ARGLIST1 [%d] must be even',length(arglist1));
end
if rem(length(arglist2),2),
  error('length of ARGLIST2 [%d] must be even',length(arglist2));
end

argnames1 = mmil_rowvec(arglist1(1:2:end));
argnames2 = mmil_rowvec(arglist2(1:2:end));

args = mmil_rowvec(arglist1);

[addl_names2,ind] = setdiff(argnames2,argnames1);

if ~isempty(ind)
  inds = sort(mmil_rowvec([ind*2-1 ind*2]));  % get inds for arglist2 name/value pairs
  args = [args arglist2(inds)];
end

if rem(length(args),2),
  error('failed to merge arglists');
end
