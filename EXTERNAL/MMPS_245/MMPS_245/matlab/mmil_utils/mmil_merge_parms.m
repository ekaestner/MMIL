function parms = mmil_merge_parms(parms1,parms2)
%function parms = mmil_merge_parms(parms1,parms2)
%
% Purpose:
%   Merge two struct arrays into one struct array
%   If a field is present in both input structs, the
%   associated value in parms1 will be saved in parms.
%
% Required Parameters:
%   parms1: a struct with one or more fields
%   parms2: a struct with one or more fields
%
% Optional Parameters:
%   args: a named_argument-value pair list (1-D cell array).
%
% Created:  10/21/11 by Don Hagler
% Last Mod: 10/18/12 by Don Hagler
%

parms = struct();
if ~mmil_check_nargs(nargin,2), return; end;

args1 = mmil_parms2args(parms1);
args2 = mmil_parms2args(parms2);
args = mmil_merge_args(args1,args2);
parms = mmil_args2parms(args,[],0);

