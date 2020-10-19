function result=mmil_isrelative(fname)
%function result=mmil_isrelative(fname)
%
% Required Input:
%   fname: full or relative path file name
%
% Returns true if fname is a relative path
%   or false if fname is a full path
%
%
% Created:  04/26/10 by Don Hagler
% Last Mod: 04/26/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
result = isempty(regexp(fname,'^/'));

