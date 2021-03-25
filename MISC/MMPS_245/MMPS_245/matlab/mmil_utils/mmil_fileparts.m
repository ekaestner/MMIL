function s = mmil_fileparts(f,fp)
%function s = mmil_fileparts(f,fp)
%
% Purpose: 
%   Works like FILEPARTS, but allows you to 
%   specify WHICH filepart(s) to return.
%
%   Requested 'parts' are concatenated together
%   'appropriately' (if no 'appropriate' concatenation is
%   found, then the parts are just concatenated in the order
%   that they were requested by param fp)
%
% Input Arguments:
%   f   : path to parse (string)
%   fp  : path parts to return (cell array)
%     'path', 'name', 'ext', or 'versn'
%
% Output Arguments:
%   s   : output string of concatenated parts
%
% See also: fileparts
%
% Created:  07/15/07 by Ben Cipollini
% Rcnt Mod: 08/21/07 by Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

if ~iscell(fp), fp = {fp}; end;

[pathstr name ext]   = fileparts(f);
s = '';

for i=1:length(fp)
  switch (fp{i})
    case {'path','pathstr'}
                    s   = fullfile(s,pathstr);
    case 'name',    s   = [s name];
    case 'ext',     s   = [s ext];
    otherwise
      error('unknown filepart: %s', fp);
  end;
end;
