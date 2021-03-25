function out = mmil_splitstr( str, sep, max_iter )
% function out = mmil_splitstr( str, sep )
%
% Purpose: 
%   cuts a string into parts based on sep,
%   returns a cell array of the parts.
%
% Input Arguments:
%   str  : string to split
%   sep  : delimiter that indicates where in the string to split.
%          Can be a string or a cell array of strings
%   max  : max # of cuts to parse
% 
% Output Arguments:
%   out  : cell array of the parts contained within str
%
% Created:  07/15/07 by Ben Cipollini
% Last Mod: 03/15/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

allowed_escapes = {'\t','\r','\n'};

% Default separator: any whitespace char
if (~exist('sep','var'))
  sep = {' ' '\t' '\r' '\n'};
end;
if (~exist('max','var'))
  max_iter = Inf;
end;

out = {};

if iscell(sep) %  Allow multiple delimiters
  for i=1:length(sep)
    tmp   =  mmil_splitstr(str, sep{i});
    out   = {out{:} tmp{:}};
  end;
  out = unique(out);
else %  Do a single delimiter
  % Find the indices of the delimiter
  if ismember(sep,allowed_escapes)
    delimIdx = regexp(str,sep);
  else
    delimIdx = regexp(str,regexptranslate('escape',sep));
  end;

  % Find the start (1) and end (2) indices for each non-separator
  % substring
  strIdx(:,1) = [1 length(sep)+delimIdx(1:end)];
  strIdx(:,2) = [delimIdx(1:end)-1 length(str)];

  strIdx = strIdx(intersect(find(strIdx(:,2)>0),...
    find(strIdx(:,1)<=length(str))), :);

  out = cell( min(size(strIdx,1), max_iter), 1);
  for i=1:length(out)
    out{i} = str(strIdx(i,1):strIdx(i,2));
  end;
  out{end} = [out{end} str(strIdx(i,2)+1:end)];
end;

