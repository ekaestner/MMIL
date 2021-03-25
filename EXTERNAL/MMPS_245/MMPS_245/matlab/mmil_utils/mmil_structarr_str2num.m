function infostruct=mmil_structarr_str2num(infostruct,tags)
%function infostruct=mmil_structarr_str2num(infostruct,[tags])
%
% Purpose: for specified fields in struct array, convert
%   from string to numeric vector if possible
%
% Required Input:
%   infostruct: struct array
%
% Optional Input:
%   tags: cell array of field names to convert from string to numeric
%     {default = all fields}
%
% Created:  02/17/11 by Don Hagler
% Last Mod: 03/23/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('tags','var'), tags = []; end;

if ~isstruct(infostruct)
  error('input infostruct is not a struct');
end;

if isempty(tags)
  tags = fieldnames(infostruct);
else
  tags = intersect(tags,fieldnames(infostruct));
end;

for i=1:length(infostruct)
  for t=1:length(tags)
    tag = tags{t};
    val = infostruct(i).(tag);
    if isempty(val), continue; end;
    if ~isnumeric(val)
      if ischar(val)
        val = str2num(val);
        if ~isempty(val) && mmil_isint(val), val = int32(val); end;
      else
        val = [];
      end;
    elseif mmil_isint(val)
      val = int32(val);
    end;
    if isempty(val)
      fprintf('%s: WARNING: failed to convert field %s (%s) to numeric vector\n',...
        mfilename,tag,infostruct(i).(tag));
    else
      infostruct(i).(tag) = val;
    end;
  end;
end;

