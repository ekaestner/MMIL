function mmil_write_tags(fid,tags,parms)
%function mmil_write_tags(fid,tags,parms)
%
% Purpose: Help create matlab script for batch processing
%   Writes ('key',value) pairs to open file specified by fid
%   Handles multiple data types, including vectors and cell arrays
%
% Created:  05/17/09 by Don Hagler
% Last Mod: 02/25/11 by Don Hagler
%

for t=1:length(tags)
  if isfield(parms,tags{t})
    value = getfield(parms,tags{t});
  else
    value = [];
  end;
  fprintf(fid,'  ''%s'',',tags{t});
  mmil_write_value(fid,value);
  if t<length(tags)
    fprintf(fid,',...\n');
  end;
end;



