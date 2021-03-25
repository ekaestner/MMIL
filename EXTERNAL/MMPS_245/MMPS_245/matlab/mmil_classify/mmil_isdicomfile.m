function res = mmil_isdicomfile(fname)
%function res = mmil_isdicomfile(fname)
%
% Created:  03/14/11 by Don Hagler (based on code from Anders Dale)
% Last Mod: 03/14/11 by Don Hagler
%

if strcmp(fname(end-7:end),'DICOMDIR')
  res = 0;
  return;
end;

if ((length(fname)>4)&strcmp(lower(fname(end-3:end)),'.dcm'))
  res = 1;
  return
end
fid = fopen(fname,'r');
if fid < 0
  res = 0;
else
  stat = fseek(fid,128,'bof'); % move to DICM string
  tmp = char(fread(fid,4,'uchar')');
  res = strcmp(tmp,'DICM');
  fclose(fid);
end

