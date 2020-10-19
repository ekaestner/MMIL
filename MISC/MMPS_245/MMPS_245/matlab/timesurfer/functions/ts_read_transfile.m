function T=ts_read_transfile(fname,translation_scale)
%function T=ts_read_transfile(fname,[translation_scale])
%
% required input:
%   fname - name of text file containing 4x4 transformation matrix
%
% optional parameters:
%   translation_scale - scaling factor applied to translations
%     purpose is to have resultant matrix with units of meters
%     for old .trans files created with freesurfer's tkmedit, use 1/1000
%       to convert from millimeters to meters
%     {default: 1}
%
% output:
%   T - 4x4 transformation matrix, units of translations should
%       be in meters (assuming translation_scale is correct)
%

if nargin<1
  help(mfilename);
  return;
end;

if ~exist('translation_scale','var')
  translation_scale = 1;
end;

T = eye(4);
fid=fopen(fname,'rt');
j = 1;
while (~feof(fid))
  temp=fgetl(fid);
  A = sscanf(temp,'%f');
  if length(A)>4, break; end;
  T(j,:)=A';
  j = j+1;
  if j>4, break; end;
end
fclose(fid);

if translation_scale~=1
  T(1,4) = T(1,4)*translation_scale;
  T(2,4) = T(2,4)*translation_scale;
  T(3,4) = T(3,4)*translation_scale;
end;

return;
