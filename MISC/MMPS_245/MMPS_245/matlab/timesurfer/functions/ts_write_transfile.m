function ts_write_transfile(fname,T,translation_scale)
%function ts_write_transfile(fname,T,[translation_scale])
%
% required input:
%   fname - name of text file containing 4x4 transformation matrix
%   T - 4x4 transformation matrix, units of translations should
%       be in meters
%
% optional parameters:
%   translation_scale - scaling factor to make translations in meters
%     assuming T has units of meters, to write with same convention as
%     old .trans files created with tkmedit, use 1000 to convert from
%     meters to millimeters
%     {default: 1}
%

if nargin < 2
  help(mfilename);
  return;
end;

if ~exist('translation_scale','var')
  translation_scale = 1;
end;

if translation_scale~=1
  T(1,4) = T(1,4)*translation_scale;
  T(2,4) = T(2,4)*translation_scale;
  T(3,4) = T(3,4)*translation_scale;
end;

fid=fopen(fname,'wt');
for j=1:4
  fprintf(fid,'%0.6f %0.6f %0.6f %0.6f\n',...
    T(j,1),T(j,2),T(j,3),T(j,4));
end;
fclose(fid);

return;
