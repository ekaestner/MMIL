function ts_extract_trans_from_fiff(alignment_fiff,outtransfile,translation_scale)
%function ts_extract_trans_from_fiff(alignment_fiff,outtransfile,[translation_scale])
%
% required input:
%   alignment_fiff - fif file containing 4x4 mri2head transformation
%   outtransfile - name of text file to write out 4x4 matrix
%
% optional parameters:
%   translation_scale - scaling factor to make translations in meters
%     assuming T has units of meters, to write with same convention as
%     old .trans files created with tkmedit, use 1000 to convert from
%     meters to millimeters
%     {default: 1}
%

if ~exist('translation_scale','var')
  translation_scale = 1;
end;

if ~exist(alignment_fiff,'file')
  fprintf('%s: alignment_fiff %s not found...quitting\n',...
    funcname,alignment_fiff);
  return;
end

T=loadtrans(alignment_fiff,'MRI','HEAD');

if translation_scale~=1
  T(1,4) = T(1,4)*translation_scale;
  T(2,4) = T(2,4)*translation_scale;
  T(3,4) = T(3,4)*translation_scale;
end;

fid=fopen(outtransfile,'wt');
for j=1:4
  fprintf(fid,'%0.6f %0.6f %0.6f %0.6f\n',...
    T(j,1),T(j,2),T(j,3),T(j,4));
end;
fclose(fid);

return;
