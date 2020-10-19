function [slice_tpattern,slice_order] = mmil_slicetiming_Philips(SeriesInfo) 
%function [slice_tpattern,slice_order] = mmil_slicetiming_Philips(SeriesInfo)
%
% Purpose: gather slice timing information from SeriesInfo 
%
% Required Input:
%   SeriesInfo: struct containing information about dicom
%     initialized by mmil_classify_dicom
%
% Output:
%   slice_tpattern: slice timing pattern string, as used by AFNI's 3dTshift
%     e.g. 'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z'
%     returns 'unknown' if required dicom header fields are missing
%   slice_order: vector of slice numbers in order of acquisition
%
% Created:  10/19/12 by Don Hagler
% Prev Mod: 10/29/12 by Don Hagler
% Last Mod: 11/01/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
slice_tpattern = 'unknown'; slice_order = [];

%slice_tpattern = 'seq+z';
%slice_order = double([1:SeriesInfo.nslices]);

if length(SeriesInfo.FileNames) < SeriesInfo.nslices
  fprintf('%s: WARNING: number of files (%d) is less than number of slices (%d)\n',...
    mfilename,length(SeriesInfo.FileNames),SeriesInfo.nslices);
  return;
end;

slice_times = zeros(1,SeriesInfo.nslices);
for i=1:SeriesInfo.nslices
  tmpinfo = dicominfo(SeriesInfo.FileNames{i});
  UID = mmil_getfield(tmpinfo,'SOPInstanceUID');
  if isempty(UID)
    fprintf('%s: WARNING: missing SOPInstanceUID field required for determining slice order\n',...
      mfilename);
    return;
  end;
  n = regexp(UID,'\.(?<year>\d{4})(?<month>\d{2})(?<day>\d{2})(?<time>\d{11}$)','names');
  if isempty(n)
    fprintf('%s: WARNING: SOPInstanceUID does not have expected pattern\n',...
      mfilename);
    return;
  else
    slice_times(i) = str2double(n.time);
  end;
end;
[tmp,slice_order] = sort(slice_times);

diff_vals = diff(slice_order);
if all(diff_vals==1)
  slice_tpattern = 'seq+z';
elseif all(diff_vals==-1)
  slice_tpattern = 'seq-z';
elseif length(find(diff_vals==2))==SeriesInfo.nslices-2
  if slice_order(1) == 1
    slice_tpattern = 'alt+z';
  else
    slice_tpattern = 'alt+z2';  
  end;
elseif length(find(diff_vals==-2))==SeriesInfo.nslices-2
  if slice_order(1) == 1
    slice_tpattern = 'alt-z';
  else
    slice_tpattern = 'alt-z2';  
  end;
elseif length(diff_vals==1)>SeriesInfo.nslices-4
  % allow for some differences in timing (for unknown reasons)
  fprintf('%s: WARNING: unexpected slice order, but assuming seq+z\n',mfilename);
  slice_tpattern = 'seq+z';
else
  fprintf('%s: WARNING: unexpected slice order\n',mfilename);
end;

return;

