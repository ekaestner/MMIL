function [slice_tpattern,slice_order] = mmil_slicetiming_GE(SeriesInfo) 
%function [slice_tpattern,slice_order] = mmil_slicetiming_GE(SeriesInfo)
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
%   slice_order: vector of slice numbers in order of acquisition
%
% NOTE: for GE, hard-coded to be alt+z 
%
% Created:  10/19/12 by Don Hagler
% Last Mod: 10/20/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
slice_tpattern = []; slice_order = [];

slice_tpattern = 'alt+z';
slice_order = double([1:2:SeriesInfo.nslices,2:2:SeriesInfo.nslices]);

return;

