function [slice_tpattern,slice_order] = mmil_slicetiming_Siemens(SeriesInfo) 
%function [slice_tpattern,slice_order] = mmil_slicetiming_Siemens(SeriesInfo)
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
% Last Mod: 10/29/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sources of information:

% http://www.nmr.mgh.harvard.edu/~greve/dicom-unpack

% https://wiki.cimec.unitn.it/tiki-index.php?page=MRIBOLDfMRI&offset=&sort_mode=hits_desc&atts_show=y#Slice_acquisition_order

% https://www.icts.uiowa.edu/confluence/plugins/viewsource/viewpagesrc.action?pageId=54756326

% http://practicalfmri.blogspot.com/2012/07/siemens-slice-ordering.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
slice_tpattern = 'unknown'; slice_order = [];

[mdhStruct,txt] = mmil_read_mdhStruct(SeriesInfo.info);
if isfield(mdhStruct,'sSliceArrayducMode')
  slice_mode = get_mdh_tagval(mdhStruct,'sSliceArrayducMode');
else
  fprintf('%s: WARNING: missing ucMode field required to determine slice order\n',...
    mfilename);
  return;
end;
if isfield(mdhStruct,'sSliceArraydasSlicel0rdsNormalddTra')
  slice_dir = ...
    get_mdh_tagval(mdhStruct,'sSliceArraydasSlicel0rdsNormalddTra');
else
  fprintf('%s: WARNING: missing sNormal.Tra field required to determine slice order\n',...
    mfilename);
  return;
end;
if isfield(mdhStruct,'sSliceArrayducImageNumbTra')
  slice_rev = get_mdh_tagval(mdhStruct,'sSliceArrayducImageNumbTra');
  if slice_rev~=0
    slice_dir = -slice_dir;
  end;
end;
slice_pos = (slice_dir > 0);
if mod(SeriesInfo.nslices,2)
  slice_odd = 1;
else
  slice_odd = 0;
end;

switch slice_mode
  case '0x1'
    if slice_pos
      slice_tpattern = 'seq+z';
      slice_order = [1:SeriesInfo.nslices];
    else
      slice_tpattern = 'seq-z';
      slice_order = [SeriesInfo.nslices:-1:1];
    end;
  case '0x2'
    if slice_pos
      slice_tpattern = 'seq-z';
      slice_order = [SeriesInfo.nslices:-1:1];
    else
      slice_tpattern = 'seq+z';
      slice_order = [1:SeriesInfo.nslices];
    end;
  case '0x4'
    if slice_pos
      if slice_odd
        slice_tpattern = 'alt+z';
        slice_order = [1:2:SeriesInfo.nslices,2:2:SeriesInfo.nslices-1];
      else
        slice_tpattern = 'alt+z2';
        slice_order = [2:2:SeriesInfo.nslices,1:2:SeriesInfo.nslices-1];
      end;
    else
      if slice_odd
        slice_tpattern = 'alt-z';
        slice_order = [SeriesInfo.nslices:-2:1,SeriesInfo.nslices-1:-2:2];
      else
        slice_tpattern = 'alt-z2';
        slice_order = [SeriesInfo.nslices:-2:1,SeriesInfo.nslices-1:-2:2];
      end;
    end;
end;
slice_order = double(slice_order);

return;

