function [rois,info] = dti_read_roimap(fname)
%function [rois,info] = dti_read_roimap(fname)
%
% Purpose: read ROI map file from DTI Studio
% 
% Input:
%   fname: full or relative path of roi map file
%
% Output:
%   rois: structure defining ROIs to be used in fiber tracking
%     example:
%     rois(1).fname = 'CST_L_RoiMap_00.dat'
%     rois(1).logic = 'OR'
%     rois(2).fname = 'CST_L_RoiMap_01.dat'
%     rois(2).logic = 'AND'
%     rois(3).fname = 'CST_L_RoiMap_02.dat'
%     rois(3).logic = 'AND'
%     rois(4).fname = 'CST_L_RoiMap_03.dat'
%     rois(4).logic = 'NOT'
%   info: image information containing these fields:
%     nx,ny,nz (pixel num x, y, z)
%     dx,dy,dz (pixel size x, y, z)
%     sliceorientation: 0=Coronal, 1=Axial, 2=Sagittal
%     slicesequencing: 0=Positive, 1=Negative
%
% Created:  03/25/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

rois = [];
info = [];

if ~mmil_check_nargs(nargin,1), return; end;

% check that input files exist
if ~exist(fname,'file')
  error('file %s not found',fname);
end;

fid = fopen(fname);
begin_flag = 0;
nrois = 0;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break; end;
  if strcmp(tline,'Begin')
    begin_flag = 1;
    continue;
  end;
  if ~begin_flag, continue; end;
  if strcmp(tline,'End'), break; end;
  tline = strtrim(tline);
  k = findstr(tline,':');
  if isempty(k), continue; end;
  k=k(1);
  name = tline(1:k-1);
  val = strtrim(tline(k+1:end));
  switch name
    case 'ImageWidth'
      info.nx = str2double(val);
    case 'ImageHeight'
      info.ny = str2double(val);
    case 'ImageSlices'
      info.nz = str2double(val);
    case 'PixelSizeWidth'
      info.dx = str2double(val);
    case 'PixelSizeHeight'
      info.dy = str2double(val);
    case 'SliceThickness'
      info.dz = str2double(val);
    case 'SliceOrientation'
      info.sliceorientation = str2double(val);
    case 'SliceSequencing'
      info.slicesequencing = str2double(val);
    case 'BinaryFile'
      nrois = nrois + 1;
      j=findstr(val,'\');
      if ~isempty(j)
        val = val(j(end)+1:end);
      end;
      rois(nrois).fname = val;
    case 'Operation'
      val = str2double(val);
      switch val
        case 0
          rois(nrois).logic = 'OR';
        case 1
          rois(nrois).logic = 'AND';
        case 3
          rois(nrois).logic = 'NOT';
        case 4
          rois(nrois).logic = 'Cut0';
        case 5
          rois(nrois).logic = 'Cut1';
        case 6
          rois(nrois).logic = 'AndClr';
        case 7
          rois(nrois).logic = 'NotClr';
        case 8
          rois(nrois).logic = 'CutClr0';
        case 9
          rois(nrois).logic = 'CutClr1';
        otherwise
          rois(nrois).logic = [];
      end;
  end;
end;
fclose(fid);

