function dti_fiber_path_to_mask(fname_path,fname_mask,countflag,maskval)
%function dti_fiber_path_to_mask(fname_path,fname_mask,[countflag],[maskval])
%
% Purpose: convert DTI Studio fiber path file to mask file
%
% Required Input:
%   fname_path: full or relative path of DTI Studio format fiber path file
%   fname_mask: full or relative path of output fiber mask file
%
% Optional Input:
%   countflag: [0|1] toggle output of number of fibers for each
%     voxel (rather than binary mask)
%     if 1, will save as float dat file
%     if 0, will save as uint8 dat file
%     {default = 0}
%   maskval: value assigned to voxels inside fiber
%     {default = 127}
%
% Created:  08/06/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('countflag','var') | isempty(countflag), countflag=0; end;
if ~exist('maskval','var') | isempty(maskval), maskval=127; end;

excl_firstlast_flag = 1;

fiber = dti_read_DTIStudio_fiber(fname_path);
if isempty(fiber)
  fprintf('%s: WARNING: failed to read fiber path %s\n',mfilename,fname_path);
  return;
end;
nx = fiber.header.nImgWidth;
ny = fiber.header.nImgHeight;
nz = fiber.header.nImgSlices;
volmask = zeros(nx,ny,nz);
for i=1:length(fiber.chain)
  npoints = fiber.chain(i).nLength;
  SelStatus = fiber.chain(i).nSelStatus;
  % exclude first and last points along chain (and those excluded by DTI Studio)
  if excl_firstlast_flag
    FirstPoint = max(2,fiber.chain(i).nSelBeginIdx);
    LastPoint = min(npoints-1,fiber.chain(i).nSelEndIdx);
  else
    FirstPoint = max(1,fiber.chain(i).nSelBeginIdx);
    LastPoint = min(npoints,fiber.chain(i).nSelEndIdx);
  end;
  pxyzChain = fiber.chain(i).pxyzChain(FirstPoint:LastPoint,:);
  pxyzChain(:,1) = min(max(pxyzChain(:,1),1),nx);
  pxyzChain(:,2) = min(max(pxyzChain(:,2),1),ny);
  pxyzChain(:,3) = min(max(pxyzChain(:,3),1),nz);
  x = round(pxyzChain(:,1)+0.5);
  y = round(pxyzChain(:,2)+0.5);
  z = round(pxyzChain(:,3)+0.5);
  ind = sub2ind(size(volmask),x,y,z);
  if countflag
    for j=1:length(ind)
      k = ind(j);
      volmask(k) = volmask(k)+1;
    end;
  else
    volmask(ind) = maskval;
  end;
end;
if countflag
  dti_save_dat(volmask,fname_mask,[],[],'float');
else
  dti_save_dat(volmask,fname_mask,[],[],'uint8');
end;



