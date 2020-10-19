function dti_write_DTIStudio_fiber(fiber,fname)
%function dti_write_DTIStudio_fiber(fiber,fname)
%
% Purpose: write DTI Studio format streamline fiber file
%
% Required Input:
%   fiber: fiber struct containing fields
%     header: struct with fields:
%       sFiberFileTag
%     	nFiberNr
%     	nFiberLenMax
%     	fFiberLenMean
%     	nImgWidth
%     	nImgHeight
%     	nImgSlices
%     	fPixelSizeWidth
%     	fPixelSizeHeight
%     	fSliceThickness
%     	enumSliceOrientation
%     	enumSliceSequencing
%     	offset
%     chain: struct array with fields:
%       nLength
%   		nSelStatus
%   		rgbFiberClr
%   		nSelBeginIdx
%   		nSelEndIdx
%   		pxyzChain
%   fname: name of DTI Studio fiber file
%
% See also: dti_read_DTIStudio_fiber
%
% Created:  05/11/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% open file
fid = fopen(fname,'wb');
if fid==-1
  error('failed to open file %s for writing',fname);
end;

% write header information
fwrite(fid,fiber.header.sFiberFileTag,'schar');
fwrite(fid,fiber.header.nFiberNr,'int');
fwrite(fid,fiber.header.nFiberLenMax,'int');
fwrite(fid,fiber.header.fFiberLenMean,'float');
fwrite(fid,fiber.header.nImgWidth,'int');
fwrite(fid,fiber.header.nImgHeight,'int');
fwrite(fid,fiber.header.nImgSlices,'int');
fwrite(fid,fiber.header.fPixelSizeWidth,'float');
fwrite(fid,fiber.header.fPixelSizeHeight,'float');
fwrite(fid,fiber.header.fSliceThickness,'float');
fwrite(fid,fiber.header.enumSliceOrientation,'int');
fwrite(fid,fiber.header.enumSliceSequencing,'int');

% write zeros to get to start of fiber information
offset = ftell(fid);
offset_bytes = 128-offset;
if offset_bytes<0
  error('header is longer than 128 bytes');
elseif offset_bytes>0
  fwrite(fid,zeros(1,offset_bytes,'int8'));
end;

% write fiber chain coordinates
for i = 1:fiber.header.nFiberNr
  fwrite(fid,fiber.chain(i).nLength,'int');
  fwrite(fid,fiber.chain(i).nSelStatus,'uchar');
  fwrite(fid,fiber.chain(i).rgbFiberClr,'uchar');
  fwrite(fid,fiber.chain(i).nSelBeginIdx,'int');
  fwrite(fid,fiber.chain(i).nSelEndIdx,'int');
  pxyzChain = fiber.chain(i).pxyzChain;
  nLength = fiber.chain(i).nLength;
  pxyzChain = single(reshape(pxyzChain',[1,3*nLength]));
  fwrite(fid,pxyzChain,'float');
end

fclose(fid);

