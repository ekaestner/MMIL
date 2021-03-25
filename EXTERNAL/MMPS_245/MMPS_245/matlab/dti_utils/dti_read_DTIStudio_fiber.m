function fiber = dti_read_DTIStudio_fiber(fname)
%function fiber = dti_read_DTIStudio_fiber(fname)
%
% Purpose: read DTI Studio format streamline fiber file
%
% Required Input:
%   fname: name of DTI Studio fiber file
%
% Output:
%   fiber: struct containing fiber streamline information
%
% Created:  08/06/07 by Sumiko Abe
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
fiber = [];

% open file
fid = fopen(fname, 'r');
if(fid == -1)
  error('failed to open fiber file %s',fname);
end

% initialize header
header = struct(...
  'sFiberFileTag',        '        ', ...
  'nFiberNr',             -1,  ...
  'nFiberLenMax',         -1, ...
  'fFiberLenMean',        0.0,  ...
  'nImgWidth',            -1, ...
  'nImgHeight',           -1,...
  'nImgSlices',           -1,...
  'fPixelSizeWidth',      0.0, ...
  'fPixelSizeHeight',     0.0, ...
  'fSliceThickness',      0.0, ...
  'enumSliceOrientation', 0.0, ...
  'enumSliceSequencing',  0.0, ...
  'offset',               -1 ...
);

% read header information
header.sFiberFileTag        = fread(fid, 8, 'schar');
header.nFiberNr             = fread(fid, 1, 'int'); 
header.nFiberLenMax         = fread(fid, 1, 'int'); 
header.fFiberLenMean        = fread(fid, 1, 'float');
header.nImgWidth            = fread(fid, 1, 'int');
header.nImgHeight           = fread(fid, 1, 'int');
header.nImgSlices           = fread(fid, 1, 'int');
header.fPixelSizeWidth      = fread(fid, 1, 'float');
header.fPixelSizeHeight     = fread(fid, 1, 'float');
header.fSliceThickness      = fread(fid, 1, 'float');
header.enumSliceOrientation = fread(fid, 1, 'int');
header.enumSliceSequencing  = fread(fid, 1, 'int');
header.offset = ftell(fid);

% check for errors in reading
if(header.offset == -1)
  fprintf('%s: WARNING: internal error with fiber file %s\n',mfilename,fname);
  return;
end;
if header.nFiberNr<=0
  fprintf('%s: WARNING: zero fibers in file %s\n',mfilename,fname);
  return;
end;

% initialize chain structure
clear chain;
for i = 1:header.nFiberNr
  chain(i) = struct(...
    'nLength',       '-1', ...
    'nSelStatus',   '-1', ...
    'rgbFiberClr',  zeros(3, 1), ...
    'nSelBeginIdx', '-1', ...
    'nSelEndIdx',   '-1', ...
    'pxyzChain',    zeros(0,3,'single'));
end

% skip to start of fiber information
fseek(fid, 128, 'bof');

% read coordinate information
for i = 1:header.nFiberNr
  nLength       = fread(fid,1,'int');
  nSelStatus    = fread(fid,1,'uchar');
  rgbFiberClr   = fread(fid,3,'uchar');
  nSelBeginIdx  = fread(fid,1,'int');
  nSelEndIdx    = fread(fid,1,'int');
  chain(i).nLength   = nLength;
  chain(i).nSelStatus = nSelStatus;
  chain(i).rgbFiberClr(:, 1)  = rgbFiberClr;
  chain(i).nSelBeginIdx = nSelBeginIdx;
  chain(i).nSelEndIdx = nSelEndIdx;
  if nLength > header.nFiberLenMax
    fprintf('%s: WARNING: nLength = %d\n',mfilename,nLength);
    return;
  end;
  pxyzChain = fread(fid, nLength*3, 'float');
  pxyzChain = reshape(pxyzChain,[3,nLength])';
  chain(i).pxyzChain = single(pxyzChain);
  clear nLength nSelStatus rgbFiberClr nSelBeginIdx nSelEndIdx pxyzChain;
end

% close file
fclose(fid);

fiber.header = header;
fiber.chain = chain;

clear header chain;

return;

