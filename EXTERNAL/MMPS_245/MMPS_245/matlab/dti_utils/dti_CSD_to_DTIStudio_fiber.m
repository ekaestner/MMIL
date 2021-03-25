function dti_CSD_to_DTIStudio_fiber(fname_in,fname_out,varargin)
%function dti_CSD_to_DTIStudio_fiber(fname_in,fname_out,[options])
%
% Purpose: write DTI Studio format streamline fiber file
%
% Required Input:
%   fname_in: name of CSD tractography mat file
%   fname_out: name of DTI Studio fiber file
%
% Optional Input:
%   'imgorient': [0|1|2] slice orientation
%     0=Coronal, 1=Axial, 2=Sagittal
%     {default = 1}
%   'imgseq': [0|1] slice sequencing
%     0=Positive, 1=Negative
%     {default = 0}
%   'M': vox2ras matrix of data
%     if supplied, will set imgorient and imgseq accordingly
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% See also: dti_CSD_tracto
%           dti_write_DTIStudio_fiber
%           dti_read_DTIStudio_fiber
%
% Created:  05/12/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'imgorient',1,[0:2],...
  'imgseq',0,[0,1,1],...
  'M',[],[],...
  'forceflag',false,[false true],...
});

if exist(fname_out,'file') && ~parms.forceflag, return; end;
if ~exist(fname_in,'file'), error('file %s not found',fname_in); end;

if ~isempty(parms.M)
  orient = fs_read_orient([],parms.M);
  switch orient(3)
    case {'L','R'}
      parms.imgorient = 2;
    case {'P','A'}
      parms.imgorient = 0;
    case {'I','S'}
      parms.imgorient = 1;
  end;
  switch orient(3)
    case {'L','P','I'}
      parms.imgseq = 1;
    otherwise
      parms.imgseq = 0;
  end;
end;

load(fname_in);

if ~exist('Tracts','var') ||...
   ~exist('TractMask','var') ||...
   ~exist('VDims','var')
  error('input file %s is missing required data',fname_in);
end;

fiber_lengths = cellfun(@length,Tracts);

% header information
header = [];
header.sFiberFileTag        = double('Fiberdat')';
header.nFiberNr             = length(fiber_lengths);
header.nFiberLenMax         = max(fiber_lengths);
header.fFiberLenMean        = mean(fiber_lengths);
header.nImgWidth            = size(TractMask,1);
header.nImgHeight	          = size(TractMask,2);
header.nImgSlices           = size(TractMask,3);
header.fPixelSizeWidth      = VDims(1);
header.fPixelSizeHeight     = VDims(2);
header.fSliceThickness      = VDims(3);
header.enumSliceOrientation = parms.imgorient;
header.enumSliceSequencing  = parms.imgseq;

% initialize chain structure
clear chain;
for i = 1:header.nFiberNr
	chain(i) = struct(...
    'nLength', 	    '-1', ...
		'nSelStatus',   '-1', ...
		'rgbFiberClr',  zeros(3,1), ...
		'nSelBeginIdx', '-1', ...
		'nSelEndIdx',   '-1', ...
		'pxyzChain',    zeros(0,3,'single'));
end

% streamline information
for i = 1:header.nFiberNr
  nLength = size(Tracts{i},1);
  pxyzChain = Tracts{i} ./ repmat(VDims,nLength,1) - 0.5;
  chain(i).nLength = nLength;
  chain(i).nSelStatus    = 1;
  chain(i).rgbFiberClr   = randi(255,3,1);
  chain(i).nSelBeginIdx  = 0;
  chain(i).nSelEndIdx    = nLength-1;;
  chain(i).pxyzChain = single(pxyzChain);
  clear nLength pxyzChain;
end

% add header and chain to fiber structure
fiber.chain = chain;
fiber.header = header;

clear chain header;

% write DTI Studio fiber file
dti_write_DTIStudio_fiber(fiber,fname_out);

