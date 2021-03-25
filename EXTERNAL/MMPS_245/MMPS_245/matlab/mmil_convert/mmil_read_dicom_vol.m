function [vol,M] = mmil_read_dicom_vol(fnames)
%function [vol,M] = mmil_read_dicom_vol(fnames)
%
% Required Input:
%   fnames: cell array of dicom files for 3D volume
%
% Output:
%   vol: 3D volume matrix
%   M: vox2ras matrix
%
% Created:  01/01/05 by Anders Dale
% Last Mod: 08/06/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
vol = []; M = eye(4);

warnstatussave = warning;
warning('off');

if ~iscell(fnames), fnames = {fnames}; end;
nfiles = length(fnames);

nii_flag = 1; %1 = convert to nifti using dcm2nii and read vol. 0 = read dicoms directly using dicomread

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = mmil_read_dicom_M(fnames);
if isempty(M)
  vol=[];
  return;
end;

if nii_flag
  vol = mmil_dcm2nii2vol(fnames);
else
  dcminfo = dicominfo(char(fnames{1}));
  nrow = dcminfo.Rows;
  ncol = dcminfo.Columns;
  nslices = nfiles;
  vol = zeros(ncol,nrow,nslices);
  for n=1:nfiles
    fname = char(fnames{n});
    if ~exist(fname,'file')
      fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
      vol=[];
      return;
    end;
    try
      im = dicomread(fname)';
    catch
      fprintf('%s: ERROR: dicom file %s may be corrupt\n',mfilename,fname);
      vol=[];
      return;
    end;
    if(isempty(im))
      fprintf('%s: ERROR: could not load pixel data from %s\n',mfilename,fname);
      vol=[];
      return;
    end
    try
      vol(:,:,n) = im;
    catch
      fprintf('%s: ERROR: problem with pixel data from %s: %s\n',...
        mfilename,fname,lasterr);
      vol=[];
      return;
    end;
  end
end
warning(warnstatussave);

