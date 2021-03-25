function M = mmil_read_dicom_M(fnames);
%function M = mmil_read_dicom_M(fnames);
%
% Required Input:
%   fnames: cell array of dicom file names
%
% Output:
%   M: vox2ras matrix
%
% Created:  01/01/05 by Anders Dale
% Last Mod: 12/22/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
M = eye(4);

if ~iscell(fnames), fnames = {fnames}; end;
nfiles = length(fnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname=char(fnames{1});
if ~exist(fname,'file')
  %% todo: error?
  fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
  M=[];
  return;
end;
dcminfo = dicominfo(fname);

dcr = dcminfo.ImageOrientationPatient(1:3);
dcc = dcminfo.ImageOrientationPatient(4:6);
impos = dcminfo.ImagePositionPatient;

M(1:3,1:2) = [dcminfo.PixelSpacing(1)*dcr...
              dcminfo.PixelSpacing(2)*dcc];

if length(fnames)>1
  fname=char(fnames{end});
  if ~exist(fname,'file')
    %% todo: error?
    fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
    M=[];
    return;
  end;
  dcminfo_end = dicominfo(fname);
  M(1:3,3) = (dcminfo_end.ImagePositionPatient-dcminfo.ImagePositionPatient)/(nfiles-1);
else
  dcs = cross(dcr,dcc);
  M(1:3,3) = dcminfo.SliceThickness*dcs;
end;

M(1:3,4) = impos-M(1:3,:)*[1 1 1 1]'; % Adjust for Matlab 1-based indexing

M = M_LPH_TO_RAS*M; % Convert from DICOM LPH to RAS coordinates

