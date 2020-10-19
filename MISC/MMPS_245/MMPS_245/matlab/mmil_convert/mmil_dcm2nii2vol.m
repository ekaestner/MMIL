function vol = mmil_dcm2nii2vol(fnames)
% function vol = mmil_dcm2nii2vol(fnames)
% 
% Purpose: convert dcm to nii then read vol data
%
% Required:
%    fnames: cell array of dcm filenames, will use it to extract root dir of dcm files.
%
% Output:
%    vol: volume pixal data
%
% Created:         09/05/2017 by Feng Xue
% Last Mod:        09/06/2017 by Feng Xue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

if ~iscell(fnames)
  %error('fnames should be a cell array');
  fprintf('ERROR: fnames should be a cell array\n');
  return
end
[dcmdir,~,~] = fileparts(char(fnames{1}));
vol = [];

%if exist(dcmdir,'dir') ~= 7, error('dcm dir does not exist'); end;
if exist(dcmdir,'dir') ~= 7
  fprintf('ERROR: dcm dir does not exist\n');
  return;
end;
tempfolder = mmil_tempfname;
%if exist(tempfolder,'dir') == 7, error('temp dir exists, pls check'); end;
if exist(tempfolder,'dir') == 7
  fprintf('ERROR: temp dir exists, pls check\n');
  return;
end;
mkdir(tempfolder);
%TODO: use -c option to save mr_parms
cmd=sprintf('dcm2niix_noflipY -z n -p y -f %s -o %s %s','%n_%k_%j_%t_%s',tempfolder,dcmdir);
unix(cmd);
niifile=dir(sprintf('%s/*.nii',tempfolder));

if length(niifile) ~=1
  fprintf('multiple nifti files found in %s, pls check\n',tempfolder);
  return
end;
nifti = fullfile(tempfolder,niifile.name);
[vol,M] = fs_load_nifti(nifti);
orientation_nii = fs_read_orient('',M);

if isempty(vol), return; end;

M_dcm=mmil_read_dicom_M(char(fnames{1}));
orientation_dcm = fs_read_orient('',M_dcm);

if isempty(M_dcm)
  fprintf('Check dcmdir\n');
  return;
end;

%TODO: test DmyDisableReorient3dToOrtho switch and comment below if works.
if ~strcmp(orientation_nii,orientation_dcm),
  [vol,~] = fs_reorient(vol,M,orientation_dcm);
end
rmdir(tempfolder,'s');
