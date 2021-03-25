function fs_save_nifti(vol,fname,M,mr_parms,orient,options)
%function fs_save_nifti(vol,fname,M,mr_parms,orient,options)
%
% fs_save_nifti(vol,fname, M, <mr_parms>);
%
% M is the 4x4 vox2ras transform such that
% y(i1,i2,i3), xyz = M*[i1 i2 i3 1] where the
% indicies are 1-based. M is converted to
% 0-based before saving.
%
% mr_parms = [tr flipangle te ti fov]
% orient: orientation string, e.g. 'RAS'
%
% See also: fs_save_mgh, fs_mri_convert
%
% Created:  11/16/17 by Don Hagler
% Last Mod: 11/17/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('mr_parms','var'), mr_parms = []; end;
if ~exist('orient','var'), orient = []; end;
if ~exist('options','var'), options = []; end;

[fpath,fstem,fext] = fileparts(fname);
if isempty(fpath)
  fpath = pwd;
elseif ~exist(fpath,'dir')
  mmil_mkdir(fpath);
end;
fname_tmp = sprintf('%s.mgh',mmil_tempfname(fstem,fpath));

fs_save_mgh(vol,fname_tmp,M,mr_parms);

fs_mri_convert(fname_tmp,fname,...
  'out_orient',orient,...
  'options',options,...
  'forceflag',1);

delete(fname_tmp);

