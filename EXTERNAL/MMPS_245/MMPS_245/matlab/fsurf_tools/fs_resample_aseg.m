function fs_resample_aseg(fname_in,fname_reg,fname_ref,fname_out,forceflag)
%function fs_resample_aseg(fname_in,fname_reg,fname_ref,fname_out,forceflag)
%
% Purpose: Resample aseg volume to space of functional volume
%
% Required Input:
%   fname_aseg: full path to input segmentation volume (e.g. aseg.mgz)
%   fname_reg: register.dat file containing 4x4 registration matrix
%     from fname_ref to FreeSurfer T1 (same space as aseg)
%   fname_ref: full path of reference volume (e.g. fMRI)
%   fname_out: full path to output segmentation volume (e.g. aseg_resBOLD.mgz)
%
% Optional Input:
%   forceflag: [0|1] toggle overwrite of existing output file
%     {default = 0}
%
% Created:  07/12/12 by Don Hagler
% Last Mod: 07/12/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,4), return; end;

if ~exist('forceflag','var') | isempty(forceflag), forceflag=0; end;

if ~exist(fname_in,'file')
  error('file %s not found',fname_aseg);
end;
if ~exist(fname_reg,'file')
  error('file %s not found',fname_reg);
end;
if exist(fname_out,'file') & ~forceflag, return; end;

% resample aseg to BOLD resolution
fprintf('%s: resampling %s to resolution of %s...\n',...
  mfilename,fname_in,fname_ref);
cmd = 'mri_vol2vol';
cmd = sprintf('%s --mov %s',cmd,fname_ref);
cmd = sprintf('%s --targ %s',cmd,fname_in);
cmd = sprintf('%s --o %s',cmd,fname_out);
cmd = sprintf('%s --reg %s',cmd,fname_reg);
cmd = sprintf('%s --inv --interp nearest',cmd);
[s,r] = mmil_unix(cmd);
if s
  error('failed to resample %s to resolution of %s:\n%s',...
    fname_in,fname_ref,r);
end;
fname_tmp = sprintf('%s.reg',fname_out);
if exist(fname_tmp,'file'), delete(fname_tmp); end;

