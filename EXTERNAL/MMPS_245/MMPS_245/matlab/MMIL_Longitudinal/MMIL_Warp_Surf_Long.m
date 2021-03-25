function errcode = MMIL_Warp_Surf_Long(FSRootDir,FSDir_source,FSDir_targ,nlreg_dir,forceflag);
%function errcode = MMIL_Warp_Surf_Long(FSRootDir,FSDir_source,FSDir_targ,nlreg_dir,[forceflag]);
%
% Required Input:
%   FSRootDir: Freesurfer subject root directory
%   FSDir_source: Source Freesurfer directory
%   FSDir_targ: Target Freesurfer directory
%   nlreg_dir: directory containing nonlinear registration volumes
%
% Optional Input:
%   forceflag: [0|1] overwrite exiting output files
%     {default = 0}
%
% Created:  06/21/07 by Don Hagler
% Rcnt Mod: 11/19/10 by Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

errcode = 0;
if nargin<4, help(mfilename); return; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag = 0; end
hemilist = {'lh','rh'};

if forceflag
  MMIL_CleanUp_Longitudinal_Warped_Surface(FSRootDir,FSDir_targ);
end;

% warp source surfaces to target
for h=1:length(hemilist)
  hemi = hemilist{h};
  fname_surf_in = sprintf('%s/%s/surf/%s.white',FSRootDir,FSDir_source,hemi);
  fname_surf_out = sprintf('%s/%s/surf/%s.orig',FSRootDir,FSDir_targ,hemi);
  fname_dx = sprintf('%s/dx.mgz',nlreg_dir);
  fname_dy = sprintf('%s/dy.mgz',nlreg_dir);
  fname_dz = sprintf('%s/dz.mgz',nlreg_dir);
  mmil_warp_freesurfer_surface(fname_surf_in,fname_surf_out,...
    fname_dx,fname_dy,fname_dz,forceflag);
  if ~exist(fname_surf_out,'file')
    fprintf('%s: ERROR: failed to warp surface to target for %s to %s\n',mfilename,FSDir_source,FSDir_targ);
    errcode = 1;
    return;
  end;
  % copy sphere, sphere.reg, aparc
  cmdlist = {...
    sprintf('cp %s/%s/surf/%s.sphere %s/%s/surf/%s.sphere',...
      FSRootDir,FSDir_source,hemi,FSRootDir,FSDir_targ,hemi),...
    sprintf('cp %s/%s/surf/%s.sphere.reg %s/%s/surf/%s.sphere.reg',...
      FSRootDir,FSDir_source,hemi,FSRootDir,FSDir_targ,hemi),...
    sprintf('cp %s/%s/label/%s.aparc.annot %s/%s/label/%s.aparc.annot',...
      FSRootDir,FSDir_source,hemi,FSRootDir,FSDir_targ,hemi),...
    sprintf('cp %s/%s/label/%s.aparc.a2005s.annot %s/%s/label/%s.aparc.a2005s.annot',...
      FSRootDir,FSDir_source,hemi,FSRootDir,FSDir_targ,hemi)...
  };
  for c=1:length(cmdlist)
    cmd = cmdlist{c};
    [status,result] = unix(cmd);
    if status
      fprintf('%s: ERROR: cmd %s failed:\n',mfilename,cmd);
      disp(result);
      return;
    end;
  end;
end;        

% prepare for final surfs recon (only) for target
%fs_mktouchfiles(FSDir_targ,'subjdir',FSRootDir,'what','longitudinalwarp');


