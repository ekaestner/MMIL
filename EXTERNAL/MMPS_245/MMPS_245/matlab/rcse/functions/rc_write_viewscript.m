function rc_write_viewscript(subj,varargin)
%function rc_write_viewscript(subj,[options])
%
% Purpose: create csh script to view retinotopy data with tksurfer
%
% Required Parameters:
%   subj: FreeSurfer recon subject name
%
% Optional Parameters:
%   'fname_out': output file name
%     {default = 'view_ret.csh'}
%   'fstem': retinotopy data file stem
%     {default = 'pol'}
%   'indir': input directory
%     {default = pwd}
%   'roi_name': file stem of ROI file (e.g. lh.roi_name.label)
%     {default = 'v123'}
%   'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%     subjdir/subj should contain the FreeSurfer subject directory
%     {default = $SUBJECTS_DIR}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  04/07/11 by Don Hagler
% Last Mod: 07/30/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_out','view_ret.csh',[],...
  'fstem','pol',[],...
  'indir',pwd,[],...
  'roi_name','v123',[],...
  'label_dir',pwd,[],...
  'subjdir',[],[],...
  'forceflag',false,[false true],...
... % retinotopy data / tksurfer
  'surf','sphere',{'white','pial','inflated','sphere'},...
  'smooth',10,[0,100],...
  'fthresh',0,[],...
  'fmid',1.5,[],...
  'fslope',3,[],...
  'revflag',false,[false true],...
  'sph_rot',{[45 0 90],[45 -20 -90]},[],...
...
  'default_hemi','lh',[],...
});

if isempty(parms.subjdir)
  parms.subjdir = deblank(getenv('SUBJECTS_DIR'));
  if isempty(parms.subjdir)
    error('Cannot find SUBJECTS_DIR environment variable');
  end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.fname_out,'file') || parms.forceflag
  fid = fopen(parms.fname_out,'wt');
  if fid==-1
    error('failed to open %s for writing',parms.fname_out);
  end;
  fprintf(fid,'#!/bin/csh -f\n');
  fprintf(fid,'# automatically generated script for viewing retinotopy data\n');
  fprintf(fid,'\n');
  fprintf(fid,'if ($#argv == 0) then\n');
  fprintf(fid,'  set hemi = %s\n',parms.default_hemi);
  fprintf(fid,'else\n');
  fprintf(fid,'  set hemi = $argv[1]\n');
  fprintf(fid,'endif\n');
  fprintf(fid,'\n');
  fprintf(fid,'if ($#argv > 1) then\n');
  fprintf(fid,'  set save_flag = $argv[2]\n');
  fprintf(fid,'else\n');
  fprintf(fid,'  set save_flag = 0\n');
  fprintf(fid,'endif\n');
  fprintf(fid,'\n');
  fprintf(fid,'setenv SUBJECTS_DIR %s\n',parms.subjdir);
  fprintf(fid,'set subj = %s\n',subj);
  fprintf(fid,'set fstem = %s\n',parms.fstem);
  fprintf(fid,'set indir = %s\n',parms.indir);
  fprintf(fid,'\n');
  fprintf(fid,'set realfile = $fstem''_r''\n');
  fprintf(fid,'set imagfile = $fstem''_i''\n');
  fprintf(fid,'\n');
  fprintf(fid,'set smooth = %d\n',parms.smooth);
  fprintf(fid,'set fthresh = %0.3f\n',parms.fthresh);
  fprintf(fid,'set fmid = %0.3f\n',parms.fmid);
  fprintf(fid,'set fslope = %0.3f\n',parms.fslope);
  fprintf(fid,'set surf = %s\n',parms.surf);
  fprintf(fid,'set sph_rot_lh = (%s)\n',sprintf('%d ',parms.sph_rot{1}));
  fprintf(fid,'set sph_rot_rh = (%s)\n',sprintf('%d ',parms.sph_rot{2}));
  fprintf(fid,'\n');
  fprintf(fid,'set label = %s/$hemi".%s.label"\n',...
    parms.label_dir,parms.roi_name);
  fprintf(fid,'if (-e $label) then\n');
  fprintf(fid,'  set labelstr = "-label $label"\n');
  fprintf(fid,'else\n');
  fprintf(fid,'  set labelstr = " "\n');
  fprintf(fid,'endif\n');
  fprintf(fid,'\n');
  fprintf(fid,'if $save_flag then\n');
%  fprintf(fid,'  set savestr = "-savetiff -outstem $fstem -offscreen"\n');
  fprintf(fid,'  set savestr = "-savetiff -outstem $fstem"\n');
  fprintf(fid,'else\n');
  fprintf(fid,'  set savestr = " "\n');
  fprintf(fid,'endif\n');
  fprintf(fid,'\n');
  fprintf(fid,'if ($surf == "sphere") then\n');
  fprintf(fid,'  set view = pos\n');
  fprintf(fid,'  if ($hemi == "lh") then\n');
  fprintf(fid,'    set rotx = $sph_rot_lh[1]\n');
  fprintf(fid,'    set roty = $sph_rot_lh[2]\n');
  fprintf(fid,'    set rotz = $sph_rot_lh[3]\n');
  fprintf(fid,'  else\n');
  fprintf(fid,'    set rotx = $sph_rot_rh[1]\n');
  fprintf(fid,'    set roty = $sph_rot_rh[2]\n');
  fprintf(fid,'    set rotz = $sph_rot_rh[3]\n');
  fprintf(fid,'  endif\n');
  fprintf(fid,'else\n');
  fprintf(fid,'  set view = med\n');
  fprintf(fid,'  set rotx = 0\n');
  fprintf(fid,'  set roty = 0\n');
  fprintf(fid,'  set rotz = 0\n');
  fprintf(fid,'endif\n');
  fprintf(fid,'\n');
  fprintf(fid,'fs_surfmgh.csh $subj $realfile $hemi \\\n');
  fprintf(fid,'  -imagfile $imagfile \\\n');
  fprintf(fid,'  -surf $surf \\\n');
  fprintf(fid,'  -indir $indir \\\n');
  fprintf(fid,'  -smooth $smooth \\\n');
  fprintf(fid,'  -fthresh $fthresh \\\n');
  fprintf(fid,'  -fmid $fmid \\\n');
  fprintf(fid,'  -fslope $fslope \\\n');
  fprintf(fid,'  -view $view \\\n');
  fprintf(fid,'  -rotx $rotx -roty $roty -rotz $rotz \\\n');
  if ~isempty(regexp(parms.fstem,'pol'))
    fprintf(fid,'  -polar');
  elseif ~isempty(regexp(parms.fstem,'ecc'))
    fprintf(fid,'  -eccen');
  end;
  fprintf(fid,' $labelstr $savestr');
  if parms.revflag
    fprintf(fid,' -revphase \\\n');
  else
    fprintf(fid,' \\\n');
  end;
  fprintf(fid,'\n');
  fprintf(fid,'rm tempsurf.tcl.*\n');
  fprintf(fid,'\n');
  fclose(fid);

  % change permissions so group can read, write, and execute
  cmd = sprintf('chmod ug+rwx %s',parms.fname_out);
  [status,result] = unix(cmd);
  if status
    fprintf('%s: WARNING: failed to set permissions for %s:\n%s\n',...
      mfilename,parms.fname_out,result);
  end;
end;

