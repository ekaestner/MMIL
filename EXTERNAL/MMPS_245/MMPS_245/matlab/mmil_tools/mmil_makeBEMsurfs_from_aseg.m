function mmil_makeBEMsurfs_from_aseg(fname_T1,fname_aseg,outdir,forceflag)
%function mmil_makeBEMsurfs_from_aseg(fname_T1,fname_aseg,outdir,forceflag)
%
% Last Mod:  09/07/11 by Don Hagler
%

if (~mmil_check_nargs(nargin,3)), return; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag=0; end;

[indir,fstem,fext] = fileparts(fname_T1);
fname_T1 = [fstem fext];
if isempty(indir), indir = pwd; end;
full_fname_T1 = [indir '/' fstem fext];
if ~exist(full_fname_T1,'file')
  error('file %s not found',full_fname_T1);
end;

[success,msg] = mkdir(outdir);
if ~success
  error('unable to create dir %s:\n%s',outdir,msg);
end;

fname_mask = 'skull_mask.mgz';
fname_surf = 'inner_skull_aseg.tri';

full_fname_mask = [outdir '/' fname_mask];
full_fname_surf = [outdir '/' fname_surf];

smooth1 = 20;
thresh1 = 0.5;
smooth2 = 40;
thresh2 = 0.2;
smooth3 = 10;

vol_mask = mmil_dilate_mask([],...
  'fname_in',fname_aseg,'fname_out',full_fname_mask,...
  'smooth1',smooth1,'thresh1',thresh1,...
  'smooth2',smooth2,'thresh2',thresh2,...
  'smooth3',smooth3,...
  'forceflag',forceflag);

masksurf = mmil_masksurf(vol_mask,...
  'fname_out',full_fname_surf,'forceflag',forceflag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_csh = sprintf('%s/view_inner_skull_aseg.csh',outdir);
fid = fopen(fname_csh,'wt');
if fid<0
  error('unable to create %s',fname_csh);
end;
fprintf(fid,'#!/bin/csh -f\n');
fprintf(fid,'# view_inner_skull_aseg.csh  %s\n',date);
fprintf(fid,'set indir = "%s"\n',indir);
fprintf(fid,'set outdir = "%s"\n',outdir);
fprintf(fid,'tkmedit -f "$indir/%s" -surface "$outdir/%s"\n',fname_T1,fname_surf);
fclose(fid);
cmd = sprintf('chmod +x %s',fname_csh);
[status,result] = unix(cmd);
if status
  fprintf('%s: WARNING: chmod exited with errors:\n%s\n',mfilename,result);
end;
