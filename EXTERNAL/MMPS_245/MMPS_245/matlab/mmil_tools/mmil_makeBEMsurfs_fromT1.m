function mmil_makeBEMsurfs_fromT1(fname_in,outdir,forceflag,options)
%function mmil_makeBEMsurfs_fromT1(fname_in,outdir,forceflag,[options])
%
% Early Mod: 02/24/10 by Don Hagler
% Last Mod:  12/29/14 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag=0; end;
if ~exist('options','var'), options=[]; end;

bemsurfs_cmd = 'mri_watershed';
matlab_cmd = 'matlab -nojvm -nosplash -r';
nverts = 2562; % same as ico4

[indir,fstem,fext] = fileparts(fname_in);
fname_in = [fstem fext];
if isempty(indir), indir = pwd; end;
full_fname_in = [indir '/' fstem fext];
if ~exist(full_fname_in,'file')
  error('file %s not found',full_fname_in);
end;

mmil_mkdir(outdir);

outdir2 = [outdir '/watershed'];
[success,msg] = mkdir(outdir2);
if ~success
  error('unable to create dir %s:\n%s',mfilename,outdir2,msg);
end;

fname_surf1 = 'lh.bem_inner_skull_surface';
fname_surf2 = 'lh.bem_outer_skull_surface';
fname_surf3 = 'lh.bem_outer_skin_surface';
fname_brain = 'brain.mgh';
fname_tri1 = 'inner_skull4.tri';
fname_tri2 = 'outer_skull4.tri';
fname_tri3 = 'outer_scalp4.tri';

full_fname_tri1 = [outdir '/' fname_tri1];
full_fname_tri2 = [outdir '/' fname_tri2];
full_fname_tri3 = [outdir '/' fname_tri3];

if exist(full_fname_tri1,'file') && exist(full_fname_tri2,'file') &&...
   exist(full_fname_tri3,'file') && ~forceflag, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_csh = sprintf('%s/make_bem_surfs.csh',outdir);
fid = fopen(fname_csh,'wt');
if fid<0
  fprintf('%s: ERROR: unable to create %s\n',mfilename,fname_csh);
  return;
end;
fprintf(fid,'#!/bin/csh -f\n');
fprintf(fid,'# make_bem_surfs.csh  %s\n',date);
fprintf(fid,'set forceflag = %d\n',forceflag);
fprintf(fid,'set cwd = `pwd`\n');
fprintf(fid,'set indir = "%s"\n',indir);
fprintf(fid,'set outdir = "%s"\n',outdir);
fprintf(fid,'set outdir2 = "%s"\n',outdir2);
fprintf(fid,'cd $outdir2\n');
fprintf(fid,'set fname_surf1 = "$outdir2/%s"\n',fname_surf1);
fprintf(fid,'set fname_surf2 = "$outdir2/%s"\n',fname_surf2);
fprintf(fid,'set fname_surf3 = "$outdir2/%s"\n',fname_surf3);
% run mri_watershed
fprintf(fid,'if ((! -e $fname_surf1) || (! -e $fname_surf2) || (! -e $fname_surf3) || $forceflag) then\n');
fprintf(fid,'  %s -surf bem -useSRAS %s "$indir/%s" "$outdir2/%s"\n',...
  bemsurfs_cmd,options,fname_in,fname_brain);
fprintf(fid,'endif\n');
fprintf(fid,'if ((! -e $fname_surf1) || (! -e $fname_surf2) || (! -e $fname_surf3) || $forceflag) then\n');
fprintf(fid,'  exit\n');
fprintf(fid,'endif\n');
% convert to tri file
fprintf(fid,'set fname_tri1 = "$outdir/%s"\n',fname_tri1);
fprintf(fid,'set fname_tri2 = "$outdir/%s"\n',fname_tri2);
fprintf(fid,'set fname_tri3 = "$outdir/%s"\n',fname_tri3);
fprintf(fid,'set nverts = %d\n',nverts);
fprintf(fid,'if ((! -e $fname_tri1) || (! -e $fname_tri2) || (! -e $fname_tri3) || $forceflag) then\n');
fprintf(fid,'  %s "try, fs_convert_fsurf2tri(''$fname_surf1'',''$fname_tri1'',$nverts); fs_convert_fsurf2tri(''$fname_surf2'',''$fname_tri2'',$nverts); fs_convert_fsurf2tri(''$fname_surf3'',''$fname_tri3'',$nverts); catch, disp(lasterr), end; exit;"\n',...
  matlab_cmd);
fprintf(fid,'endif\n');
fprintf(fid,'cd $cwd\n');
fclose(fid);

cmd = sprintf('source %s | tee %s/stdout.log',fname_csh,outdir);
[status,result] = mmil_unix(cmd);
if status
  fprintf('%s: WARNING: %s exited with errors:\n%s\n',...
    mfilename,bemsurfs_cmd,result);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_csh = sprintf('%s/view_inner_skull.csh',outdir);
fid = fopen(fname_csh,'wt');
if fid<0
  error('unable to create %s',fname_csh);
end;
fprintf(fid,'#!/bin/csh -f\n');
fprintf(fid,'# view_inner_skull.csh  %s\n',date);
fprintf(fid,'set indir = "%s"\n',indir);
fprintf(fid,'set outdir = "%s"\n',outdir);
fprintf(fid,'tkmedit -f "$indir/%s" -surface "$outdir/%s"\n',fname_in,fname_tri1);
fclose(fid);

fname_csh = sprintf('%s/view_outer_skull.csh',outdir);
fid = fopen(fname_csh,'wt');
if fid<0
  error('unable to create %s',fname_csh);
end;
fprintf(fid,'#!/bin/csh -f\n');
fprintf(fid,'# view_outer_skull.csh  %s\n',date);
fprintf(fid,'set indir = "%s"\n',indir);
fprintf(fid,'set outdir = "%s"\n',outdir);
fprintf(fid,'tkmedit -f "$indir/%s" -surface "$outdir/%s"\n',fname_in,fname_tri2);
fclose(fid);

fname_csh = sprintf('%s/view_outer_scalp.csh',outdir);
fid = fopen(fname_csh,'wt');
if fid<0
  error('unable to create %s',fname_csh);
end;
fprintf(fid,'#!/bin/csh -f\n');
fprintf(fid,'# view_outer_scalp.csh  %s\n',date);
fprintf(fid,'set indir = "%s"\n',indir);
fprintf(fid,'set outdir = "%s"\n',outdir);
fprintf(fid,'tkmedit -f "$indir/%s" -surface "$outdir/%s"\n',fname_in,fname_tri3);
fclose(fid);

cmd = sprintf('chmod +x %s/*.csh',outdir);
[status,result] = mmil_unix(cmd);
if status
  fprintf('%s: WARNING: chmod exited with errors:\n%s\n',mfilename,result);
end;

