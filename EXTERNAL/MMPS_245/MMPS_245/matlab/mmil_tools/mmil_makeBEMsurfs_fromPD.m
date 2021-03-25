function mmil_makeBEMsurfs_fromPD(fname_in,outdir,forceflag)
%function mmil_makeBEMsurfs_fromPD(fname_in,outdir,forceflag)
%
% Created:  01/19/09 by Don Hagler
% Last Mod: 05/20/14 by Don Hagler
%

if ~mmil_check_nargs(nargin, 2), return; end;
if ~exist('forceflag','var'), forceflag=0; end;

bemsurfs_cmd = 'find_bem_surfs';

tri_file_dir = [getenv('MMPS_DIR') '/cpp/find_bem_surfs/BEM_TRI_FILES/'];
if ~exist(tri_file_dir,'dir')
  error('tri file dir %s not found',tri_file_dir);
end;

ico = 4;
switch ico
  case 4
    ISk_nMove = 300;
    ISk_nRest = 30;
    ISk_nRelax = 1;
    ISk_cTang = 0.07;
    ISk_cNorm = 0.07;
    ISk_cMRI = 0.7;
    ISk_cRepell = 0.0;
    ISk_thresh = 80.0;
    OSk_nMove = 150;
    OSk_nRest = 30;
    OSk_nRelax = 2;
    OSk_cTang = 0.18;
    OSk_cNorm = 0.18;
    OSk_cMRI = 0.7;
    OSk_cRepell = 1.0;
    OSk_thresh = 90.0;
    OSc_nMove = 200;
    OSc_nRest = 30;
    OSc_nRelax = 1;
    OSc_cTang = 0.04;
    OSc_cNorm = 0.02;
    OSc_cMRI = 0.5;
    OSc_cRepell = 0.5;
    OSc_thresh = 35.0;
  case 5
    ISk_nMove = 300;
    ISk_nRest = 30;
    ISk_nRelax = 2;
    ISk_cTang = 0.1;
    ISk_cNorm = 0.1;
    ISk_cMRI = 0.3;
    ISk_cRepell = 0.0;
    ISk_thresh = 80.0;
    OSk_nMove = 100;
    OSk_nRest = 30;
    OSk_nRelax = 3;
    OSk_cTang = 0.2;
    OSk_cNorm = 0.2;
    OSk_cMRI = 0.2;
    OSk_cRepell = 1.0;
    OSk_thresh = 70.0;
    OSc_nMove = 200;
    OSc_nRest = 30;
    OSc_nRelax = 1;
    OSc_cTang = 0.04;
    OSc_cNorm = 0.02;
    OSc_cMRI = 0.5;
    OSc_cRepell = 0.5;
    OSc_thresh = 35.0;
  otherwise
    return;
end;

[indir,fstem,fext] = fileparts(fname_in);
fname_in = [fstem fext];
if isempty(indir), indir = pwd; end;
full_fname_in = [indir '/' fstem fext];
if ~exist(full_fname_in,'file')
  fprintf('%s: ERROR: %s not found\n',mfilename,full_fname_in);
  return;
end;

mmil_mkdir(outdir);

if strcmp(fext,'.mgz')
  fname_tmp = [outdir '/' fstem '.mgh'];
  if ~exist(fname_tmp,'file') || forceflag
    fs_copy_mgh(full_fname_in,fname_tmp);
  end;
  fname_in = [fstem '.mgh'];
  indir = outdir;
end;

fname_tri1 = sprintf('inner_skull%d.tri',ico);
fname_tri2 = sprintf('outer_skull%d.tri',ico);
fname_tri3 = sprintf('outer_scalp%d.tri',ico);

full_fname_tri1 = [outdir '/' fname_tri1];
full_fname_tri2 = [outdir '/' fname_tri2];
full_fname_tri3 = [outdir '/' fname_tri3];

if exist(full_fname_tri1,'file') && exist(full_fname_tri2,'file') &&...
   exist(full_fname_tri3,'file') && ~forceflag, return; end;

fname_csh = sprintf('%s/make_bem_surfs.csh',outdir);
fid = fopen(fname_csh,'wt');
if fid<0
  error('unable to create %s\n',fname_csh);
end;
fprintf(fid,'#!/bin/csh -f\n');
fprintf(fid,'# make_bem_surfs.csh  %s\n',date);
fprintf(fid,'set forceflag = %d\n',forceflag);
fprintf(fid,'set indir = "%s"\n',indir);
fprintf(fid,'set outdir = "%s"\n',outdir);
fprintf(fid,'set fname_tri1 = "$outdir/%s"\n',fname_tri1);
fprintf(fid,'set fname_tri2 = "$outdir/%s"\n',fname_tri2);
fprintf(fid,'set fname_tri3 = "$outdir/%s"\n',fname_tri3);
fprintf(fid,'if ((! -e $fname_tri1) || (! -e $fname_tri2) || (! -e $fname_tri3) || $forceflag) then\n');
fprintf(fid,'  %s "$indir/%s" \\\n',bemsurfs_cmd,fname_in);
fprintf(fid,'    -outdir $outdir \\\n');
fprintf(fid,'    -tri_file_dir %s \\\n',tri_file_dir);
fprintf(fid,'    -inner_skull_thresh %0.2f \\\n',ISk_thresh);
fprintf(fid,'    -outer_skull_thresh %0.2f \\\n',OSk_thresh);
fprintf(fid,'    -outer_scalp_thresh %0.2f \\\n',OSc_thresh);
fprintf(fid,'    -inner_skull_cTang %0.2f \\\n',ISk_cTang);
fprintf(fid,'    -inner_skull_cNorm %0.2f \\\n',ISk_cNorm);
fprintf(fid,'    -inner_skull_cMRI %0.2f \\\n',ISk_cMRI);
fprintf(fid,'    -outer_skull_cTang %0.2f \\\n',OSk_cTang);
fprintf(fid,'    -outer_skull_cNorm %0.2f \\\n',OSk_cNorm);
fprintf(fid,'    -outer_skull_cMRI %0.2f \\\n',OSk_cMRI);
fprintf(fid,'    -outer_scalp_cTang %0.2f \\\n',OSc_cTang);
fprintf(fid,'    -outer_scalp_cNorm %0.2f \\\n',OSc_cNorm);
fprintf(fid,'    -outer_scalp_cMRI %0.2f \\\n',OSc_cMRI);
fprintf(fid,'    -inner_skull_nMove %d \\\n',ISk_nMove);
fprintf(fid,'    -inner_skull_nRest %d \\\n',ISk_nRest);
fprintf(fid,'    -inner_skull_nRelax %d \\\n',ISk_nRelax);
fprintf(fid,'    -inner_skull_cRepell %0.2f \\\n',ISk_cRepell);
fprintf(fid,'    -outer_skull_nMove %d \\\n',OSk_nMove);
fprintf(fid,'    -outer_skull_nRest %d \\\n',OSk_nRest);
fprintf(fid,'    -outer_skull_nRelax %d \\\n',OSk_nRelax);
fprintf(fid,'    -outer_skull_cRepell %0.2f \\\n',OSk_cRepell);
fprintf(fid,'    -outer_scalp_nMove %d \\\n',OSc_nMove);
fprintf(fid,'    -outer_scalp_nRest %d \\\n',OSc_nRest);
fprintf(fid,'    -outer_scalp_nRelax %d \\\n',OSc_nRelax);
fprintf(fid,'    -outer_scalp_cRepell %0.2f \\\n',OSc_cRepell);
fprintf(fid,'    -ico %d\n',ico);
fprintf(fid,'\n  rm $outdir/*.par\n');
fprintf(fid,'endif\n');
if ico==5
  fprintf(fid,'triFileDecimator %s/inner_skull5.tri',outdir);
  fprintf(fid,'triFileDecimator %s/outer_skull5.tri',outdir);
  fprintf(fid,'triFileDecimator %s/outer_scalp5.tri',outdir);
end;
fclose(fid);

cmd = sprintf('source %s | tee %s/stdout.log',fname_csh,outdir);
[status,result] = unix(cmd);
if status
  fprintf('%s: WARNING: %s exited with errors:\n%s\n',...
    mfilename,bemsurfs_cmd,result);
end;

fname_csh = sprintf('%s/view_inner_skull.csh',outdir);
fid = fopen(fname_csh,'wt');
if fid<0
  error('unable to create %s',fname_csh);
end;
fprintf(fid,'#!/bin/csh -f\n');
fprintf(fid,'# view_inner_skull.csh  %s\n',date);
fprintf(fid,'set indir = "%s"\n',indir);
fprintf(fid,'set outdir = "%s"\n',outdir);
fprintf(fid,'tkmedit -f "$indir/%s" -surface "$outdir/%s"\n',...
  fname_in,fname_tri1);
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
fprintf(fid,'tkmedit -f "$indir/%s" -surface "$outdir/%s"\n',...
  fname_in,fname_tri2);
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
fprintf(fid,'tkmedit -f "$indir/%s" -surface "$outdir/%s"\n',...
  fname_in,fname_tri3);
fclose(fid);

cmd = sprintf('chmod +x %s/*.csh',outdir);
[status,result] = unix(cmd);
if status
  fprintf('%s: WARNING: chmod exited with errors:\n%s\n',mfilename,result);
end;

