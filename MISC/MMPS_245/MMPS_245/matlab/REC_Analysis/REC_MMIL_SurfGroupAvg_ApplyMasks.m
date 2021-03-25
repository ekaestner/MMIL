function REC_MMIL_SurfGroupAvg_ApplyMasks(ProjID)
%function REC_MMIL_SurfGroupAvg_ApplyMasks(ProjID)
%
% Required Inputs:
% ProjID - name of Recharge project to run (default = run all projects)
%
% Last Mod:  08/16/10 by Cooper Roddey
%

mask_roinums = [1,5];

if isempty(ProjID)
   error('Empty ProjID');
end

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);
homedir = RootDirs.home;

for p=1:length(ProjInfo) % loop thru projects
    outdir = sprintf('%s/MetaData/%s/SurfGroupAvgs',homedir,ProjInfo(p).ProjID);
    if ~exist(outdir,'dir')
      fprintf('%s: ERROR: input/output dir %s not found\n',...
        mfilename,outdir);
      return;
    end;

    subjdir = getenv('FREESURFER_HOME');
    if isempty(subjdir)
      fprintf('%s: ERROR: FREESURFER_HOME not defined as an environment variable\n',mfilename);
      return;
    end;
    subjdir = sprintf('%s/subjects',subjdir);
    subj = 'fsaverage';
    fullsubj = sprintf('%s/%s',subjdir,subj);
    if ~exist(fullsubj,'dir')
      fprintf('%s: ERROR: average subject %s not found\n',mfilename,subj);
      return;
    end;

    flist = dir(sprintf('%s/*.mgh',outdir));
    for f=1:length(flist)
      infile = sprintf('%s/%s',outdir,flist(f).name);
      k=findstr(infile,'h.mgh');
      if isempty(k), continue; end;
      if ~isempty(regexp(flist(f).name,'mbmask')), continue; end;
      k=k(end);
      instem = infile(1:k-3);
      hemi = infile(k-1:k);
      outfile = sprintf('%s-mbmask-%s.mgh',instem,hemi);
      fs_mask_surfmgh_with_aparc(subj,hemi,infile,outfile,subjdir,mask_roinums);
      if ~exist(outfile,'file');
        fprintf('%s: ERROR: failed to apply mask\n',mfilename);
        return;
      end;
    end;
end


