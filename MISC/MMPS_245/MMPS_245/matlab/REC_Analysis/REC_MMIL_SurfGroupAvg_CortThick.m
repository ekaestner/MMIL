function REC_MMIL_SurfGroupAvg_CortThick(ProjID,analysis_outdir)
%function REC_MMIL_SurfGroupAvg_CortThick(ProjID,analysis_outdir)
%
% Optional Parameters:
%   analysis_outdir: FSrecon subdir where stats are 
%     {default: 'analysis'}
%
% Created:  06/02/09 by Alain Koyama
% Last Mod: 11/08/12 Don Hagler
%

if ~exist('ProjID','var'), ProjID = []; end;
if ~exist('analysis_outdir','var'), analysis_outdir = 'analysis'; end;

if isempty(ProjID)
   error('Empty ProjID');
end

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);

hemilist = {'lh','rh'};
sphsmooth = 176;

dirlist_fs = dir(sprintf('%s/FSURF*',RootDirs.fsurf));
outdir = sprintf('%s/MetaData/%s/SurfGroupAvgs',RootDirs.home,ProjID);
[success,message]=mkdir(outdir);
if ~success
    fprintf('%s: ERROR: unable to create output dir %s\n',...
        mfilename,outdir);
    return;
end;
fname = sprintf('%s/SurfGroupAvg_CortThick.log',outdir);
flog = fopen(fname,'wt');

fname = sprintf('%s/MetaData/%s/REC_%s_StudyInfo.mat',RootDirs.home,ProjID,ProjID);
if exist(fname,'file') % load subject meta data if it exists
    load(fname);
else
    fprintf('%s: ERROR: %s not found\n',mfilename,fname);
    return;
end
% find group data
Diagsall = {StudyInfo.Group};
Diagsbak = Diagsall;
Diagsall(cellfun(@isempty,Diagsall)) = []; % if some subjects are missing metadata info
Diags = unique(Diagsall);

% diffs btw all groups will be made by default
contrastlist = nchoosek(1:length(Diags),2);
total_sum = 0;
fprintf('%s: loading data...\n',mfilename);
fprintf(flog,'Studies included in surface group average:\n');
for d=1:length(Diags) % loop thru each diag    
  diag = Diags{d};
  S_ind = find(strcmp(Diagsbak,diag));
  for h = 1:length(hemilist) % loop thru each hemi
    hemi = hemilist{h};
    fprintf(flog,'diagnosis %s, hemi %s\n',diag,hemi);
    info{d}.info{h}.nsum = 0;
    info{d}.info{h}.vecsum = 0;
    info{d}.info{h}.vecsum2 = 0;
    for s = 1:length(S_ind) % loop thru each subject of specific diag
      SubjID = StudyInfo(S_ind(s)).SubjID;
      n = regexp({dirlist_fs.name},['FSURF_' SubjID '_\d{8}.+_1$'],'match');
      n = [n{:}];
      if isempty(n)
          fprintf('%s: WARNING: no Recon container found for subject %s\n',mfilename,SubjID);
          continue;
      end
      ContainerDir = n{1};
      ContainerPath = sprintf('%s/%s',RootDirs.fsurf,ContainerDir);

      fname = sprintf('%s/%s/thickness-mbmask-sphere-sm%d-%s.mgh',...
          ContainerPath,analysis_outdir,sphsmooth,hemi);
      if ~exist(fname,'file')
          fprintf('%s: WARNING: %s not found\n',mfilename,fname);
          continue;
      end;
      vec = mmil_rowvec(fs_load_mgh(fname));
      if isempty(vec)
          fprintf('%s: WARNING: unable to correctly read %s\n',mfilename,fname);
          continue;
      end;
      info{d}.info{h}.nsum = info{d}.info{h}.nsum + 1;
      info{d}.info{h}.vecsum = info{d}.info{h}.vecsum + vec;
      info{d}.info{h}.vecsum2 = info{d}.info{h}.vecsum2 + vec.^2;
      fprintf(flog,'%s\n',ContainerDir);
    end
    nsum = info{d}.info{h}.nsum;
    fprintf(flog,'N for diagnosis %s, hemi %s: %d\n',diag,hemi,nsum);
    if nsum<=1
        fprintf('%s: WARNING: <=1 subjects included for diag %s, hemi %s\n',...
            mfilename,diag,hemi);
        nsum=2;
    end;
    info{d}.info{h}.vecmean = info{d}.info{h}.vecsum / (eps+nsum);
    info{d}.info{h}.vecstd = sqrt((nsum*info{d}.info{h}.vecsum2 - info{d}.info{h}.vecsum.^2)./(eps+nsum*(nsum-1)));
    info{d}.info{h}.vecstderr = info{d}.info{h}.vecstd/sqrt(eps+nsum);
  end
  total_sum = total_sum + nsum;
end
fprintf('%s: TOTAL n = %d\n',mfilename,total_sum);
fclose(flog);

fprintf('%s: saving group averages...\n',mfilename);
for d=1:length(Diags)
  for h = 1:length(hemilist)
    hemi = hemilist{h};
    vecmean = info{d}.info{h}.vecmean;
    vecstderr = info{d}.info{h}.vecstderr;
    vectval = vecmean./vecstderr;
    fname = sprintf('%s/CortThick_mean_%s-sm%d-%s.mgh',...
        outdir,Diags{d},sphsmooth,hemi);
    fs_save_mgh(vecmean,fname,eye(4));
    fname = sprintf('%s/CortThick_tval_%s-sm%d-%s.mgh',...
        outdir,Diags{d},sphsmooth,hemi);
    fs_save_mgh(vectval,fname,eye(4));
    fname = sprintf('%s/CortThick_stderr_%s-sm%d-%s.mgh',...
        outdir,Diags{d},sphsmooth,hemi);
    fs_save_mgh(vecstderr,fname,eye(4));
  end
end

for c = 1:size(contrastlist,1)
  i1 = contrastlist(c,1);
  i2 = contrastlist(c,2);
  for h = 1:length(hemilist)
    hemi = hemilist{h};

    % absolute difference
    vecdiff = info{i1}.info{h}.vecmean - info{i2}.info{h}.vecmean;
    vecstderr = sqrt(info{i1}.info{h}.vecstderr.^2+info{i2}.info{h}.vecstderr.^2);
    vectval = vecdiff./vecstderr;
    fname = sprintf('%s/CortThick_diff_mean_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vecdiff,fname,eye(4));
    fname = sprintf('%s/CortThick_diff_tval_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vectval,fname,eye(4));
    fname = sprintf('%s/CortThick_diff_stderr_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vecstderr,fname,eye(4));

    % fractional difference
    vecdiff = (info{i1}.info{h}.vecmean - info{i2}.info{h}.vecmean)./info{1}.info{h}.vecmean; % Difference in percent of NL
    vecstderr = sqrt((info{i1}.info{h}.vecstderr./info{1}.info{h}.vecmean).^2+...
        (info{i2}.info{h}.vecstderr./info{1}.info{h}.vecmean).^2);
    vectval = vecdiff./vecstderr;
    fname = sprintf('%s/CortThick_pc_diff_mean_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vecdiff,fname,eye(4));
    fname = sprintf('%s/CortThick_pc_diff_tval_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vectval,fname,eye(4));
    fname = sprintf('%s/CortThick_pc_diff_stderr_%s_%s-sm%d-%s.mgh',...
        outdir,Diags{i1},Diags{i2},sphsmooth,hemi);
    fs_save_mgh(vecstderr,fname,eye(4));
  end
end
