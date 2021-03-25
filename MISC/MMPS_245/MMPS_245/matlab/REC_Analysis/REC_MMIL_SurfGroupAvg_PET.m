function REC_MMIL_SurfGroupAvg_PET(ProjID,normflag,qcflag)
%function REC_MMIL_SurfGroupAvg_PET(ProjID,normflag,qcflag)
%
% Required Input:
%   ProjID: name of Recharge project to run (default = run all projects)
%
% Optional Input:
%   normflag: normalize to pons
%    (0 = absolute values, 1 = norm to pons, 2 = both)
%    {default = 1}
%   qcflag: average only good qc
%    {default = 1}
%
% Created:  01/13/09 by Alain Koyama
% Last Mod: 09/13/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('normflag','var'), normflag = 1; end
if ~exist('qcflag','var'), qcflag = 1; end
if isempty(ProjID), error('empty ProjID'); end

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);
homedir = RootDirs.home;

hemilist = {'lh','rh'};
sphsmooth = 176;

dirlist_fs = dir(sprintf('%s/FSURF*',RootDirs.fsurf));
dirlist_pet = dir(sprintf('%s/PETPROC*',RootDirs.proc_pet));
outdir = sprintf('%s/MetaData/%s/SurfGroupAvgs',homedir,ProjID);
[success,message]=mkdir(outdir);
if ~success
  fprintf('%s: ERROR: unable to create output dir %s\n',...
    mfilename,outdir);
  return;
end;
fname = sprintf('%s/SurfGroupAvg_PET.log',outdir);
flog = fopen(fname,'wt');

fname = sprintf('%s/MetaData/%s/REC_%s_StudyInfo.mat',homedir,ProjID,ProjID);
if exist(fname,'file') % load subject meta data if it exists
  load(fname);
else
  fprintf('%s: ERROR: %s not found\n',mfilename,fname);
  return;
end
% find group data
Diagsall = {StudyInfo.Group};
Diags = unique(Diagsall(~cellfun(@isempty,Diagsall)));

% diffs btw all groups will be made by default
contrastlist = nchoosek(1:length(Diags),2);

switch normflag
  case 0, normflags = 0;
  case 1, normflags = 1;
  case 2, normflags = [0 1];
end
for nn=1:length(normflags)
  total_sum = 0;
  normflag = normflags(nn);

  fprintf('%s: loading data (normflag = %d)...\n',mfilename,normflag);
  fprintf(flog,'Studies included in surface group average (normflag = %d)\n',normflag);
  for d=1:length(Diags)
    diag = Diags{d};
    S_ind = find(strcmp(Diagsall,diag));
    for h = 1:length(hemilist)
      createinfo = 0;
      hemi = hemilist{h};
      fprintf(flog,'diagnosis %s, hemi %s\n',diag,hemi);
      for s = 1:length(S_ind)
        SubjID = StudyInfo(S_ind(s)).SubjID;
        if qcflag
          PETRegQC = StudyInfo(S_ind(s)).PETRegQC;
          if strcmpi(PETRegQC,'Bad')
            fprintf('%s: WARNING: bad PETRegQC for subject %s\n',mfilename,SubjID);
            continue;
          end
        end

        % does this subject have a recon?
        n = regexp({dirlist_fs.name},['FSURF_' SubjID '_\d{8}.+'],'match');
        n = [n{:}];

        if isempty(n)
          fprintf('%s: WARNING: no FreeSurfer recon container found for subject %s\n',mfilename,SubjID);
          continue;
        end

        % does this subject have PET data?
        n = regexp({dirlist_pet.name},['PETPROC_' SubjID '_\d{8}.+'],'match');
        n = [n{:}];
        if isempty(n)
          fprintf('%s: WARNING: no PET container found for subject %s\n',mfilename,SubjID);
          continue;
        end
        PETContainer = n{1};
        PETContainerPath = sprintf('%s/%s',,PETContainer);

        volnames = {};

        % search for single or multiple PET series
        PET_files = dir(sprintf('%s/PET_reg*txt',PETContainerPath));
        if isempty(PET_files)
          fprintf('%s: ERROR: PET reg file not found\n',mfilename);
          return;
        else
          for j=1:length(PET_files)
            volnames{end+1} = char(regexp(PET_files(j).name,'(?<=PET_reg_).+(?=\.txt)','match'));
          end
        end
        if ~createinfo % create info var at first iteration
          for numvols=1:length(volnames)
            info{d}.info{h}.nsum{numvols} = 0;
            info{d}.info{h}.vecsum{numvols} = 0;
            info{d}.info{h}.vecsum2{numvols} = 0;
          end
        end
        createinfo = 1;
        for numvols=1:length(volnames) % loop thru PET series
          fname_roi = sprintf('%s/PET_%s_normroi_data.mat',PETContainerPath,volnames{numvols});
          if ~exist(fname_roi,'file') & normflag
            fprintf('%s: WARNING: %s not found\n',mfilename,fname_roi);
            continue;
          end;
          fname = sprintf('%s/PET_reg_%s-paint-mbmask-sphere-sm%d-%s.mgh',...
            PETContainerPath,volnames{numvols},sphsmooth,hemi);
          if ~exist(fname,'file')
            fprintf('%s: WARNING: %s not found\n',mfilename,fname);
            continue;
          end;
          vec = mmil_rowvec(fs_load_mgh(fname));
          roi_data = [];
          if normflag
            load(fname_roi);
            if isempty(roi_data)
              fprintf('%s: ERROR: empty roi_data for %s\n',mfilename,fname_roi);
              continue;
            end;
            normval = roi_data(1).avg;
            if isempty(normval) | normval==0
              fprintf('%s: WARNING: normval is empty or zero for %s\n',...
                mfilename,fname_roi);
              continue;
            end;
            vec = vec / normval;
          end
          info{d}.info{h}.nsum{numvols} = info{d}.info{h}.nsum{numvols} + 1;
          info{d}.info{h}.vecsum{numvols} = info{d}.info{h}.vecsum{numvols} + vec;
          info{d}.info{h}.vecsum2{numvols} = info{d}.info{h}.vecsum2{numvols} + vec.^2;
        end
        fprintf(flog,'%s\n',PETContainer);
      end
      for numvols=1:length(volnames)
        nsum = info{d}.info{h}.nsum{numvols};
        fprintf(flog,'N for diagnosis %s, %s hemi %s: %d\n',diag,volnames{numvols}(2:end),hemi,nsum);
        if nsum<=1
          fprintf('%s: WARNING: <=1 subjects included for hemi %s\n',...
            mfilename,hemi);
          nsum=2;
        end
        info{d}.info{h}.vecmean{numvols} = info{d}.info{h}.vecsum{numvols} / (eps+nsum);
        info{d}.info{h}.vecstd{numvols} = sqrt((nsum*info{d}.info{h}.vecsum2{numvols} - info{d}.info{h}.vecsum{numvols}.^2)./(eps+nsum*(nsum-1)));
        info{d}.info{h}.vecstderr{numvols} = info{d}.info{h}.vecstd{numvols}/sqrt(eps+nsum);
      end
    end
    total_sum = total_sum + nsum;
  end
  fprintf('%s: TOTAL n = %d\n',mfilename,total_sum)

  fprintf('%s: saving group averages...\n',mfilename);
  types = {'mean','tval','stderr'};

  for d=1:length(Diags)
    for h = 1:length(hemilist)
      for numvols=1:length(volnames)
        hemi = hemilist{h};
        vecmean = info{d}.info{h}.vecmean{numvols};
        vecstderr = info{d}.info{h}.vecstderr{numvols};
        vectval = vecmean./vecstderr;
        for i=1:length(types)
          if normflag
            fname = sprintf('%s/PET_%s_norm_%s_%s-sm%d-%s.mgh',...
            outdir,volnames{numvols},types{i},Diags{d},sphsmooth,hemi);
          else
            fname = sprintf('%s/PET_%s_%s_%s-sm%d-%s.mgh',...
            outdir,volnames{numvols},types{i},Diags{d},sphsmooth,hemi);
          end
          fs_save_mgh(eval(['vec' types{i}]),fname,eye(4));
        end
      end
    end
  end

  for c = 1:size(contrastlist,1)
    i1 = contrastlist(c,1);
    i2 = contrastlist(c,2);
    for h = 1:length(hemilist)
      hemi = hemilist{h};
      % fractional difference
      vecmean = (info{i1}.info{h}.vecmean{numvols} - info{i2}.info{h}.vecmean{numvols})./info{1}.info{h}.vecmean{numvols}; % Difference in percent of NL
      vecstderr = sqrt((info{i1}.info{h}.vecstderr{numvols}./info{1}.info{h}.vecmean{numvols}).^2+...
          (info{i2}.info{h}.vecstderr{numvols}./info{1}.info{h}.vecmean{numvols}).^2);
      vectval = vecmean./vecstderr;
      for i=1:length(types)
        if normflag
          fname = sprintf('%s/PET_%s_norm_pc_diff_%s_%s_%s-sm%d-%s.mgh',...
              outdir,volnames{numvols},types{i},Diags{i1},Diags{i2},sphsmooth,hemi);
          fs_save_mgh(eval(['vec' types{i}]),fname,eye(4));
        else
          fname = sprintf('%s/PET_%s_pc_diff_%s_%s_%s-sm%d-%s.mgh',...
              outdir,volnames{numvols},types{i},Diags{i1},Diags{i2},sphsmooth,hemi);
          fs_save_mgh(eval(['vec' types{i}]),fname,eye(4));
        end
      end
    end
  end
end
fclose(flog);
