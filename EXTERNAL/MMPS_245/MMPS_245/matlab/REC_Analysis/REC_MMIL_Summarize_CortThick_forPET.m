function REC_MMIL_Summarize_CortThick_forPET(ProjID,normflag)
%function REC_MMIL_Summarize_CortThick_forPET(ProjID,normflag)
%
% Required Parameters:
%   ProjID: project ID to run (e.g. 'REC_TEST')
%
% Optional Parameters:
%   normflag - [0|1|2] whether to normalize to pons
%    (0 = absolute values, 1 = norm to pons, 2 = both)
%    {default = 1}
%
% Created:  11/26/08 by Alain Koyama
% Rcnt Mod: 09/09/12 by Don Hagler
% Last Mod: 11/05/12 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
if isempty(ProjID), error('empty ProjID'); end;
if ~exist('normflag','var'), normflag = 1; end

RootDirs = REC_MMIL_RootDirs(ProjID);
outdir = sprintf('%s/MetaData/%s/ROI_Summaries',RootDirs.home,ProjID);
if ~exist(outdir,'dir'), mkdir(outdir); end;
dirlist_pet = dir(sprintf('%s/PETPROC*',RootDirs.proc_pet));
dirlist_fs = dir(sprintf('%s/FSURF*',RootDirs.fsurf));
if isempty(dirlist_pet)
    fprintf('%s: Project %s has no PET Processed directory\n',mfilename,ProjID);
    continue;
end
fname = sprintf('%s/MetaData/%s/REC_%s_StudyInfo.csv',RootDirs.home,ProjID,ProjID);
if exist(fname,'file') % load subject meta data if it exists
    metadata = mmil_readtext(fname,',','','"');
else
    metadata = {};
end
createinfo = 0;
for i = 1:length(dirlist_pet) % loop thru directories
  PETContainerPath = sprintf('%s/%s',RootDirs.proc_pet,dirlist_pet(i).name);

  % search for single or multiple PET series
  PET_files = dir(sprintf('%s/PET_reg*txt',PETContainerPath));
  if isempty(PET_files)
      fprintf('%s: ERROR: PET reg file not found\n',mfilename);
      return;
  else
      volnames = regexp({PET_files.name},'(?<=PET_reg_).+(?=\.txt)','match');
  end

  if ~createinfo % because of possible multiple series, create output and header at first iteration
      for numvols=1:length(volnames)
          outfile{numvols} = sprintf('%s/REC_%s_%s_CortThick_forPET.csv',outdir,ProjID,char(volnames{numvols}));
          fid=fopen(outfile{numvols},'wt');
          fclose(fid);
      end
      wrote_header_flag=zeros(1,length(volnames));
  end
  createinfo = 1;

  % cortical thickness
  for numvols=1:length(volnames)
    fid=fopen(outfile{numvols},'at');
    SubjID = char(regexp(dirlist_pet(i).name,'(?<=PETPROC_).+(?=_\d{8})','match'));
    StudyDate = char(regexp(dirlist_pet(i).name,'(?<=PETPROC_.+)\d{8}(?=\.)','match'));
    metaheaders = cell(1,size(metadata,2));
    if ~isempty(metadata) % if subject metadata exists, load headers
        for j=2:size(metadata,2)
            metaheaders{j} = metadata{1,j};
        end
    end
    S_ind = [];
    for j=2:size(metadata,1)
        if strcmp(SubjID,metadata{j,1}), S_ind = j; break; end;
    end

    studymatch = {}; % find corresponding fs recon
    for j=1:length(dirlist_fs) %find matching studies and get earliest one as baseline
        if regexp(dirlist_fs(j).name,SubjID)
            studymatch{end+1} = dirlist_fs(j).name;
        end
    end
    if ~isempty(studymatch)
        studymatch = sort(studymatch);
        ContainerPath = sprintf('%s/%s',RootDirs.fsurf,studymatch{1});
    else
        fprintf('%s: WARNING: no FS recon for baseline MRI... skipping %s\n',mfilename,SubjID);
        continue;
    end

    fname_roi = sprintf('%s/PET_%s_normroi_data.mat',PETContainerPath,char(volnames{numvols}));
    if ~exist(fname_roi,'file')
        fprintf('%s: WARNING: %s not found\n',mfilename,fname_roi);
        continue;
    end;
    roi_data = [];
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

    matfile = sprintf('%s/analysis/aparc_stats.mat',ContainerPath);
    if exist(matfile,'file')
      aparc_stats = [];
      load(matfile);
      if ~isempty(aparc_stats)
        if ~wrote_header_flag(numvols)
          nrois = length(aparc_stats);
          fprintf(fid,'"SubjectID","StudyDate"');
          for j=2:size(metaheaders,2)
              fprintf(fid,',"%s"',metaheaders{j});
          end
          for k=1:length(aparc_stats)
              fprintf(fid,',"%s"',aparc_stats(k).roiname);
          end;
          fprintf(fid,'\n');
          wrote_header_flag(numvols)=1;
        end;
        fprintf(fid,'"%s","%s"',SubjID,StudyDate);
        for j=2:size(metaheaders,2)
          if ~isempty(S_ind) % copy subject metadata
              fprintf(fid,',"%s"',num2str(metadata{S_ind,j}));
          else % blank values if no metadata
              fprintf(fid,',');
          end
        end
        if length(aparc_stats)~=nrois
          fprintf('%s: WARNING: number of ROIs (%d) for %s does not match first (%d)\n',...
            mfilename,length(aparc_stats),ContainerPath,nrois);
          for k=1:nrois
            fprintf(fid,',NaN');
          end;
        else
          for k=1:nrois
            fprintf(fid,',%0.6f',...
              aparc_stats(k).thickavg);
          end;
        end;
        fprintf(fid,'\n');
      end
    else
      fprintf('%s: WARNING: aparc_stats.mat not found for %s\n',...
        mfilename,dirlist(i).name);
    end
    fclose(fid);
  end
end

