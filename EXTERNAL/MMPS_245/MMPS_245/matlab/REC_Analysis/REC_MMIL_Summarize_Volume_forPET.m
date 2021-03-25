function REC_MMIL_Summarize_Volume_forPET(ProjID)
%function REC_MMIL_Summarize_Volume_forPET(ProjID)
%
% Required Parameters:
%   ProjID: project ID to run (e.g. 'REC_TEST')
%
% Created:  10/31/08 by Alain Koyama
% Rcnt Mod: 03/24/10 by Don Hagler
% Last Mod: 09/09/12 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
if isempty(ProjID), error('empty ProjID'); end;

RootDirs = REC_MMIL_RootDirs(ProjID);

roilist = [1:28,40:60,77:79,20001:20003];
exclude_roilist = [1,6,9,19:23,27,40,45,48,55,56,59];
roilist = setdiff(roilist,exclude_roilist);

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
for i = 1:length(dirlist_pet)
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
      outfile{numvols} = sprintf('%s/REC_%s_%s_SubCortVolume_forPET.csv',outdir,ProjID,char(volnames{numvols}));
      fid=fopen(outfile{numvols},'wt');
      fclose(fid);
    end
    wrote_header_flag=zeros(1,length(volnames));
  end
  createinfo = 1;

  % sub-cortical volume
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
    for j=1:length(dirlist_fs) %find matching studies
      if regexp(dirlist_fs(j).name,SubjID)
        studymatch{end+1} = dirlist_fs(j).name;
      end
    end
    if ~isempty(studymatch)
      studymatch = sort(studymatch);
      ContainerPath = sprintf('%s/%s',RootDirs.fsurf,studymatch{1});
    else
      fprintf('%s: WARNING: no FS recon found... skipping %s\n',mfilename,SubjID);
      continue;
    end

    matfile = sprintf('%s/stats/segstats.mat',ContainerPath);
    if exist(matfile,'file')
      aseg_stats = [];
      load(matfile);
      if ~isempty(aseg_stats)
        % reduce to subset of ROIs
        roicodes = cell2mat({aseg_stats.roicode});
        i_roicodes = find(ismember(roicodes,roilist));
        aseg_stats = aseg_stats(i_roicodes);
        if ~wrote_header_flag
          nrois = length(aseg_stats);
          fprintf(fid,'"SubjectID","StudyDate"');
          for j=2:size(metaheaders,2)
            fprintf(fid,',"%s"',metaheaders{j});
          end
          for k=1:nrois
            fprintf(fid,',"%s"',aseg_stats(k).roiname);
          end;
          fprintf(fid,'\n');
          wrote_header_flag=1;
        end;
        fprintf(fid,'"%s","%s"',SubjID,StudyDate);
        for j=2:size(metaheaders,2)
          if ~isempty(S_ind) % copy subject metadata
            fprintf(fid,',"%s"',num2str(metadata{S_ind,j}));
          else % blank values if no metadata
            fprintf(fid,',');
          end
        end
        if length(aseg_stats)~=nrois
          fprintf('%s: WARNING: number of ROIs (%d) for %s does not match first (%d)\n',...
            mfilename,length(aseg_stats),ContainerPath,nrois);
          for k=1:nrois
            fprintf(fid,',NaN');
          end
        else
          for k=1:nrois
            fprintf(fid,',%0.6f',...
              aseg_stats(k).volume);
          end
        end
        fprintf(fid,'\n');
      end
    else
      fprintf('%s: Warning: segstats.mat not found for %s\n',mfilename,dirlist_fs(i).name);
    end
  end
  fclose(fid);
end
