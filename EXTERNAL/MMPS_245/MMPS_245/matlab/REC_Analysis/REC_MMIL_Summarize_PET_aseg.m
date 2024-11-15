function REC_MMIL_Summarize_PET_aseg(ProjID,normflag)
%function REC_MMIL_Summarize_PET_aseg(ProjID,normflag)
%
% Required Parameters:
%   ProjID: project ID to run (e.g. 'REC_TEST')
%
% Optional Parameters:
%   normflag - [0|1|2] whether to normalize to pons
%    (0 = absolute values, 1 = norm to pons, 2 = both)
%    {default = 1}
%
% Early Mod: 06/12/09 by Alain Koyama
% Last Mod:  03/24/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
if isempty(ProjID), error('empty ProjID'); end;
if ~exist('normflag','var'), normflag = 1; end

RootDirs = REC_MMIL_RootDirs(ProjID);

roilist = [1:28,40:60];
exclude_roilist = [1,3,6,9,19:23,27,40,42,45,48,55,56,59];
roilist = setdiff(roilist,exclude_roilist);

outdir = sprintf('%s/MetaData/%s/ROI_Summaries',RootDirs.home,ProjID);
if ~exist(outdir,'dir'), mkdir(outdir); end;
dirlist_pet = dir(sprintf('%s/PETPROC*',RootDirs.proc_pet));
if isempty(dirlist_pet)
  fprintf('%s: Project %s has no PET Processed directory\n',mfilename,ProjID);
  continue;
end
switch normflag
  case 0, normflags = 0;
  case 1, normflags = 1;
  case 2, normflags = [0 1];
end
for nn=1:length(normflags)
  createinfo = 0;    
  normflag = normflags(nn);
  for i = 1:length(dirlist_pet) % loop thru directories
    ContainerPath = sprintf('%s/%s',RootDirs.proc_pet,dirlist_pet(i).name);

    % search for single or multiple PET series
    PET_files = dir(sprintf('%s/PET_reg*txt',ContainerPath));
    if isempty(PET_files)
      fprintf('%s: ERROR: PET reg file not found\n',mfilename);
      return;
    else
      volnames = regexp({PET_files.name},'(?<=PET_reg_).+(?=\.txt)','match');
    end

    if ~createinfo % because of possible multiple series, create output and header at first iteration
      for numvols=1:length(volnames)
        if normflag
          outfile{numvols} = sprintf('%s/REC_%s_PET_norm_%s_aseg.csv',outdir,ProjID,char(volnames{numvols}));
        else
          outfile{numvols} = sprintf('%s/REC_%s_PET_%s_aseg.csv',outdir,ProjID,char(volnames{numvols}));
        end
        fid=fopen(outfile{numvols},'wt');
        fclose(fid);
      end
      wrote_header_flag=zeros(1,length(volnames));
    end

    createinfo = 1;
    for numvols=1:length(volnames) % loop thru PET series
      fid=fopen(outfile{numvols},'at');
      SubjID = char(regexp(dirlist_pet(i).name,'(?<=PETPROC_).+(?=_\d{8})','match'));
      StudyDate = char(regexp(dirlist_pet(i).name,'(?<=PETPROC_.+)\d{8}(?=\.)','match'));
      if isempty(ContainerPath)
        fprintf('%s: Error - %s directory is empty\n',mfilename,ContainerPath);
        continue;
      end
      if normflag
        matfile = sprintf('%s/PET_aseg_roi_data_norm_%s.mat',ContainerPath,char(volnames{numvols}));
      else
        matfile = sprintf('%s/PET_aseg_roi_data_%s.mat',ContainerPath,char(volnames{numvols}));
      end
      if ~exist(matfile,'file')
        fprintf('%s: Error - %s not found\n',mfilename,matfile);
        continue;
      end
      roi_data = [];
      load(matfile);
      if ~isempty(roi_data)
        % reduce to subset of ROIs
        roicodes = cell2mat({roi_data.roicode});
        i_roicodes = find(ismember(roicodes,roilist));
        roi_data = roi_data(i_roicodes);
        if ~wrote_header_flag(numvols)
          fprintf(fid,'"SubjectID","StudyDate"');
          for k=1:length(roi_data)
            fprintf(fid,',"%s"',roi_data(k).roiname);
          end
          fprintf(fid,'\n');
          wrote_header_flag(numvols)=1;
        end;
        fprintf(fid,'"%s","%s"',SubjID,StudyDate);
        for k=1:length(roi_data)
          if normflag
            fprintf(fid,',%0.6f',roi_data(k).norm_avg);
          else
            fprintf(fid,',%0.6f',roi_data(k).avg);
          end
        end
        fprintf(fid,'\n');
      end
      fclose(fid);
    end
  end
end
