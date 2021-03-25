function [RAWContainerDir,errcode] = MMIL_Unpack_NonDicoms(...
  origRootDir,origSourceDir,RAWContainerRootDir,VisitID,fname_info,forceflag)
%function [RAWContainerDir,errcode] = MMIL_Unpack_NonDicoms(...
% origRootDir,origSourceDir,RAWContainerRootDir,VisitID,fname_info,forceflag)
%
% created:  05/28/09 by Alain Koyama
% Rcnt Mod: 09/09/12 by Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

%% todo: no unpackstatus file

if ~exist('forceflag','var'), forceflag = false; end;

errcode = 0;

origPath = sprintf('%s/%s',origRootDir,origSourceDir);

statfname = sprintf('%s/.unpackstatus',origPath);
errfname = sprintf('%s/.errors',origPath);

% read unpackstatus to get output container
%   any problems, do unpacking again
if exist(statfname,'file') && ~forceflag
  fid_stat = fopen(statfname,'rt');
  if fid_stat~=-1
    tmp = fgetl(fid_stat);
    fclose(fid_stat);
    if ischar(tmp)
      ContainerOutPath = tmp;
      % ignore rootdir provided
      try
        [tmp_path,tmp_fstem,tmp_ext] = fileparts(ContainerOutPath);
        ContainerOutDir = [tmp_fstem,tmp_ext];
        ContainerOutPath = [RAWContainerRootDir '/' ContainerOutDir];
        if exist(ContainerOutPath,'dir'), return; end;
      catch
      end;
    end;    
  end;
end;

if exist(statfname,'file'), delete(statfname); end
if exist(errfname,'file'), delete(errfname); end
d = recursive_dir(sprintf('%s',origPath));
if isempty(d)
  fid_err = fopen(errfname,'a');
  fprintf(fid_err,'%s: no files in %s\n',mfilename,origPath);
  fclose(fid_err);
  fprintf('%s: ERROR: no files found in %s\n',mfilename,origPath);
  errcode = 1;
  return;
end;

% load info from input csv file
tmpinfo = MMIL_Read_StudyInfo(fname_info);

% create ContainerInfo
ContainerInfo.nondicomflag = 1;
ContainerInfo.VisitID = VisitID;
ContainerInfo.StudyDate = num2str(tmpinfo(1).SeriesDate);
ContainerInfo.StudyTime = num2str(tmpinfo(1).SeriesTime);
ContainerInfo.Manufacturer = tmpinfo(1).Manufacturer;
ContainerInfo.MagneticFieldStrength = tmpinfo(1).MagneticFieldStrength;
origdirs = {tmpinfo.SeriesDirectory};

% set up RAW dir
ContainerUID = '1';
RAWContainerDir = sprintf('MRIRAW_%s_%s.%s_%s',VisitID,ContainerInfo.StudyDate,ContainerInfo.StudyTime,ContainerUID);
ContainerOutPath = sprintf('%s/%s',RAWContainerRootDir,RAWContainerDir);
if ~exist(ContainerOutPath,'dir'), mkdir(ContainerOutPath); end
ContainerInfoFile = sprintf('%s/ContainerInfo.mat',ContainerOutPath);
if exist(ContainerInfoFile,'file') & ~forceflag, return; end;

% copy each series to seperate RAW subdir
dirlistorig = dir(origPath);
dirlistorig = dirlistorig(cellfun('isempty', regexp({dirlistorig.name},'^\.'))); % get rid of . and .. in dirlist
numseries = 0;
for j=1:length(dirlistorig) 
  if ~dirlistorig(j).isdir, continue; end;
  numseries = numseries + 1;
  origdir = dirlistorig(j).name;
  S_ind=find(strcmp(origdir,origdirs));
  if isempty(S_ind), continue; end;
  RAWsubdir = sprintf('ser%04d',numseries);
  ContainerInfo.SeriesInfo(numseries).SeriesDirPath = RAWsubdir;
  ContainerInfo.SeriesInfo(numseries).SeriesType = tmpinfo(S_ind).SeriesType;
  [gradwarpinfo,errmsg] = ctx_get_gradwarpinfo(tmpinfo);
  if ~isempty(errmsg)
    fprintf('%s: WARNING: %s',mfilename,errmsg);
  end
  ContainerInfo.SeriesInfo(numseries).gradwarpinfo = gradwarpinfo;
  cmd = sprintf('cp -r %s/%s %s/%s',origPath,dirlistorig(j).name,ContainerOutPath,RAWsubdir);
  [status, result] = system(cmd);
  if status
    fprintf('%s: Error - %s',mfilename,result);
    errcode = 1;
    return;
  end
end

if exist('ContainerOutPath','var')
  save(ContainerInfoFile,'ContainerInfo');
  fid_stat = fopen(statfname,'w');
  if fid_stat == -1
    fprintf('%s: ERROR: failed to write file %s... check permissions\n',...
      mfilename,statfname);
    return;
  end;
  fprintf(fid_stat,'unpacked by %s at %f (%s) to\n%s\n',...
    mfilename,now,datestr(now),ContainerOutPath); 
  fclose(fid_stat);
else
  fprintf('%s: WARNING: unpacking failed\n',mfilename);
end


