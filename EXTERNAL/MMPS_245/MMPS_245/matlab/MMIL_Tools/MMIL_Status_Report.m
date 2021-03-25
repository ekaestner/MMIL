function MMIL_Status_Report(RootDirs,varargin)
%function MMIL_Status_Report(RootDirs,varargin)
%
% Created:  04/14/09 by Alain Koyama
% Last Mod: 02/12/13 Don Hagler
%

% todo: do DTI, BOLD, diff QC's (raw, aseg, surf)

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'ProjID','',[],...
  'analysis_outdir','analysis','',... % to account for analyses done w/ older and newer Analyze scripts
  'sourcedir','orig',{'orig','raw'},... % if orig were deleted
  'pet_flag',0,[0 1],...
  'dti_flag',0,[0 1],...
  'bold_flag',0,[0 1],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

homedir = RootDirs.homedir;

outdir = sprintf('%s/MetaData/%s',homedir,parms.ProjID);
outfileMRI = sprintf('%s/MetaData/%s/%s_Status_Report_MRI.csv',homedir,parms.ProjID,parms.ProjID);
outfilePET = sprintf('%s/MetaData/%s/%s_Status_Report_PET.csv',homedir,parms.ProjID,parms.ProjID);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

fout = fopen(outfileMRI,'wt');
if fout==-1
  fprintf('%s: ERROR: unable to create output file\n',mfilename);
  return;
end;

fprintf(fout,'"SubjID","StudyDate","VisitNumber","Temp",');
fprintf(fout,'"RAW Container","classified","PROC Container","processed",');
fprintf(fout,'"FS Container","recon complete","roi stats",');
fprintf(fout,'\n');

dirlist_orig = dir(RootDirs.orig); 
n = regexp({dirlist_orig.name},'^\.'); % get rid of . and .. in dirlist
n = cellfun('isempty', n);
dirlist_orig = dirlist_orig(n);

dirlist_raw=dir(sprintf('%s/MRIRAW*',RootDirs.raw));
dirlist_proc=dir(sprintf('%s/MRIPROC*',RootDirs.proc));
dirlist_fs=dir(sprintf('%s/FSURF*',RootDirs.fsurf));

switch parms.sourcedir
    case 'orig', dirlist_source = dirlist_orig;
    case 'raw', dirlist_source = dirlist_raw;
end

for i=1:length(dirlist_source)
    switch parms.sourcedir
        case 'orig'
            SubjID = dirlist_source(i).name;
        case 'raw'
            SubjID = char(regexp(dirlist_source(i).name,'(?<=MRIRAW_).+(?=_\d{8})','match'));
    end
    StudyDate = '';
    visit = char(regexp(SubjID,'(?<=v)\d{1}$','match'));
    if ~isempty(visit)
        VisitNumber = str2double(visit);
    else
        VisitNumber = 1;
    end

    downloaded_flag = 1;
    RAW_container_flag = 0;
    classified_flag = 0;
    PROC_container_flag = 0;
    processed_flag = 0;
    FS_container_flag = 0;
    recon_flag = 0;
    roistats_flag = 0;
    
    % unpacked?
    n = regexp({dirlist_raw.name},['MRIRAW_' SubjID '_\d{8}.+_1$'],'match');
    n = char([n{:}]);
    if ~isempty(n)
        RAW_container_flag = 1;
        StudyDate = char(regexp(n,'(?<=_)\d{8}(?=\.)','match'));
        RAW_ContainerPath = sprintf('%s/%s',RootDirs.raw,n);
        % classified?
        fname = sprintf('%s/ContainerInfo.mat',RAW_ContainerPath);
        if exist(fname,'file')
            try
                load(fname);
                if isfield(ContainerInfo,'ClassificationCompleted')
                    classified_flag = 1;
                end
            catch
                sprintf('%s: WARNING: failed to load %s\n',mfilename,fname);
                classified_flag = 0;
            end
        end
    end
    
    % processed?
    n = regexp({dirlist_proc.name},['MRIPROC_' SubjID '_\d{8}.+_1$'],'match');
    n = char([n{:}]);
    PROC_ContainerPath = sprintf('%s/%s',RootDirs.proc,n);
    if ~isempty(n)
        PROC_container_flag = 1;
        %fprintf('%s  %s\n',SubjID,proc_dir);
        % MPR_res or hiFA_res exists?
        fname1 = sprintf('%s/MPR_res%s',PROC_ContainerPath,parms.ext);
        fname2 = sprintf('%s/hiFA_res%s',PROC_ContainerPath,parms.ext);
        if exist(fname1,'file') | exist(fname2,'file'), processed_flag = 1; end;
    end

    % recon'd?
    n = regexp({dirlist_fs.name},['FSURF_' SubjID '_\d{8}.+_1$'],'match');
    n = char([n{:}]);
    if ~isempty(n)
        FS_ContainerPath = sprintf('%s/%s',RootDirs.fsurf,n);
        FS_container_flag = 1;
        % recon status
        [status,message] = MMIL_Get_FSReconStatus(FS_ContainerPath);
        if status==2 | status==5
            recon_flag = 1;
        end
        % stats analyzed?
        fname_asegstats = sprintf('%s/%s/aseg_stats.mat',FS_ContainerPath,parms.analysis_outdir);
        fname_aparcstats = sprintf('%s/%s/aparc_stats.mat',FS_ContainerPath,parms.analysis_outdir);
        if exist(fname_asegstats,'file') & exist(fname_aparcstats,'file')
          roistats_flag = 1;
        end
    end
    
    % write to log
    fprintf(fout,'"%s","%s",%d,',SubjID,StudyDate,VisitNumber);
    fprintf(fout,'%d,%d,%d,%d,%d,',downloaded_flag,RAW_container_flag,...
        classified_flag,PROC_container_flag,processed_flag);
    fprintf(fout,'%d,%d,%d',...
      FS_container_flag,recon_flag,roistats_flag);
    fprintf(fout,'\n');
end

fclose(fout);

%%%%%%%%%%%%%%%%%%%%%%%% do PET

if parms.pet_flag
  fout = fopen(outfilePET,'wt');
  if fout==-1
    fprintf('%s: ERROR: unable to create output file\n',mfilename);
    return;
  end;

  fprintf(fout,'"SubjID","StudyDate","VisitNumber","PETRAW Container",');
  fprintf(fout,'"vols.mat","PETPROC container","PET_reg%s",',parms.ext);
  fprintf(fout,'"FS Container","recon complete",');
  fprintf(fout,'"roi stats","painted"');
  fprintf(fout,'\n');

  dirlist_petorig=dir(sprintf('%s/*',RootDirs.orig_pet));
  n = regexp({dirlist_petorig.name},'^\.'); % get rid of . and .. in dirlist
  n = cellfun('isempty', n);
  dirlist_petorig = dirlist_petorig(n);

  dirlist_petraw=dir(sprintf('%s/*',RootDirs.raw_pet));
  dirlist_petproc=dir(sprintf('%s/*',RootDirs.proc_pet));


  for i=1:length(dirlist_petorig)
    SubjID = dirlist_petorig(i).name;
    StudyDate = '';
    visit = char(regexp(SubjID,'(?<=v)\d{1}$','match'));
    if ~isempty(visit)
      VisitNumber = str2double(visit);
    else
      VisitNumber = 1;
    end

    RAW_container_flag = 0;
    vols_flag = 0;
    PROC_container_flag = 0;
    petreg_flag = 0;
    FS_container_flag = 0;
    recon_flag = 0;
    roistats_flag = 0;
    painted = 0;

    % unpacked?
    for j=1:length(dirlist_petraw)
      raw_dir = dirlist_petraw(j).name;
      RAW_ContainerPath = sprintf('%s/%s',RootDirs.raw_pet,raw_dir);
      if dirlist_petraw(j).isdir & ~isempty(RAW_ContainerPath) & regexpi(raw_dir,SubjID)
        RAW_container_flag = 1;
        StudyDate = char(regexp(raw_dir,'(?<=_)\d{8}(?=\.)','match'));
        % vols.mat?
        files_raw = dir(RAW_ContainerPath);
        for k=1:length(files_raw)
          curfile = files_raw(k).name;
          if regexp(curfile,'^vols') & regexp(curfile,'\.mat$')
            vols_flag = 1;
            break
          end
        end
        break;
      end
    end

    % processed?
    for j=1:length(dirlist_petproc)
      proc_dir = dirlist_petproc(j).name;
      PROC_ContainerPath = sprintf('%s/%s',RootDirs.proc_pet,proc_dir);
      if dirlist_petproc(j).isdir & ~isempty(PROC_ContainerPath) & regexpi(proc_dir,SubjID)
        PROC_container_flag = 1;
        % pet_reg.mgh?
        files_proc = dir(PROC_ContainerPath);
        for k=1:length(files_proc)
          curfile = files_proc(k).name;
          if regexp(curfile,'^PET_reg') & regexp(curfile,sprintf('\\%s$',parms.ext))
            petreg_flag = 1; break;
          end
        end
        % thickness painted?
        lh_paint = 0; rh_paint = 0;
        for k=1:length(files_proc)
          curfile = files_proc(k).name;
          if ~lh_paint & regexp(curfile,'^PET_reg') & regexp(curfile,sprintf('-paint-lh\\%s$',parms.ext))
            lh_paint = 1;
          end
          if ~rh_paint & regexp(curfile,'^PET_reg') & regexp(curfile,sprintf('-paint-rh\\%s$',parms.ext))
            rh_paint = 1;
          end
          if lh_paint && rh_paint
            painted = 1;
            break;
          end
        end
        break;
      end
    end

    % recon'd?
    for j=1:length(dirlist_fs)
      FS_dir = dirlist_fs(j).name;
      FS_ContainerPath = sprintf('%s/%s',RootDirs.fsurf,FS_dir);
      if dirlist_fs(j).isdir & ~isempty(FS_ContainerPath) & regexpi(FS_dir,SubjID)
        FS_container_flag = 1;
        % recon status
        [status,message] = MMIL_Get_FSReconStatus(FS_ContainerPath);
        if status==2 | status==5
          recon_flag = 1;
        end
        % stats analyzed?
        fname_asegstats = sprintf('%s/%s/aseg_stats.mat',...
          FS_ContainerPath,parms.analysis_outdir);
        fname_aparcstats = sprintf('%s/%s/aparc_stats.mat',...
          FS_ContainerPath,parms.analysis_outdir);
        if exist(fname_asegstats,'file') & exist(fname_aparcstats,'file')
          roistats_flag = 1;
        end
        break;
      end
    end

    % write to log
    fprintf(fout,'"%s","%s",%d,',...
      SubjID,StudyDate,VisitNumber);
    fprintf(fout,'%d,%d,%d,%d,',...
      RAW_container_flag,vols_flag,...
      PROC_container_flag,petreg_flag);
    fprintf(fout,'%d,%d,%d,%d',...
      FS_container_flag,recon_flag,roistats_flag,painted);
    fprintf(fout,'\n');
  end


  fclose(fout);
end




