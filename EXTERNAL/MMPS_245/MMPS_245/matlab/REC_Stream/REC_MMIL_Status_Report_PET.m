function REC_MMIL_Status_Report_PET(ProjID)
%function REC_MMIL_Status_Report_PET(ProjID)
%
%
% Created:  09/09/08 by Alain Koyama
% Rcnt Mod: 03/24/10 by Don Hagler
% Last Mod: 09/09/12 by Don Hagler
%

if ~exist('ProjID','var'), ProjID = []; end;
if isempty(ProjID), error('Empty ProjID'); end

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);

if isempty(RootDirs.orig_pet)
  fprintf('%s: Project %s does not appear to have any PET data\n',mfilename,ProjID);
  continue;
end
metadatadir = sprintf('%s/MetaData/%s',RootDirs.home,ProjID);
if ~exist(metadatadir,'dir'), mkdir(metadatadir); end

qcfile = sprintf('%s/REC_%s_PET_QC.csv',metadatadir,ProjID);
if ~exist(qcfile,'file')
    fprintf('%s: WARNING: unable to find %s\n',mfilename,qcfile);
    fprintf('%s: All QC values will be 0\n',mfilename,qcfile);
end

outfile = sprintf('%s/REC_%s_Status_Report_PET.csv',metadatadir,ProjID);
fout = fopen(outfile,'wt');
if fout==-1
  fprintf('%s: ERROR: unable to create output file\n',mfilename);
  return;
end;

fprintf(fout,'"SubjID","StudyDate","VisitNumber","PETRAW Container",');
fprintf(fout,'"vols.mat","PETPROC container","PET_reg.mgh",');
fprintf(fout,'"FS Container","recon complete",');
fprintf(fout,'"roi stats","painted",');
fprintf(fout,'"aseg QC","PETreg QC"');

fprintf(fout,'\n');

dirlist_orig = dir(RootDirs.orig_pet); dirlist_orig = dirlist_orig(3:end);
dirlist_raw=dir(sprintf('%s/PETRAW*',RootDirs.raw_pet));
dirlist_proc=dir(sprintf('%s/PETPROC*',RootDirs.proc_pet));
dirlist_fs=dir(sprintf('%s/FSURF*',RootDirs.fsurf));

for i=1:length(dirlist_orig)   
    SubjID = dirlist_orig(i).name;
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
    asegQC = 0;
    petregQC = 0;

    % unpacked?
    for j=1:length(dirlist_raw)
        raw_dir = dirlist_raw(j).name;
        RAW_ContainerPath = sprintf('%s/%s',RootDirs.raw_pet,raw_dir);
        if dirlist_raw(j).isdir & ~isempty(RAW_ContainerPath) & regexpi(raw_dir,SubjID)
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
    for j=1:length(dirlist_proc)
        proc_dir = dirlist_proc(j).name;
        PROC_ContainerPath = sprintf('%s/%s',RootDirs.proc_pet,proc_dir);
        if dirlist_proc(j).isdir & ~isempty(PROC_ContainerPath) & regexpi(proc_dir,SubjID)
            PROC_container_flag = 1;
            % pet_reg.mgh?
            files_proc = dir(PROC_ContainerPath);
            for k=1:length(files_proc)
                curfile = files_proc(k).name;
                if regexp(curfile,'^PET_reg') & regexp(curfile,'\.mgh$')
                    petreg_flag = 1; break;
                end
            end
            % thickness painted?
            lh_paint = 0; rh_paint = 0;
            for k=1:length(files_proc)
                curfile = files_proc(k).name;
                if ~lh_paint & regexp(curfile,'^PET_reg') & regexp(curfile,'-paint-lh\.mgh$')
                    lh_paint = 1;
                end
                if ~rh_paint & regexp(curfile,'^PET_reg') & regexp(curfile,'-paint-rh\.mgh$')
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
            fname = sprintf('%s/stats/segstats.mat',FS_ContainerPath);
            if exist(fname,'file')
                roistats_flag = 1;
            end
            break;
        end
    end

    % qc
    petqc = mmil_readtext(qcfile);
    petqc_SubjID = petqc(:,1);
    petqc_petregQC = petqc(:,2);
    petqc_asegQC = petqc(:,3);
    ind = strcmp(petqc_SubjID,SubjID);
    if ~unique(ind)
        fprintf('%s: WARNING - subject %s not found in %s\n',mfilename,SubjID,qcfile);
    else
        switch upper(petqc_petregQC{ind})
            case 'BAD', petregQC = -1;
            case 'AVERAGE', petregQC = 1;
            case 'GOOD', petregQC = 1;
        end
        switch upper(petqc_asegQC{ind})
            case 'BAD', asegQC = -1;
            case 'AVERAGE', asegQC = 1;
            case 'GOOD', asegQC = 1;
        end
    end

    % write to log
    fprintf(fout,'"%s","%s",%d,',...
        SubjID,StudyDate,VisitNumber);
    fprintf(fout,'%d,%d,%d,%d,',...
        RAW_container_flag,vols_flag,...
        PROC_container_flag,petreg_flag);
    fprintf(fout,'%d,%d,%d,%d,%d,%d,%d',...
        FS_container_flag,recon_flag,roistats_flag,painted,...
        asegQC,petregQC);
    fprintf(fout,'\n');
end
fclose(fout);
