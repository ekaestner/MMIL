function REC_MMIL_Status_Report_MRI(ProjID,varargin)
%function REC_MMIL_Status_Report_MRI(ProjID,varargin)
%
% Created:  05/18/09 by Alain Koyama
% Last Mod: 11/06/12 by Don Hagler
%

% todo: for dti analysis, check for mat's in fiber_masks_from_atlas
% todo: add a Status_Report_Summary or something, or add this to main file

if ~exist('ProjID','var'), ProjID = []; end;
if isempty(ProjID)
  error('empty ProjID');
end

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);
homedir = RootDirs.home;

parms = mmil_args2parms(varargin, { ...
  'infix','ecc_reg_mc_B0uw_gruw_iso',[],...
  'atlas_flag',[2],[],...
  'weighted_avg_flag',1,[0 1],...
  'FA_thresh',[0,0.15],[],...
  'prob_thresh',[0,0.07,0.08],[],...
  'fibers',[101:110,115:123,133:138,141:150],[],...
  'atl_tensor_smooth_sigma',5,[0,100],...
  'measlist',{'FA','ADC'},[],...
  'first_only_flag',01,[],...
  'erode_flag',1,[0 1],...
  'ext','.mgz',{'.mgh','.mgz'},...
});
    
outdir = sprintf('%s/MetaData/%s',homedir,ProjID);
outfile = sprintf('%s/MetaData/%s/REC_%s_Status_Report_MRI.csv',...
  homedir,ProjID,ProjID);
mmil_mkdir(outdir);

fout = fopen(outfile,'wt');
if fout==-1
  fprintf('%s: ERROR: unable to create output file\n',mfilename);
  return;
end;

fprintf(fout,'"SubjID","StudyDate","VisitNumber","Temp",');
fprintf(fout,'"RAW Container","classified","PROC Container","processed",');
fprintf(fout,'"FS Container","recon complete","roi stats","FS Container (v4)","recon complete (v4)","roi stats (v4)",');
fprintf(fout,'"DTI Data","DTI Fiber Analysis","DTI Aseg Analysis","FiberTracks"');
fprintf(fout,'\n');

dirlist_tmp = dir(ProjInfo.TMP_MRI_RootDir); dirlist_tmp = dirlist_tmp(3:end);
dirlist_raw=dir(sprintf('%s/MRIRAW*',RootDirs.raw));
dirlist_proc=dir(sprintf('%s/MRIPROC*',RootDirs.proc));
dirlist_fs=dir(sprintf('%s/FSURF*',RootDirs.fsurf));

for i=1:length(dirlist_tmp)   
  SubjID = dirlist_tmp(i).name;
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
  dti_container_flag = 0;
  dti_fiber = 0;
  dti_aseg = 0;
  fibers_flag = 0;

  dtiflag = ProjInfo.DTIflag;

  % unpacked?
  n = regexp({dirlist_raw.name},['MRIRAW_' SubjID '_\d{8}.+'],'match');
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
  n = regexp({dirlist_proc.name},['MRIPROC_' SubjID '_\d{8}.+'],'match');
  n = char([n{:}]);
  PROC_ContainerPath = sprintf('%s/%s',RootDirs.proc,n);
  if ~isempty(n)
    PROC_container_flag = 1;
    % MPR_res or hiFA_res exists?
    fname1 = sprintf('%s/MPR_res%s',PROC_ContainerPath,parms.ext);
    fname2 = sprintf('%s/hiFA_res%s',PROC_ContainerPath,parms.ext);
    if exist(fname1,'file') | exist(fname2,'file'), processed_flag = 1; end;
    % dti?
    cfile = fullfile(RootDirs.proc,n,'ContainerInfo.mat');
    if exist(cfile,'file')
        load(cfile);
    end
    if isfield(ContainerInfo,'DTI_cntr') && ContainerInfo.DTI_cntr >= 1, dti_container_flag = 1; end;

    % dti analysis?
    dti_dir = sprintf('%s/%s/DTanalysis',RootDirs.proc,n);
    dti_fiber = 1;
    dti_aseg = 1;
    if ~dti_container_flag || ~exist(dti_dir,'dir')
      dti_fiber = 0; dti_aseg = 0;
    else % track fiber analysis mat files generated in DTanalysis
      dti_files = sprintf('%s/DTI*_%s_*.mat',dti_dir,parms.infix);
      dti_aseg_files = sprintf('%s/DTI*_%s_*_aseg*.mat',dti_dir,parms.infix);
      if isempty(dir(dti_files))
        fprintf('%s: File %s not found for %s\n',mfilename,dti_files,SubjID);
        dti_fiber = 0;
      end
      if isempty(dir(dti_aseg_files))
        fprintf('%s: File %s not found for %s\n',mfilename,dti_aseg_files,SubjID);
        dti_aseg = 0;
      end
    end

    % fiber tracks?
    fiberdir = sprintf('%s/fiber_masks_from_atlas',PROC_ContainerPath);
    fiber_files = {};
    fibers_flag = 1;
    if (~dti_container_flag | ~exist(fiberdir,'dir')) & ~max(ismember([1 3],parms.atlas_flag))
      fibers_flag = 0;
    else
      for f=parms.fibers
        for atlas_flag = parms.atlas_flag
          switch atlas_flag
            case 0 % manual
              if weighted_avg_flag
                  fiberstem = sprintf('%s/fiber_%02d_count',fiberdir,f);
              else
                  fiberstem = sprintf('%s/fiber_%02d_mask',fiberdir,f);
              end;
              fiberext = '.mgh';
            case 1 % loc only, count atlas
              fiberstem = sprintf('%s/fiber_%02d_fiber_from_countatlas',...
                  fiberdir,f);
              fiberext = '.mgz';
            case 2 % loc+dir, count atlas
              fiberstem = sprintf('%s/fiber_%02d_sm%0.2f_probV0_countatlas',...
                  fiberdir,f,parms.atl_tensor_smooth_sigma);
              fiberext = '.mgz';
            case 3 % loc only, mask atlas
              fiberstem = sprintf('%s/fiber_%02d_fiber_from_atlas',...
                  fiberdir,f);
              fiberext = '.mgz';
            case 4 % loc+dir, mask atlas
              fiberstem = sprintf('%s/fiber_%02d_sm%0.2f_probV0',...
                  fiberdir,f,parms.atl_tensor_smooth_sigma);
              fiberext = '.mgz';
          end;
          maskstem = fiberstem;
          fname_mask = [maskstem fiberext];
          fiber_files{end+1} = fname_mask;
        end
        if parms.first_only_flag
          fname_tensor = sprintf('%s/fiber_%02d_tensorV0_from_atlas_sm%0.2f_norm.mgz',...
              fiberdir,f,parms.atl_tensor_smooth_sigma);
          fname_prob_dir = sprintf('%s/fiber_%02d_sm%0.2f_prob_dirV0.mgz',...
              fiberdir,f,parms.atl_tensor_smooth_sigma);
        else
          fname_tensor = sprintf('%s/fiber_%02d_tensor_from_atlas_sm%0.2f_norm.mgz',...
              fiberdir,f,parms.atl_tensor_smooth_sigma);
          fname_prob_dir = sprintf('%s/fiber_%02d_sm%0.2f_prob_dir.mgz',...
              fiberdir,f,parms.atl_tensor_smooth_sigma);
        end
        fiber_files{end+1} = fname_tensor;
        if dtiflag, fiber_files{end+1} = fname_prob_dir; end; % if fiber tracks from loc data
      end
      for k=1:length(fiber_files)
        if ~exist(fiber_files{k},'file')
          % fprintf('%s: File %s not found for %s\n',mfilename,fiber_files{k},SubjID);
          fibers_flag = 0;
          break;
        end
      end
    end
  end

  % recon
  n = regexp({dirlist_fs.name},['FSURF_' SubjID '_\d{8}.+'],'match');
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
    fname_segstats = sprintf('%s/stats/segstats.mat',FS_ContainerPath);
    fname_asegstats = sprintf('%s/stats/aseg_stats.mat',FS_ContainerPath);
    fname_aparcstats = sprintf('%s/stats/aparc_stats.mat',FS_ContainerPath);
    if exist(fname_segstats,'file') || ...
        (exist(fname_asegstats,'file') & exist(fname_aparcstats,'file'))
      roistats_flag = 1;
    end
  end

  % write to log
  fprintf(fout,'"%s","%s",%d,',SubjID,StudyDate,VisitNumber);
  fprintf(fout,'%d,%d,%d,%d,%d,',...
    downloaded_flag,RAW_container_flag,classified_flag,...
    PROC_container_flag,processed_flag);
  fprintf(fout,'%d,%d,%d,%d,%d,%d,%d',...
    FS_container_flag,recon_flag,roistats_flag,dti_container_flag,...
    dti_fiber,dti_aseg,fibers_flag);
  fprintf(fout,'\n');
end
fclose(fout);

