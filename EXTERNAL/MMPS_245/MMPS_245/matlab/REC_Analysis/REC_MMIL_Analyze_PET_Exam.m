function REC_MMIL_Analyze_PET_Exam(ContainerRootDir,ContainerDir,normflag,forceflag)
%function REC_MMIL_Analyze_PET_Exam(ContainerRootDir,ContainerDir,normflag,forceflag)
%
% Rcnt Mod: 08/15/12 by Don Hagler
% Last Mod: 11/03/12 by Don Hagler
%

if ~exist('normflag','var'), normflag = 1; end
if ~exist('forceflag','var'), forceflag = 0; end

fprintf('%s(''%s'',''%s'',%d,%d)\n',mfilename,ContainerRootDir,ContainerDir,normflag,forceflag);

homedir = '/home/mmilrec';
volnames = {}; % if multiple series (e.g. GMR and Patlak for SAX_OCD)

hemilist = {'lh' 'rh'};
smoothsteps_list = [12,50,112,199,311,699];
% approx FWHM (mm) = 4, 8, 12, 16, 20, 30
% slope of FWHM vs. sqrt(N) is ~1.13 for fsaverage
% (FWHM = full-width-half-max smoothing kernel, N = num smoothing steps)
norm_mask_name = 'pons_mask';
norm_mask_file = sprintf('%s/MetaData/REC/AtlasMasks/%s.mgh',...
  homedir,norm_mask_name);
normroi = 16; % brain-stem
%roilist = [0:1099,2000:2099];
roilist = [2 3 7 8 10 11 12 13 16 17 18 26 28 41 42 46 47 49 50 51 52 53 54 58 60 1001 1002 1003 1004 1005 1006 1007 1008 1009 1010 1011 1012 1013 1014 1015 1016 1017 1018 1019 1020 1021 1022 1023 1024 1025 1026 1027 1028 1029 1030 1031 1032 1033 1034 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030 2031 2032 2033 2034];

ContainerPath = sprintf('%s/%s',ContainerRootDir,ContainerDir);

if ~exist(ContainerPath,'file')
  fprintf('%s: %s not found...skipping\n',mfilename,ContainerPath);
  return;
end;

% get necessary files
PET_files = dir(sprintf('%s/PET_reg*txt',ContainerPath));
if isempty(PET_files)
  fprintf('%s: ERROR: PET reg file not found\n',mfilename);
  return;
else % gather all PET series
  for i=1:length(PET_files)
    volnames{end+1} = char(regexp(PET_files(i).name,'(?<=PET_reg_).+(?=\.txt)','match'));
  end
end

switch normflag
    case 0, normflags = 0;
    case 1, normflags = 1;
    case 2, normflags = [0 1];
end
for nn=1:length(normflags)
    normflag = normflags(nn);
for numvols=1:length(volnames)
  funcname = sprintf('%s/PET_reg_%s.mgh',ContainerPath,volnames{numvols});
  fprintf('%s: Start analysis for %s, normflag=%d\n',mfilename,funcname,normflag);
  regfile = sprintf('%s/PET_reg_%s.txt',ContainerPath,volnames{numvols});
  fid = fopen(regfile,'rt');
  if fid==-1
    fprintf('%s: ERROR: unable to open PET registration file %s\n',...
      mfilename,regfile);
    return;
  end;
  tline = fgetl(fid);
  fclose(fid);
  k = findstr(tline,'/');
  mridir = tline(k(1):k(end)-1);
  segname = sprintf('%s/aparc+aseg.mgz',mridir);
  if ~exist(segname,'file')
    fprintf('%s: ERROR: %s not found\n',mfilename,segname);
    return;
  end;

  subjdir = tline(k(1):k(end-2)-1);
  subjname = tline(k(end-2)+1:k(end-1)-1);

  % warp atlas roi to subj
  if normflag
    matfile = sprintf('%s/PET_%s_normroi_data.mat',ContainerPath,volnames{numvols});
    if ~exist(matfile,'file') | forceflag
      if ~exist(norm_mask_file,'file')
        fprintf('%s: ERROR: %s not found\n',mfilename,norm_mask_file);
        return;
      end;
      subj_norm_mask_file = sprintf('%s/%s.mgh',ContainerPath,norm_mask_name);
      fs_warp_vol2atlas(subjname,norm_mask_file,subj_norm_mask_file,...
        'subjdir',subjdir,'overwrite_flag',forceflag,'inverse_flag',1);
      if ~exist(subj_norm_mask_file,'file')
        fprintf('%s: ERROR: unable to create %s\n',mfilename,subj_norm_mask_file);
        return;
      end;
      % get roi avg for brainstem AND atlas roi mask
      fprintf('%s: getting ROI data...\n',mfilename);
      roi_data = mmil_aseg_roi(funcname,segname,'roilist',normroi,'fname_mask',subj_norm_mask_file);
      if isempty(roi_data)
        fprintf('%s: ERROR: aseg ROI analysis failed\n',mfilename);
        return;
      end;
      save(matfile,'roi_data');
    else
      load(matfile);
    end;
    normroi_data = roi_data;
  end
  
  if normflag
    matfile = sprintf('%s/PET_aseg_roi_data_norm_%s.mat',ContainerPath,volnames{numvols});
  else
    matfile = sprintf('%s/PET_aseg_roi_data_%s.mat',ContainerPath,volnames{numvols});
  end
  if ~exist(matfile,'file') | forceflag
    % run ROI analysis
    fprintf('%s: getting ROI data...\n',mfilename);
    roi_data = mmil_aseg_roi(funcname,segname,'roilist',roilist);
    if isempty(roi_data)
      fprintf('%s: ERROR: aseg ROI analysis failed\n',mfilename);
      return;
    end;
    if normflag
      normval = normroi_data(1).avg;
    else
      normval = [];
    end;
    if ~isempty(normval) & normval~=0
      for i=1:length(roi_data)
        roi_data(i).norm_avg = roi_data(i).avg/normval;
        roi_data(i).norm_stdv = roi_data(i).stdv/normval;
      end;
    else
      for i=1:length(roi_data)
        roi_data(i).norm_avg = 0;
        roi_data(i).norm_stdv = 0;
      end;
    end;

    % save results
    save(matfile,'roi_data');

    % convert to csv
    if normflag
        outfile = sprintf('%s/PET_aseg_roi_data_norm_%s.csv',ContainerPath,volnames{numvols});
    else
        outfile = sprintf('%s/PET_aseg_roi_data_%s.csv',ContainerPath,volnames{numvols});
    end
    fid = fopen(outfile,'wt');
    fprintf(fid,'"ROI","Avg","Stdv","normAvg","normStdv","nvals","nvals valid"\n');
    for i=1:length(roi_data)
      fprintf(fid,'"%s",%0.6f,%0.6f,%0.6f,%0.6f,%d,%d\n',...
        roi_data(i).roiname,...
        roi_data(i).avg,...
        roi_data(i).stdv,...
        roi_data(i).norm_avg,...
        roi_data(i).norm_stdv,...
        roi_data(i).nvals,...
        roi_data(i).nvals_valid);
    end;
    fclose(fid);
  else
    load(matfile);
  end;

  for i=1:length(smoothsteps_list)
    sphsmoothsteps = smoothsteps_list(i);
    % paint PET data to ico sphere surface, smooth on ico
    fs_paint(subjname,funcname,'outtype','mgh','sphere_flag',1,...
      'subjdir',subjdir,'infix','-paint',...
      'sphsmoothsteps',sphsmoothsteps,'overwrite_flag',forceflag);
  end;

  % get stats from surface ROIs
  for h=1:length(hemilist)
    hemi = hemilist{h};
    aparcname = sprintf('%s/../label/%s.aparc.annot',mridir,hemi);

    if ~exist(aparcname,'file')
      fprintf('%s: ERROR: %s not found\n',mfilename,aparcname);
      return;
    end;
    surffuncname = sprintf('%s/PET_reg_%s-paint-%s.mgh',...
      ContainerPath,volnames{numvols},hemi);
    if ~exist(surffuncname,'file')
      fprintf('%s: ERROR: %s not found\n',...
        mfilename,surffuncname);
      return;
    end;
    
    if normflag
        matfile = sprintf('%s/PET_aparc_roi_data_norm_%s_%s.mat',ContainerPath,hemi,volnames{numvols}); 
    else
        matfile = sprintf('%s/PET_aparc_roi_data_%s_%s.mat',ContainerPath,hemi,volnames{numvols});
    end
    if ~exist(matfile,'file') | forceflag
      % run ROI analysis
      fprintf('%s: getting aparc ROI data...\n',mfilename);
      roi_data = mmil_surf_roi(surffuncname,'fname_aparc',aparcname);
      if isempty(roi_data)
        fprintf('%s: ERROR: failed to get aparc ROI data\n',mfilename);
        return;
      end;

      if normflag
        normval = normroi_data(1).avg;
      else
        normval = [];
      end;

      if ~isempty(normval) & normval~=0
        for i=1:length(roi_data)
          roi_data(i).norm_avg = roi_data(i).avg/normval;
          roi_data(i).norm_stdv = roi_data(i).stdv/normval;
        end;
      else
        for i=1:length(roi_data)
          roi_data(i).norm_avg = 0;
          roi_data(i).norm_stdv = 0;
        end;
      end;

      % save results
      save(matfile,'roi_data');

      % convert to csv
      if normflag
          outfile = sprintf('%s/PET_aparc_roi_data_norm_%s.csv',ContainerPath,volnames{numvols});
      else
          outfile = sprintf('%s/PET_aparc_roi_data_%s.csv',ContainerPath,volnames{numvols});
      end
      fid = fopen(outfile,'wt');
      fprintf(fid,'"ROI","Avg","Stdv","normAvg","normStdv","nvals","nvals valid"\n');
      for i=1:length(roi_data)
        fprintf(fid,'"%s",%0.6f,%0.6f,%0.6f,%0.6f,%d,%d\n',...
          roi_data(i).roiname,...
          roi_data(i).avg,...
          roi_data(i).stdv,...
          roi_data(i).norm_avg,...
          roi_data(i).norm_stdv,...
          roi_data(i).nvals,...
          roi_data(i).nvals_valid);
      end;
      fclose(fid);

      if normflag & (isempty(normval) | normval==0)
        fprintf('%s: Error - normval empty or zero value\n',mfilename);
        return;
      end
      for i=1:length(smoothsteps_list)
        sphsmoothsteps = smoothsteps_list(i);
        % create surface map
        fname_in = sprintf('%s/PET_reg_%s-paint-mbmask-sphere-sm%d-%s.mgh',...
          ContainerPath,volnames{numvols},sphsmoothsteps,hemi);
        fname_out = sprintf('%s/PET_reg_%s-paint-mbmask-sphere-sm%d-norm-%s.mgh',...
          ContainerPath,volnames{numvols},sphsmoothsteps,hemi);
        if normflag & (~exist(fname_out,'file') | forceflag)
          if ~exist(fname_in,'file')
            fprintf('%s: ERROR: %s not found\n',...
              mfilename,fname_in);
            return;
          end;
          [vol,M] = fs_load_mgh(fname_in);
          vol = vol/normval;
          fs_save_mgh(vol,fname_out,M);
        end;
      end;
    end;
  end;

  % warp PET and aseg volumes to atlas
  fname_in = funcname;
  fname_out = sprintf('%s/PET_reg_%s_warp.mgh',ContainerPath,volnames{numvols});
  fs_warp_vol2atlas(subjname,fname_in,fname_out,...
    'subjdir',subjdir,'overwrite_flag',forceflag);

  fname_in = segname;
  fname_out = sprintf('%s/aparc+aseg_warp.mgh',mridir);
  fs_warp_vol2atlas(subjname,fname_in,fname_out,...
    'subjdir',subjdir,'overwrite_flag',forceflag);

  fname_in = sprintf('%s/nu.mgz',mridir);
  fname_out = sprintf('%s/nu_warp.mgh',mridir);
  fs_warp_vol2atlas(subjname,fname_in,fname_out,...
    'subjdir',subjdir,'overwrite_flag',forceflag);
end
end
