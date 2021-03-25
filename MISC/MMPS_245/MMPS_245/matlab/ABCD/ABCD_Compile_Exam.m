function ABCD_Compile_Exam(SeUIDs,VisitIDs,VisitID,root_indirs,varargin)
%function ABCD_Compile_Exam(SeUIDs,VisitIDs,VisitID,[options])
%
% Purpose: combine sets of orig and raw dirs
%          into compiled orig and raw event dirs
%
% Required Input:
%   SeUIDs: cell array of SeriesInstanceUID for each series
%   VisitIDs: cell array of VisitIDs for each series
%   VisitID: output VisitID for compiled event
%
% Optional Input: ('key', value pairs)
%   'ProjID': project ID string
%     if supplied, used to update ContainerInfo
%     {default = []}
%   'root_indir': input directory containing individual orig and raw series dirs
%     {default = '/space/syn05/1/data/MMILDB/DAL_ABCD_QC'}
%   'root_outdir': output directory to be filled with orig and raw event dirs
%     {default = '/space/syn05/1/data/MMILDB/DAL_ABCD'}
%   'linkflag': [0|1] create symbolic links to directories instead of copying
%     {default = 0}
%   'follow_links_flag': [0|1] when copying directories, turn links into copies
%     ignored if linkflag = 1
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}     
%
% Created:  12/03/16 by Don Hagler
% Last Mod: 01/24/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: option to use mmil_copy_dicoms (see mmil_unpack_dicoms)
%%         instead of compile_SeriesInfo and compile_raw_dirs
%%       reason would be if we want to follow links for orig 
%%         (e.g. if we want to archive QC)
%%         and then create new links for raw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(SeUIDs,VisitIDs,VisitID,root_indirs,varargin);

% locate orig and raw dirs corresponding to VisitIDs
parms = locate_indirs(parms);

% create output orig directory
parms = create_orig_outdir(parms);

% copy contents of orig dirs into a single orig dir
parms = compile_orig_dirs(parms);

% create output raw directory
parms = create_raw_outdir(parms);

% create concatenated OrigSeriesInfo in raw dir
compile_OrigSeriesInfo(parms);

% create concatenated SeriesInfo in raw dir
compile_SeriesInfo(parms);

% create ContainerInfo with concatenated SeriesInfo in raw dir
compile_ContainerInfo(parms);

% summarize SeriesInfo in csv file
summarize_SeriesInfo(parms)

% copy contents of input raw dirs into a single output raw dir
compile_raw_dirs(parms);

% create unpacking.txt file
write_progress_file(parms,'unpacking');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(SeUIDs,VisitIDs,VisitID,root_indirs,options)
  parms = mmil_args2parms(options,{...
    'SeUIDs',SeUIDs,[],...
    'VisitIDs',VisitIDs,[],...
    'VisitID',VisitID,[],...
    'root_indirs',root_indirs,[],...
    ...
    'ProjID',[],[],...
    'root_outdir','/space/syn05/1/data/MMILDB/DAL_ABCD',[],...
    'linkflag',false,[false true],...
    'follow_links_flag',false,[false true],...
    'forceflag',false,[false true],...
  });

  % determine cmd for copying/linking files
  if parms.linkflag
    parms.copy_cmd = 'ln -s';
  elseif parms.follow_links_flag
    parms.copy_cmd = 'cp -rpL';
  else
    parms.copy_cmd = 'cp -rp';
  end;

  % get number of series
  parms.nseries = length(parms.SeUIDs);

  % check root input directories
  parms.RootDirs = root_indirs;
%  parms.RootDirs.orig = [parms.root_indir '/orig'];
%  parms.RootDirs.raw = [parms.root_indir '/raw'];
%  parms.RootDirs = MMIL_Check_RootDirs(parms.RootDirs);

  % create root output directories
  parms.RootOutDirs = [];
  parms.RootOutDirs.orig = [parms.root_outdir '/orig'];
  parms.RootOutDirs.raw = [parms.root_outdir '/raw'];
  mmil_mkdir(parms.RootOutDirs.orig);
  mmil_mkdir(parms.RootOutDirs.raw);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = locate_indirs(parms)
  fprintf('%s: locating input directories...\n',mfilename);
  parms.orig_indirs = cell(parms.nseries,1);
  parms.raw_indirs = cell(parms.nseries,1);
  for i=1:parms.nseries
    VisitID = parms.VisitIDs{i};
    RootDirs = [];
    RootDirs.orig = [parms.RootDirs{i} '/orig'];
    RootDirs.raw = [parms.RootDirs{i} '/raw'];
    RootDirs = MMIL_Check_RootDirs(RootDirs);
    parms.orig_indirs{i} = MMIL_Get_Container(RootDirs,VisitID,'orig');
    if isempty(parms.orig_indirs{i})
      error('failed to locate input orig dir for %s',VisitID);
    end;
    parms.raw_indirs{i} = MMIL_Get_Container(RootDirs,VisitID,'raw');
    if isempty(parms.raw_indirs{i})
      error('failed to locate input raw dir for %s',VisitID);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_orig_outdir(parms)
  fprintf('%s: creating orig output directory...\n',mfilename);
  parms.orig_outdir = [parms.RootOutDirs.orig '/' parms.VisitID];
  mmil_mkdir(parms.orig_outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = compile_orig_dirs(parms);
  fname_contents = sprintf('%s/contents.mat',parms.orig_outdir);
  if ~exist(fname_contents,'file') || parms.forceflag
    if parms.linkflag
      fprintf('%s: linking orig series directories to %s...\n',...
        mfilename,parms.orig_outdir);
    else
      fprintf('%s: copying orig series directories to %s...\n',...
        mfilename,parms.orig_outdir);
    end;
    % loop over series
    orig_dirlist = cell(parms.nseries,1);
    for i=1:parms.nseries
      SeUID = parms.SeUIDs{i};
      VisitID = parms.VisitIDs{i};
      % load raw ContainerInfo for this series
      raw_indir = parms.raw_indirs{i};
      [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(raw_indir);
      if errcode
        error('failed to read ContainerInfo for %s',parms.VisitIDs{1});
      end;
      % determine series index for this series
      SeUIDs = {ContainerInfo.SeriesInfo.SeriesInstanceUID};
      k = find(strcmp(SeUID,SeUIDs));
      if isempty(k)
        error('series %s not found for visit %s',SeUID,VisitID);
      end;
      % determine source directory
      fname = ContainerInfo.SeriesInfo(k).OrigFileNames{1};
      fdir = fileparts(fname);
      [tmp,fdir,fext] = fileparts(fdir);
      fdir = [fdir fext];
      % copy from indir to outdir
      fname_in = sprintf('%s/%s',parms.orig_indirs{i},fdir);
      if ~exist(fname_in)
        error('input %s not found',fname_in);
      end;
      fname_out = sprintf('%s/%s',parms.orig_outdir,fdir);
      if ~exist(fname_out) || parms.forceflag
        if exist(fname_out) && parms.linkflag
          delete(fname_out);
        end;
        cmd = sprintf('%s %s %s',...
          parms.copy_cmd,fname_in,fname_out);
        fprintf('%s: %s\n',mfilename,cmd);
        [s,r] = unix(cmd);
        if s, error('cmd %s failed:\n%s',cmd,r); end;
      end;
      orig_dirlist{i} = fdir;
    end;
    VisitIDs = parms.VisitIDs;
    SeUIDs = parms.SeUIDs;
    save(fname_contents,'VisitIDs','SeUIDs','orig_dirlist');
  else
    load(fname_contents);
  end;
  parms.orig_dirlist = orig_dirlist;

  fname_contents = sprintf('%s/contents.txt',parms.orig_outdir);
  if ~exist(fname_contents,'file') || parms.forceflag
    % write file with list of input directories
    fid = fopen(fname_contents,'w');
    if fid<0, error('failed to open %s for writing',fname_contents); end;
    for i=1:parms.nseries
      fprintf(fid,'%s %s %s\n',...
        parms.VisitIDs{i},parms.SeUIDs{i},parms.orig_dirlist{i});
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_raw_outdir(parms)
  fprintf('%s: creating raw output directory...\n',mfilename);
  % sort parms.orig_dirlist
  [tmp,ind_sort] = sort(parms.orig_dirlist);
  % use first element to determine StudyDate and StudyTime
  k = ind_sort(1);
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(parms.raw_indirs{k});
  if errcode
    error('failed to read ContainerInfo in %s for %s',...
      parms.raw_indirs{k},parms.VisitIDs{k});
  end;
  ContainerUID = '1'; % as done in MMIL_Unpack_Dicoms
  ContainerOutDir = sprintf('MRIRAW_%s_%s.%s_%s',parms.VisitID,...
    ContainerInfo.StudyDate,ContainerInfo.StudyTime,ContainerUID);
  parms.raw_outdir = [parms.RootOutDirs.raw '/' ContainerOutDir];
  mmil_mkdir(parms.raw_outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compile_OrigSeriesInfo(parms)
  fname_out = sprintf('%s/OrigSeriesInfo.mat',parms.raw_outdir);
  if ~exist(fname_out,'file') || parms.forceflag
    fprintf('%s: compiling OrigSeriesInfo...\n',mfilename);
    SeriesInfo = [];
    for i=1:parms.nseries
      SeUID = parms.SeUIDs{i};
      VisitID = parms.VisitIDs{i};
      indir = parms.raw_indirs{i};
      fname_in = sprintf('%s/OrigSeriesInfo.mat',indir);
      tmp = load(fname_in);      
      tmp_info = tmp.SeriesInfo;
      % only include specific series from this visit
      tmp_info = select_series(tmp_info,SeUID,VisitID);
      % update file names
      tmp_info = update_OrigSeriesInfo_FileNames(tmp_info,i,parms);
      % concatenate info across series
      SeriesInfo = [SeriesInfo,tmp_info];
    end;
    save(fname_out,'SeriesInfo');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SeriesInfo = select_series(SeriesInfo,SeUID,VisitID)
  SeUIDs = {SeriesInfo.SeriesInstanceUID};
  k = find(strcmp(SeUID,SeUIDs));
  if isempty(k)
    error('series %s not found for visit %s',SeUID,VisitID);
  end;
  SeriesInfo = SeriesInfo(k);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SeriesInfo = update_OrigSeriesInfo_FileNames(SeriesInfo,v,parms)
  old_orig_dir = parms.orig_indirs{v};
  new_orig_dir = parms.orig_outdir;
  % replace FileNames
  for j=1:length(SeriesInfo)
    SeriesInfo(j).FileNames = regexprep(SeriesInfo(j).FileNames,old_orig_dir,new_orig_dir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compile_SeriesInfo(parms)
  fname_out = sprintf('%s/SeriesInfo.mat',parms.raw_outdir);
  if ~exist(fname_out,'file') || parms.forceflag
    fprintf('%s: compiling SeriesInfo...\n',mfilename);
    SeriesInfo = [];
    for i=1:parms.nseries
      SeUID = parms.SeUIDs{i};
      VisitID = parms.VisitIDs{i};
      indir = parms.raw_indirs{i};
      fname_in = sprintf('%s/SeriesInfo.mat',indir);
      tmp = load(fname_in);
      tmp_info = tmp.SeriesInfo;
      % only include specific series from this visit
      tmp_info = select_series(tmp_info,SeUID,VisitID);
      % update file names
      tmp_info = update_SeriesInfo_FileNames(tmp_info,i,parms);
      % concatenate info across series
      SeriesInfo = [SeriesInfo,tmp_info];
    end;
    save(fname_out,'SeriesInfo');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SeriesInfo = update_SeriesInfo_FileNames(SeriesInfo,v,parms)
  SUIDs = get_SUIDs(parms);
  old_orig_dir = parms.orig_indirs{v};
  new_orig_dir = parms.orig_outdir;
  % replace FileNames
  for j=1:length(SeriesInfo)
    SeriesNumber = SeriesInfo(j).SeriesNumber;
    StudyNumber = find(strcmp(SeriesInfo(j).StudyInstanceUID,SUIDs));
    old_raw_dir = SeriesInfo(j).SeriesDirPath;
    new_raw_dir = sprintf('%s/st%03d_ser%04d',...
      parms.raw_outdir,StudyNumber,SeriesNumber);
    SeriesInfo(j).OrigFileNames = ...
      regexprep(SeriesInfo(j).OrigFileNames,old_orig_dir,new_orig_dir);
    SeriesInfo(j).FileNames = ...
      regexprep(SeriesInfo(j).FileNames,old_raw_dir,new_raw_dir);
    SeriesInfo(j).SeriesDirPath = ...
      regexprep(SeriesInfo(j).SeriesDirPath,old_raw_dir,new_raw_dir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compile_ContainerInfo(parms)
  ContainerInfo = [];
  fname_out = sprintf('%s/ContainerInfo.mat',parms.raw_outdir);
  if ~exist(fname_out,'file') || parms.forceflag
    fprintf('%s: compiling ContainerInfo...\n',mfilename);
    SeriesInfo = [];
    for i=1:parms.nseries
      SeUID = parms.SeUIDs{i};
      VisitID = parms.VisitIDs{i};
      indir = parms.raw_indirs{i};
      fname_in = sprintf('%s/ContainerInfo.mat',indir);
      tmp = load(fname_in);
      if isempty(ContainerInfo)
        ContainerInfo = tmp.ContainerInfo;
      end;
      tmp_info = tmp.ContainerInfo.SeriesInfo;
      % only include specific series from this visit
      tmp_info = select_series(tmp_info,SeUID,VisitID);
      % update file names
      tmp_info = update_SeriesInfo_FileNames(tmp_info,i,parms);
      % concatenate info across series
      SeriesInfo = [SeriesInfo,tmp_info];
    end;    
    % sort series by StudyDate, StudyTime, SeriesDate, and SeriesTime
    SeriesInfo = mmil_sort_serinfo(SeriesInfo);
    % update summary info
    ContainerInfo.ContainerCreationDate = datestr(now);
    ContainerInfo.VisitID = parms.VisitID;
    ContainerInfo.SourceRootDir = parms.RootOutDirs.orig;
    ContainerInfo.SourceDir = parms.VisitID;
    if ~isempty(parms.ProjID)
      ContainerInfo.ProjID = parms.ProjID;
    end;
    ContainerInfo.StudyInstanceUIDlist = unique({SeriesInfo.StudyInstanceUID});
    ContainerInfo.SeriesInstanceUIDlist = unique({SeriesInfo.SeriesInstanceUID});
    ContainerInfo.SeriesInfo = SeriesInfo;
    ContainerInfo.ClassificationCompleted = datestr(now);
    save(fname_out,'ContainerInfo');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function summarize_SeriesInfo(parms)
  fname_out = sprintf('%s/SeriesInfo.csv',parms.raw_outdir);
  if ~exist(fname_out,'file') || parms.forceflag
    fprintf('%s: summarizing SeriesInfo in csv file...\n',mfilename);
    [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(parms.raw_outdir);
    if errcode
      error('failed to read ContainerInfo for %s',parms.VisitID);
    end;
    mmil_summarize_serinfo(ContainerInfo.SeriesInfo,fname_out,parms.forceflag);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compile_raw_dirs(parms)
  fname_contents = sprintf('%s/contents.mat',parms.raw_outdir);
  if ~exist(fname_contents,'file') || parms.forceflag
    fprintf('%s: copying raw series directories to %s...\n',...
      mfilename,parms.raw_outdir);
    SUIDs = get_SUIDs(parms);
    for i=1:parms.nseries
      SeUID = parms.SeUIDs{i};
      VisitID = parms.VisitIDs{i};
      indir = parms.raw_indirs{i};
      fname_in = sprintf('%s/ContainerInfo.mat',indir);
      tmp = load(fname_in);
      tmp_info = tmp.ContainerInfo.SeriesInfo;
      % only include specific series from this visit
      tmp_info = select_series(tmp_info,SeUID,VisitID);
      % copy subdirectory for this series
      SeriesNumber = tmp_info.SeriesNumber;
      StudyNumber = find(strcmp(tmp_info.StudyInstanceUID,SUIDs));
      old_dir = tmp_info.SeriesDirPath;
      new_dir = sprintf('%s/st%03d_ser%04d',...
        parms.raw_outdir,StudyNumber,SeriesNumber);
      if ~exist(new_dir) || parms.forceflag
        if exist(new_dir) && parms.linkflag
          delete(new_dir);
        end;
        cmd = sprintf('%s %s %s',parms.copy_cmd,old_dir,new_dir);
        fprintf('%s: %s\n',mfilename,cmd);
        [s,r] = unix(cmd);
        if s, error('cmd %s failed:\n%s',cmd,r); end;
      end;
    end;
    VisitIDs = parms.VisitIDs;
    SeUIDs = parms.SeUIDs;
    orig_dirlist = parms.orig_dirlist;
    save(fname_contents,'VisitIDs','SeUIDs','orig_dirlist');
  else
    load(fname_contents);
  end;

  fname_contents = sprintf('%s/contents.txt',parms.raw_outdir);
  if ~exist(fname_contents,'file') || parms.forceflag
    % write file with list of input directories
    fid = fopen(fname_contents,'w');
    if fid<0, error('failed to open %s for writing',fname_contents); end;
    for i=1:parms.nseries
      fprintf(fid,'%s %s %s\n',...
        parms.VisitIDs{i},parms.SeUIDs{i},parms.orig_dirlist{i});
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_progress_file(parms,step)
  fname = sprintf('%s/%s.txt',parms.raw_outdir,step);
  fid = fopen(fname,'a');
  if fid<0, error('failed to open %s for writing',fname); end;
  fprintf(fid,'%s completed: %s\n',step,datestr(now));
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SUIDs = get_SUIDs(parms)
  % load OrigSeriesInfo to get set of StudyInstanceUIDs
  SUIDs = [];
  SeriesInfo = [];
  fname_in = sprintf('%s/OrigSeriesInfo.mat',parms.raw_outdir);
  load(fname_in);
  SUIDs = unique({SeriesInfo.StudyInstanceUID});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

