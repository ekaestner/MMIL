function VisitInfo = MMIL_Select_DTI_Series(ProjID,varargin)
%function VisitInfo = MMIL_Select_Struct_Series(ProjID,[options])
%
% Purpose: visually compare and identify usable structural DTI series
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%
% Optional Input:
%  'fnamestem': stem of file name for processed data
%     {default = 'DTI'}
%  'type': series type as named in ContainerInfo.SeriesInfo
%     {default = 'DTI_ipp'}
%  'struct_type': structural series type as named in ContainerInfo.SeriesInfo
%     e.g. MPR, FLASHhi
%     need to have one of these scans for complete processing if
%       infix includes 'regT1'
%     {default = 'MPR'}
%  'infix': string attached to file name indicating processing steps
%     e.g. 'corr_regT1'
%     {default = 'corr_regT1'}
%  'ext': file extension
%     e.g. '.mgh','.mgz'
%     {default = '.mgz'}
%  'multi_flag': [0|1] only review data for visits with > 1 series
%     otherwise review all data
%     {default = 1}
%  'attempt_flag': [0|1] review data for multiple attempts together
%     VisitInfo must contain SubjID
%     {default = 1}
%  'attempt_window': number of days less than which two visits are considered
%     to be attempts for the same visit
%     {default = 90}
%  'add_missing_flag': [0|1] add to VisitInfo and SeriesSelection
%     visits in RootDirs.orig not found in existing csv files
%     {default = 1}
%  'SubjID_pattern': regexp pattern used for getting SubjID from VisitID
%     {default = '(?<SubjID>^[^_]+)_\w+'}
%  'rm_rejects_flag': [0|1] delete processing output for rejected series
%     from raw, proc, proc_dti, proc_bold, fsurf, fsico
%     for each subject with rejects in series_select that are not in visitinfo
%     afterwards, update VisitInfo file
%     {default = 1}
%  'imagine_flag': [0|1] use imagine to allow browsing of images
%     otherwise display static montages with multiple slices
%     {default = 0}
%  'b0_movie_flag': [0|1] display series of b=0 image slices
%     {default = 0}
%  'min_ndiffdirs': minimum number of diffusion directions
%     {default = 30}
%  'max_ndiffdirs': maximum number of diffusion directions
%     {default = 30}
%  'presave_flag': [0|1] load image data save save to mat file for later use
%     {default = 0}
%  'presave_dir': output directory for pre-saved image data
%     {default = pwd}
%
% Created:  09/25/13 by Matt Erhart
% Last Mod: 04/24/17 by Don Hagler
%

%% todo: change check_selection to only skip an attempt if all of the DTI
%%         series have either been accepted or rejected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all force
VisitInfo = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
fprintf('%s: checking input...\n',mfilename);
parms = check_input(ProjID,varargin);

% initialize SeriesSelection struct
fprintf('%s: initializing SeriesSelection...\n',mfilename);
SeriesSelection = init_SeriesSelection(parms);

% initialize VisitInfo struct
fprintf('%s: initializing VisitInfo...\n',mfilename);
VisitInfo = init_VisitInfo(parms);

% get SubjID from VisitID
fprintf('%s: getting SubjIDs from VisitIDs...\n',mfilename);
VisitInfo = get_SubjIDs(parms,VisitInfo);

% display instructions
display_instructions;

% loop over visits, reviewing data when more than one series
[SeriesSelection,VisitInfo] = review_data(parms,SeriesSelection,VisitInfo);

if parms.rm_rejects_flag
  % delete files depending on rejects and write VisitInfo file
  VisitInfo = remove_rejects(parms,SeriesSelection,VisitInfo);
  % write VisitInfo csv file
  if parms.really_rm_flag
    % sort rows by VisitID
    VisitInfo = sort_VisitInfo(parms,VisitInfo);
    save_VisitInfo(parms,VisitInfo);
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function display_instructions()
  fprintf('%s: at the OK? prompt, type:\n',mfilename);
  fprintf('  0:   restart subject\n');
  fprintf('  1:   comfirm and go to next\n');
  fprintf('  1q:  comfirm and quit\n');
  fprintf('  d:   discard choice and go to next\n');
  fprintf('  q:   discard choice and quit\n');
  fprintf('  1n:  comfirm, log notes, and go to next\n');
  fprintf('  1nq: comfirm, log notes, and quit\n');
  fprintf('  h:   display this message and restart subject\n');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
  ... % optional
    'fnamestem','DTI',{'DTI'},...
    'type','DTI_ipp',{'DTI_ipp','DTI','DTI_flex'},...
    'struct_type','MPR',{'MPR','FLASHhi'},...
    'infix','corr_regT1',[],...
    'ext','.mgz',{'.mgh','.mgz'},...
    'multi_flag',true,[false true],...
    'attempt_flag',true,[false true],...
    'attempt_window',90,[1,1e5],...
    'add_missing_flag',true,[false true],...
    'SubjID_pattern','(?<SubjID>^[^_]+)_\w+',[],...
    'rm_rejects_flag',true,[false true],...
    'imagine_flag',false,[false true],...
    'b0_movie_flag',false,[false true],...
    'min_ndiffdirs',30,[0 100],...
    'max_ndiffdirs',30,[0 100],...
    'presave_flag',false,[false true],...
    'presave_dir',pwd,[],...
  ... % hidden
    'really_rm_flag',true,[false true],...
    'outfix','DTI',[],...
    'img_dir','DTI_SeriesSelect',[],...
    'img_size',[6 4],[],...
    'img_dpi',200,[],...
    'img_titles',{'EXCLUDE','KEEP'},[],...
    'img_cropy',[0.1,0.9],[0,1],...
    'img_cropz',[0,1],[0,1],...
    'img_fract_slice',0.55,[],...
    'img_fract_slices',[0.40:.05:.70],[],...
    'img_title_fractx',0.05,[0,1],...
    'img_title_fracty',0.05,[0,1],...
    'img_b0_scale_range',[5 99],[],...
    'img_low_scale_range',[4 99.9],[],...
    'img_high_scale_range',[4 99.9],[],...
    'RootDirs',[],[],...
    'StudyInfo',[],[],...
    'VisitIDs',[],[],...
    'SubjIDs',[],[],...
    'user',[],[],...
    'rmdirlist',{'raw','proc_dti'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  % set ignore_csv_flag
  parms.ignore_VisitInfo_flag = parms.add_missing_flag;
  % use ProjID to get StudyInfo and RootDirs
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  try
    [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  catch me
    fprintf('%s: WARNING: %s\n',mfilename,me.message);
    ProjInfo = [];
  end;
  % check StudyInfo, get VisitIDs
  if length(parms.StudyInfo)==0
    error('no valid subjects');
  end;
  parms.VisitIDs = {parms.StudyInfo.VisitID};
  % construct file suffix from infix and ext
  if isempty(parms.infix)
    parms.suffix = parms.infix;
  else
    parms.suffix = ['_' parms.infix];
  end;
  parms.suffix = [parms.suffix parms.ext];
  % set series selection file name
  parms.fname_SeriesSelection = ...
    sprintf('%s/ProjInfo/%s/%s_%s_Selection.csv',...
      parms.RootDirs.home,ProjID,ProjID,parms.outfix);
  % set VisitInfo file name
  parms.fname_VisitInfo = sprintf('%s/ProjInfo/%s/%s_VisitInfo.csv',...
    parms.RootDirs.home,ProjID,ProjID);
  % set output dir for images
  if mmil_isrelative(parms.img_dir)
    parms.img_dir = sprintf('%s/MetaData/%s/%s',...
      parms.RootDirs.home,parms.ProjID,parms.img_dir);
  end;
  mmil_mkdir(parms.img_dir);
  if parms.imagine_flag
    close all force
  end;
  if parms.min_ndiffdirs > parms.max_ndiffdirs
    error('max_ndiffdirs is less than min_ndiffdirs');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SeriesSelection = init_SeriesSelection(parms)
  % create SeriesSelection csv file if does not exist
  if ~exist(parms.fname_SeriesSelection,'file')
    headers = {'AcceptSeries','RejectSeries','notes'};
    data = cell(length(parms.VisitIDs),length(headers));
    mmil_write_csv(parms.fname_SeriesSelection,data,'firstcol_label','VisitID',...
      'row_labels', parms.VisitIDs, 'col_labels', headers);
  end;
  SeriesSelection = mmil_csv2struct(parms.fname_SeriesSelection);
  % add new subjects to SeriesSelection
  if parms.add_missing_flag
    ind_new = ... index of new subjects to add
      parms.VisitIDs(~ismember(parms.VisitIDs,{SeriesSelection.VisitID}));
    for v = 1:length(ind_new) 
      SeriesSelection(end+1).VisitID = ind_new{v};
    end
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitInfo = init_VisitInfo(parms)
  % create VisitInfo csv file if does not exist
  if ~exist(parms.fname_VisitInfo,'file')
    headers = {'DCM_RejectSeries'};
    data = cell(length(parms.VisitIDs),length(headers));
    mmil_write_csv(parms.fname_VisitInfo,data,'firstcol_label','VisitID',...
      'row_labels', parms.VisitIDs, 'col_labels', headers);
  end;
  VisitInfo = mmil_csv2struct(parms.fname_VisitInfo);
  % add new subjects to VisitInfo
  if parms.add_missing_flag
    ind_new = ... index of new subjects to add
      parms.VisitIDs(~ismember(parms.VisitIDs,{VisitInfo.VisitID}));
    for v = 1:length(ind_new) 
      VisitInfo(end+1).VisitID = ind_new{v};
    end
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitInfo = get_SubjIDs(parms,VisitInfo)
  nvisits = length(VisitInfo);
  for v=1:nvisits
    if ~isfield(VisitInfo,'SubjID'), VisitInfo(v).SubjID = []; end;
    if isempty(VisitInfo(v).SubjID)
      % get SubjID from VisitID using regexp
      n = regexp(VisitInfo(v).VisitID,'(?<SubjID>^[^_]+)_\w+','names');
      if ~isempty(n)
        VisitInfo(v).SubjID = n.SubjID;
      else
        fprintf('%s: WARNING: VisitID %s does not match SubjID_pattern\n',...
          mfilename,VisitInfo(v).VisitID);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SeriesSelection,VisitInfo] = ...
                                review_data(parms,SeriesSelection,VisitInfo)
  quit_flag = false;
  nvisits = length(SeriesSelection);

  % disable attempt_flag if all SubjIDs are empty
  if all(cellfun(@isempty,{VisitInfo.SubjID}))
    parms.attempt_flag = false;
  end;

  %% todo: find a way to identify the groups of attempts first
  %%    then loop over groups of attempts
  %%    to avoid duplicate messages

  for ind_select=1:nvisits
    clear vol qmat censor volb0;
    tic
    VisitID = SeriesSelection(ind_select).VisitID;

    % find matching visit in StudyInfo or skip
    ind_study = find(strcmp(VisitID,parms.VisitIDs));
    if isempty(ind_study), continue; end;
    
    % determine number of attempts for this visit
    [nattempts,ind_study_attempts,ind_select_attempts] = ...
      get_nattempts(parms,SeriesSelection,VisitInfo,ind_study,ind_select);

    % if multiple attempts, set outstem = concatenated VisitIDs
    outstem = [];
    for j=1:nattempts
      if ~isempty(outstem), outstem = [outstem '_+_']; end;
      outstem = [outstem parms.StudyInfo(ind_study_attempts(j)).VisitID];
    end;

    % skip if selection already made
    skip_flag = check_selection(SeriesSelection,ind_select_attempts);
    if skip_flag, continue; end;

    % find scans for each attempt
    [ind_series,ind_scan,...
     ind_study_series,ind_select_series,...
     auto_rejects,zero_reason] =...
       check_scans(parms,SeriesSelection,VisitInfo,ind_study_attempts,ind_select_attempts);
    nscans = length(ind_select_series);

    % accept if only one scan per visit if multiflag == 1
    result = [];
    if nscans == 1 && parms.multi_flag
      SeriesSelection(ind_select_series).AcceptSeries = ind_series;
    elseif nscans>0
      % load data
      [vol,volb0,qmat,censor,all_exist_flag] = ...
        load_data(parms,ind_study_series,ind_scan,outstem);
      if ~all_exist_flag, continue; end;

      if parms.presave_flag
        fname_presave = sprintf('%s/%s.mat',parms.presave_dir,outstem);
        fprintf('%s: saving %s...\n',mfilename,fname_presave);
        mmil_mkdir(parms.presave_dir);
        save(fname_presave, 'vol','volb0','qmat','censor');
        clear vol volb0 qmat censor
        continue;
      end

      fprintf('%s: creating images...\n',mfilename);

      % show b=0 movie to check for jitter/motion between slices
      if parms.b0_movie_flag
        show_b0_movies(volb0);
      end;

      % show motion cesoring images
      show_censor_images(censor);

      % show motion stats
      show_motion_stats(qmat,censor);

      if parms.imagine_flag
        % show data images
        imagine(vol{:});
      else
        % show multi-slce b=0 images
        show_b0_images(volb0,parms);
        % show multi-frame diffusion-weighted images
        show_dwi_images(vol,qmat,parms);
      end;

      % prompt user for selections
      restart_flag = true;
      discard_flag = false;
      while restart_flag | isempty(result)
        restart_flag = false;
        for j=1:nscans
          try
            % prompt user
            result(j,1) = input(sprintf('use %s %d? [1 | 0] ',parms.type,j));
          catch me
            restart_flag = true;
            continue;
          end
        end
        if length(result) == nscans
          % show and save image of selected scan
          save_selected_image(vol,qmat,result,outstem,parms)
          % confirm QC, quit, add notes
          [quit_flag,restart_flag,discard_flag,notes] = ...
            check_result(result,quit_flag,restart_flag,discard_flag);
          for j=1:nattempts
            ind_select = ind_select_attempts(j);
            SeriesSelection(ind_select).notes = notes;
          end;
          if discard_flag
            if quit_flag, return; end;
            continue;
          end;
        else
          restart_flag = true;
          continue;
        end
      end
      toc
      if parms.imagine_flag
        close all force
      end;
      clear vol qmat censor volb0;
    end;

    % update SeriesSelection
    SeriesSelection = update_selection(SeriesSelection,ind_series,result,...
      ind_select_attempts,ind_select_series,auto_rejects,zero_reason);

    % update SeriesSelection csv file
    mmil_write_csv(parms.fname_SeriesSelection,...
      squeeze(struct2cell(SeriesSelection))',...
      'col_labels',fieldnames(SeriesSelection)')

    % remove pre-saved data
    fname_presave = sprintf('%s/%s.mat',parms.presave_dir,outstem);
    if exist(fname_presave,'file')
      fprintf('%s: deleting %s...\n',mfilename,fname_presave);
      unix(['rm -v ' fname_presave]);
      fname_presave = [];
    end
     
    fprintf('%s: wrote selection for %s\n',mfilename,outstem)
    if quit_flag; return; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nattempts,ind_study_attempts,ind_select_attempts] = ...
           get_nattempts(parms,SeriesSelection,VisitInfo,ind_study,ind_select)
  SS_VisitIDs = {SeriesSelection.VisitID};
  VI_VisitIDs = {VisitInfo.VisitID};
  VI_SubjIDs = {VisitInfo.SubjID};
  VisitID = SeriesSelection(ind_select).VisitID;
  StudyDate = parms.StudyInfo(ind_study).StudyDate;
  ind_study_attempts = ind_study;
  ind_select_attempts = ind_select;
  if parms.attempt_flag && ~isempty(StudyDate)
    ind_visit = find(strcmp(VisitID,VI_VisitIDs));
    SubjID = VI_SubjIDs{ind_visit};
    ind_attempts = find(strcmp(SubjID,VI_SubjIDs));
    if length(ind_attempts) > 1
      % set ind_study_attempts and ind_select_attempts
      for j=1:length(ind_attempts)
        tmp_VisitID = VI_VisitIDs{ind_attempts(j)};
        if strcmp(tmp_VisitID,VisitID), continue; end;
        ind_study_tmp = find(strcmp(tmp_VisitID,parms.VisitIDs));
        if isempty(ind_study_tmp), continue; end;
        ind_select_tmp = find(strcmp(tmp_VisitID,SS_VisitIDs));
        if isempty(ind_select_tmp), continue; end;
        % include attempts within the attempt_window
        tmpStudyDate = parms.StudyInfo(ind_study_tmp).StudyDate;
        if ~isempty(tmpStudyDate)
          ndays = abs(daysact(datenum(num2str(StudyDate),'yyyymmdd'),...
                              datenum(num2str(tmpStudyDate),'yyyymmdd')));
          if ndays <= parms.attempt_window
            ind_study_attempts = cat(2,ind_study_attempts,ind_study_tmp);
            ind_select_attempts = cat(2,ind_select_attempts,ind_select_tmp);
          end;
        end;
      end;
    end;
  end;
  nattempts = length(ind_study_attempts);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function skip_flag = check_selection(SeriesSelection,ind_select_attempts)
  skip_flag = 0;
  nattempts = length(ind_select_attempts);
  selection_done = zeros(nattempts,1);
  for j=1:nattempts
    ind_select = ind_select_attempts(j);
    if ~isempty(SeriesSelection(ind_select).AcceptSeries) ||...
       ~isempty(SeriesSelection(ind_select).RejectSeries) 
       selection_done(j) = 1;
    end
  end;
  % skip if a series has been accepted or rejected for each attempt
  if all(selection_done)
    skip_flag = 1;
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_series,ind_scan,...
          ind_study_series,ind_select_series,...
          auto_rejects,zero_reason] =...
      check_scans(parms,SeriesSelection,VisitInfo,ind_study_attempts,ind_select_attempts);

  VI_VisitIDs = {VisitInfo.VisitID};
  nattempts = length(ind_select_attempts);
  ind_series = [];
  ind_scan = [];
  ind_study_series = [];
  ind_select_series = [];
  auto_rejects = cell(nattempts,1);
  zero_reason = cell(nattempts,1);
  for j=1:nattempts
    ind_study = ind_study_attempts(j);
    ind_select = ind_select_attempts(j);
    VisitID = SeriesSelection(ind_select).VisitID;
    % check for more than one scan
    if isempty(parms.StudyInfo(ind_study).raw), 
      fprintf('%s: no raw container for %s\n',mfilename,VisitID);
      continue;
    end;
    % check raw for dti
    RawPath = sprintf('%s/%s',...
      parms.RootDirs.raw,parms.StudyInfo(ind_study).raw);
    [RawInfo, errcode] = MMIL_Load_ContainerInfo(RawPath);
    if errcode
      fprintf('%s: WARNING: error loading raw ContainerInfo for %s in %s\n',...
        mfilename,VisitID,RawPath);
      continue;
    end;
    if ~isfield(RawInfo.SeriesInfo,'SeriesType'); 
      fprintf('%s: raw ContainerInfo.SeriesInfo.SeriesType not defined for %s\n',...
        mfilename,VisitID);
      continue;
    end 
    % auto reject visits without any T1 scans according to rawinfo
    %   if infix includes 'regT1'
    rawT1 = strcmp({RawInfo.SeriesInfo.SeriesType},parms.struct_type);
    ind_visit = find(strcmp(VisitID,VI_VisitIDs));
    STRUCT_VisitID = VisitInfo(ind_visit).STRUCT_VisitID;
    if ~isempty(regexp(parms.infix,'regT1')) && ~any(rawT1) &&...
        isempty(STRUCT_VisitID)
      fprintf('%s: WARNING: missing T1 scan for %s - use STRUCT_VisitID if possible\n',...
        mfilename,VisitID);
      if ~any(rawT1)
        zero_reason{j} = sprintf('no %s',parms.struct_type);
        continue;
      end
    end;
    % auto reject visits without any DTI scans according to rawinfo
    rawDTI = strcmp({RawInfo.SeriesInfo.SeriesType},parms.type);
    if ~any(rawDTI)
      zero_reason{j} = sprintf('no %s',parms.type);
      continue;
    end
    ContainerPath = sprintf('%s/%s',...
      parms.RootDirs.proc_dti,parms.StudyInfo(ind_study).proc_dti);
    [ContainerInfo, errcode] = MMIL_Load_ContainerInfo(ContainerPath);
    if errcode
      fprintf('%s: WARNING: error loading proc_dti ContainerInfo for %s\n',...
        mfilename,VisitID);
      continue;
    end;
    nscans = ContainerInfo.DTI_cntr;
    % redundant check
    if nscans == 0;
      fprintf('%s: no %s scans for %s\n',mfilename,parms.type,VisitID);
      continue;
    end;
    ind_series_tmp = [ContainerInfo.ScanInfo.DTI.SeriesIndex]';
    ind_scan_tmp = [1:nscans]';
    isDTI = strcmp({ContainerInfo.ScanInfo.DTI.ScanType},parms.type);

    % exclude invalid DTI series
    isValid = [ContainerInfo.ScanInfo.DTI.valid];
    if any(isDTI & ~isValid)
      fprintf('%s: WARNING: rejecting invalid scans for %s\n',...
        mfilename,VisitID);
      auto_rejects{j} = ind_series_tmp(isDTI & ~isValid);
    end;
    isDTI = isDTI & isValid;
    if ~any(isDTI)
      zero_reason{j} = sprintf('no valid %s',parms.type);
      continue;
    end;
    ind_series_tmp = ind_series_tmp(isDTI);
    ind_scan_tmp = ind_scan_tmp(isDTI);

    % exclude short scans (e.g. aborted)
    ndtiframes = [ContainerInfo.ScanInfo.DTI(isDTI).ndiffdirs];
    nframesok = (ndtiframes >= parms.min_ndiffdirs &...
                 ndtiframes <= parms.max_ndiffdirs);
    auto_rejects{j} = union(auto_rejects{j},ind_series_tmp(~nframesok));                 
    if all(~nframesok)
      zero_reason{j} = sprintf('no full length %s',parms.type);
      continue;
    end;
    if any(~nframesok)
      fprintf('%s: WARNING: rejecting short scan(s) for %s (number of frames = [%s])\n',...
        mfilename,VisitID,strtrim(sprintf('%d ',ndtiframes(~nframesok))));
    end;
    ind_series_tmp = ind_series_tmp(nframesok);
    ind_scan_tmp = ind_scan_tmp(nframesok);
    
    nscans = length(find(nframesok)); % only dti scans
    % exclude scans that have already been rejected
    ind_reject = SeriesSelection(ind_select).RejectSeries;
    if ~isempty(ind_reject)
      [ind_series_tmp,ind_keep] = setdiff(ind_series_tmp,ind_reject);
      ind_scan_tmp = ind_scan_tmp(ind_keep);
      nscans = length(ind_scan_tmp);
    end;
    % add scans for selection
    if ~isempty(ind_series_tmp)
      ind_series = cat(1,ind_series,ind_series_tmp);
      ind_scan = cat(1,ind_scan,ind_scan_tmp);
      ind_study_series = ...
        cat(1,ind_study_series,ind_study*ones(nscans,1)); 
      ind_select_series = ...
        cat(1,ind_select_series,ind_select*ones(nscans,1));
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol,volb0,qmat,censor,all_exist_flag] = ...
    load_data(parms,ind_study_series,ind_scan,outstem)
  vol = []; volb0 = []; qmat = []; censor = []; all_exist_flag = 1;
  nscans = length(ind_study_series);
  fname_presave = sprintf('%s/%s.mat',parms.presave_dir,outstem);
  if exist(fname_presave,'file')
    fprintf('%s: loading pre-saved mat file for %s...\n',mfilename,outstem);
    load(fname_presave);
  else
    fprintf('%s: checking data files for %s...\n',mfilename,outstem);
    for j=1:nscans
      ind_study_tmp = ind_study_series(j);
      ind_scan_tmp = ind_scan(j);
      tmp_VisitID = parms.StudyInfo(ind_study_tmp).VisitID;
      ContainerPath = sprintf('%s/%s',...
        parms.RootDirs.proc_dti,parms.StudyInfo(ind_study_tmp).proc_dti);
      [fname_data,fname_b0for,fname_b0rev,fname_censor,fname_qmat] =...
          set_data_fnames(ContainerPath,parms,ind_scan_tmp);
      fname_list = {fname_data,fname_b0for,fname_b0rev,fname_censor,fname_qmat};
      for f=1:length(fname_list)
        tmp_fname = fname_list{f};
        if ~exist(tmp_fname)
          fprintf('%s: WARNING: for VisitID %s, %s not found\n',...
            mfilename,tmp_VisitID,tmp_fname);
          all_exist_flag = 0;
          break;
        end;
      end;
      if ~all_exist_flag, return; end;
    end;

    fprintf('%s: loading data for %s...\n',mfilename,outstem);
    vol = {}; volb0 = {};
    for j=1:nscans 
      ind_study_tmp = ind_study_series(j);
      ind_scan_tmp = ind_scan(j);
      tmp_VisitID = parms.StudyInfo(ind_study_tmp).VisitID;
      ContainerPath = sprintf('%s/%s',...
        parms.RootDirs.proc_dti,parms.StudyInfo(ind_study_tmp).proc_dti);
      [fname_data,fname_b0for,fname_b0rev,fname_censor,fname_qmat] =...
          set_data_fnames(ContainerPath,parms,ind_scan_tmp);
      try
        %load motion information
        if isempty(qmat)
          clear qmat censor
        end;
        qmat(j) = load(fname_qmat);
        censor(j) = load(fname_censor);
        % load raw b=0
        [volb0_forward,M] = fs_load_mgh(fname_b0for,[],1);
        volb0_forward = fs_reorient(volb0_forward,M,'PRI');
        [volb0_rev,M] = fs_load_mgh(fname_b0rev,[],1);
        volb0_rev = fs_reorient(volb0_rev,M,'PRI');
        %% todo: keep forward and rev separate? (for cropping z)
        volb0{j} = cat(3,volb0_forward,volb0_rev);
        % load main DTI
        [vol_tmp,M] = fs_load_mgh(fname_data);
        vol{j} = fs_reorient(vol_tmp,M,'PRI');
        clear vol_tmp volb0_forward volb0_rev M;
      catch me
        fprintf('%s: WARNING: for VisitID %s, failed to load %s:\n%s\n',...
          mfilename,tmp_VisitID,fname_data,me.message);
        all_exist_flag = 0;
        break;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_data,fname_b0for,fname_b0rev,fname_censor,fname_qmat] =...
    set_data_fnames(ContainerPath,parms,ind_scan_tmp)
  %% todo: load raw data (to show dark slices)?
  % processed data
  fname_data = sprintf('%s/%s%d%s',...
    ContainerPath,parms.fnamestem,ind_scan_tmp,parms.suffix);
  %% todo: load processed b=0 images (to show B0 correction)?
  % raw b=0 images
  fname_b0for = sprintf('%s/%s%d%s',...
    ContainerPath,parms.fnamestem,ind_scan_tmp,parms.ext);
  fname_b0rev = sprintf('%s/%s%d_rev%s',...
    ContainerPath,parms.fnamestem,ind_scan_tmp,parms.ext);
  % censor file
  fname_censor = sprintf('%s/%s%d_censor.mat',...
    ContainerPath,parms.fnamestem,ind_scan_tmp);
  % qmat file with motion estimates
  if isempty(parms.infix)
    tmp_suffix = parms.infix;
  else
    tmp_suffix = ['_' parms.infix];
  end;
  tmp_suffix = [tmp_suffix '_qmat.mat'];
  fname_qmat = sprintf('%s/%s%d%s',...
    ContainerPath,parms.fnamestem,ind_scan_tmp,tmp_suffix);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_b0_movies(vols,fnum)
  if ~exist('fnum','var') || isempty(fnum), fnum=100; end;
  nscans = length(vols);
  vals = [];
  for j=1:nscans
    tmp_vals = mmil_rowvec(vols{j});
    vals = cat(2,vals,tmp_vals(1:10:end));
  end;
  clim = prctile(vals,parms.img_b0_scale_range);
  nslices = unique(cell2mat(cellfun(@(x) size(x,3), vols,'uni',0)));
  if length(nslices) > 1; 
    fprintf('%s: WARNING: b=0 images have %s slices\n',...
      fprintf('%d ',num2str(cell2mat(cellfun(@(x) size(x,3), vols,'uni',0)))));
  end;
  f = figure(fnum); clf;
  set(f,'position',[500 500 500*nscans 600]);
  for z = 1:nslices
    for ns = 1:nscans
      figure(fnum);
      % check if we have n > nslices for this scan, if so, display zeros
      if z<=size(vols{ns},3)
        im = vols{ns}(:,:,z);
      else
        im = zeros(size(vols{ns},1),size(vols{ns},2));
      end;
      subplot_tight(1,length(vols),ns);
      imagesc(im);
      caxis(clim); colormap(bone);
    end;
    if nscans == 1
      pause(.1);
    end
  end;
  close(f);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_censor_images(censor,fnum)
  if ~exist('fnum','var') || isempty(fnum), fnum=1; end;
  nscans = length(censor);
  figure(fnum); clf;
  set(gcf,'name','censor err');
  for c = 1:nscans
    subplot(nscans,2,(c*2)-1); imagesc(censor(c).censor_mat); 
    title(['Scan ' num2str(c)]); caxis([0 3]);
    subplot(nscans,2,c*2); imagesc(censor(c).censor_err);
    title(['Scan ' num2str(c)]); caxis([0 3]);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_motion_stats(qmat,censor,fnum)
  if ~exist('fnum','var') || isempty(fnum), fnum=2; end;
  nscans = length(qmat);
  figure(fnum); clf;
  set(gcf,'name', 'Motion Stats');
  motion = [];
  for ns = 1:nscans
    motion.frame_slices(ns) = sum(sum(censor(ns).censor_mat));
    motion.frames(ns) = sum(sum(censor(ns).censor_mat) > 0);
    motion.mean_motion(ns) = qmat(ns).mean_motion;
    nbadframes = sum(censor(ns).censor_mat,2);
    motion.max_bad_frames(ns) = max(nbadframes);
  end
  motion_fields = fieldnames(motion);
  for mp = 1:length(motion_fields)
    figure(2);
    subplot(1,length(motion_fields),mp); 
    scatter(repmat(1,[1 nscans]), motion.(motion_fields{mp}),'filled');
    text(repmat(1,[1 nscans])+.1, motion.(motion_fields{mp}),...
      cellfun(@num2str, num2cell([1:nscans]),'uni',0), 'FontSize',20)
    title(regexprep(motion_fields{mp},'_',' '));
    %% todo: set ylim to reasonable values for each?
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_b0_images(vols,parms,name,fnum)
  if ~exist('name','var') || isempty(name), name='b=0 images (for vs rev)'; end;
  if ~exist('fnum','var') || isempty(fnum), fnum=3; end;
  nscans = length(vols);
  nslices = length(parms.img_fract_slices);
  % set scaling range
  vals = [];
  for j=1:nscans
    tmp_vals = mmil_rowvec(vols{j});
    vals = cat(2,vals,tmp_vals(1:10:end));
  end;
  clim = prctile(vals,parms.img_b0_scale_range);
  % plot images
  figure(fnum); clf;
  set(gcf,'name',name);
  i = 1;
  for j=1:nscans
    [nx,ny,nz] = size(vols{j});
    for k=1:nslices
      img_slice = round(parms.img_fract_slices(k)*ny);
      x0 = min(max(1,round(parms.img_cropy(1)*nx)),nx);
      x1 = min(max(1,round(parms.img_cropy(2)*nx)),nx);
      z0 = min(max(1,round(parms.img_cropz(1)*nz)),nz);
      z1 = min(max(1,round(parms.img_cropz(2)*nz)),nz);
      subplot_tight(nscans,nslices,i);
      imagesc(flipud(rot90(squeeze(vols{j}(x0:x1,img_slice,z0:z1)))));
      axis off;
      colormap(gray); caxis(clim);
      tstr = sprintf('%s%d s%d',parms.fnamestem,j,img_slice);
      x = parms.img_title_fractx*nx;
      y = parms.img_title_fracty*nz;
      text(x,y,tstr,'color','w','fontsize',12);
      i = i + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_dwi_images(vols,qmat,parms,name,fnum)
  if ~exist('name','var') || isempty(name), name='diffusion weighted images'; end;
  if ~exist('fnum','var') || isempty(fnum), fnum=4; end;
  nscans = length(vols);
  nslices = length(parms.img_fract_slices);
  % set scaling ranges
  vals_b0 = [];
  vals_low = [];
  vals_high = [];
  max_nf = 0;
  for j=1:nscans
    bvals = qmat(j).bvals;
    nf = length(bvals);
    if nf > max_nf
      max_nf = nf;
    end;
    ind_b0 = find(bvals<=100);
    ind_low = find(bvals>100 & bvals<=1500);
    ind_high = find(bvals>1500);
    tmp_vals = mmil_rowvec(vols{j}(:,:,:,ind_b0));
    vals_b0 = cat(2,vals_b0,tmp_vals(1:10:end));
    tmp_vals = mmil_rowvec(vols{j}(:,:,:,ind_low));
    vals_low = cat(2,vals_low,tmp_vals(1:50:end));
    tmp_vals = mmil_rowvec(vols{j}(:,:,:,ind_high));
    vals_high = cat(2,vals_high,tmp_vals(1:50:end));
  end;
  clim_b0 = prctile(vals_b0,parms.img_b0_scale_range);
  clim_low = prctile(vals_low,parms.img_low_scale_range);
  clim_high = prctile(vals_high,parms.img_high_scale_range);
  
  % determine subplot dimensions
  nc = floor(sqrt(max_nf*nscans));
  nr = nscans*ceil(max_nf/nc);
  % show images of each frame
  figure(fnum); clf;
  set(gcf,'name',name);
  for j=1:nscans
    i = 1 + nc*(nr/nscans)*(j-1);
    [nx,ny,nz,nf] = size(vols{j});
    bvals = qmat(j).bvals;
    ind_b0 = find(bvals<=100);
    ind_low = find(bvals>100 & bvals<=1500);
    ind_high = find(bvals>1500);
    img_slice = round(parms.img_fract_slice*ny);
    x0 = min(max(1,round(parms.img_cropy(1)*nx)),nx);
    x1 = min(max(1,round(parms.img_cropy(2)*nx)),nx);
    z0 = min(max(1,round(parms.img_cropz(1)*nz)),nz);
    z1 = min(max(1,round(parms.img_cropz(2)*nz)),nz);
    % loop over frames and show a single slice
    for k=1:nf
      subplot_tight(nr,nc,i);
      imagesc(flipud(rot90(squeeze(vols{j}(x0:x1,img_slice,z0:z1,k)))));
      axis off;
      colormap(gray);
      if ismember(k,ind_b0)
        caxis(clim_b0);
      elseif ismember(k,ind_low)
        caxis(clim_low);
      else
        caxis(clim_high);
      end;
      tstr = sprintf('%s %d frame %d',...
        parms.fnamestem,j,k);
      x = parms.img_title_fractx*nx;
      y = parms.img_title_fracty*nz;
      text(x,y,tstr,'color','w','fontsize',12);
      i = i + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_selected_image(vols,qmat,result,outstem,parms,fnum)
  if ~exist('fnum','var') || isempty(fnum), fnum=4; end;
  nscans = length(vols);
  scan2plot = find(result);
  if isempty(scan2plot), return; end;
  figure(fnum); clf;
  for f = 1:length(scan2plot)
    j = scan2plot(f);
    [nx,ny,nz,nf] = size(vols{j});
    % determine scaling range
    bvals = qmat(j).bvals;
    ind_b0 = find(bvals<=100);
    ind_low = find(bvals>100 & bvals<=1500);
    ind_high = find(bvals>1500);
    tmp_vals = mmil_rowvec(vols{j}(:,:,:,ind_b0));
    vals_b0 = tmp_vals(1:10:end);
    tmp_vals = mmil_rowvec(vols{j}(:,:,:,ind_low));
    vals_low = tmp_vals(1:50:end);
    tmp_vals = mmil_rowvec(vols{j}(:,:,:,ind_high));
    vals_high = tmp_vals(1:50:end);
    clim_b0 = prctile(vals_b0,parms.img_b0_scale_range);
    clim_low = prctile(vals_low,parms.img_low_scale_range);
    clim_high = prctile(vals_high,parms.img_high_scale_range);
    % determine subplot dimensions
    nc = floor(sqrt(nf));
    nr = ceil(nf/nc);
    % plot images
    set(gcf,'name',sprintf('%s scan %d ACCEPTED',parms.fnamestem,j));
    img_slice = round(parms.img_fract_slice*ny);
    x0 = min(max(1,round(parms.img_cropy(1)*nx)),nx);
    x1 = min(max(1,round(parms.img_cropy(2)*nx)),nx);
    z0 = min(max(1,round(parms.img_cropz(1)*nz)),nz);
    z1 = min(max(1,round(parms.img_cropz(2)*nz)),nz);
    % loop over frames and show a single slice
    for k=1:nf
      subplot_tight(nr,nc,k);
      imagesc(flipud(rot90(squeeze(vols{j}(x0:x1,img_slice,z0:z1,k)))));
      axis off;
      colormap(gray);
      if ismember(k,ind_b0)
        caxis(clim_b0);
      elseif ismember(k,ind_low)
        caxis(clim_low);
      else
        caxis(clim_high);
      end;
      tstr = sprintf('%s%d f%d',...
        parms.fnamestem,j,k);
      x = parms.img_title_fractx*nx;
      y = parms.img_title_fracty*nz;
      text(x,y,tstr,'color','w','fontsize',12);
    end;
  end
  % save image
  set(gcf, 'PaperUnits', 'inches');
  set(gcf, 'PaperSize', [parms.img_size(1) parms.img_size(2)]);
  set(gcf, 'PaperPositionMode', 'manual');
  set(gcf, 'PaperPosition', [0 0 parms.img_size(1) parms.img_size(2)]);
  print(gcf,'-dtiff',...
    [parms.img_dir '/' outstem '.tif'],sprintf('-r %d',parms.img_dpi));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SeriesSelection = update_selection(SeriesSelection,ind_series,result,...
                ind_select_attempts,ind_select_series,auto_rejects,zero_reason);

  nattempts = length(ind_select_attempts);
  for j=1:nattempts
    ind_select = ind_select_attempts(j);
    if isempty(result)
      tmp_accept = SeriesSelection(ind_select).AcceptSeries;
    else
      tmp_accept = ind_series(result==1 & ind_select_series==ind_select);
    end;
    if ~isempty(tmp_accept)
      SeriesSelection(ind_select).AcceptSeries = tmp_accept;
    else
      SeriesSelection(ind_select).AcceptSeries = [];
    end;
    tmp_reject = union(...
      SeriesSelection(ind_select).RejectSeries,...
      ind_series(result==0 & ind_select_series==ind_select));
    % if scan short add to rejects
    if ~isempty(auto_rejects{j})
      tmp_reject = union(tmp_reject,auto_rejects{j});
    end
    if ~isempty(tmp_reject)
      SeriesSelection(ind_select).RejectSeries = ... append rejects
        union(tmp_reject, SeriesSelection(ind_select).RejectSeries);
    else
      SeriesSelection(ind_select).RejectSeries = [];
    end;
    if ~isempty(zero_reason{j}) &&...
       isempty(SeriesSelection(ind_select).RejectSeries)
      SeriesSelection(ind_select).RejectSeries = 0;
      SeriesSelection(ind_select).notes = zero_reason{j};
      fprintf('%s: %s has %s scans: marking as 0 and skipping\n',...
        mfilename,SeriesSelection(ind_select).VisitID,zero_reason{j});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [quit_flag,restart_flag,discard_flag,notes] = ...
                       check_result(result,quit_flag,restart_flag,discard_flag)
  notes = [];
  check_flag = ...
    input([num2str(result') ' OK? [0 | 1 | 1q | d | q | 1n | 1nq | h] '],'s');
  switch check_flag 
    case '0'
      fprintf('%s: restarting subject...\n',mfilename);
      restart_flag = true;
    case '1'
      fprintf('%s: selection confirmed, going to next subject...\n',...
        mfilename);
    case '1q'
      fprintf('%s: selection confirmed, quitting now\n',...
        mfilename);
      quit_flag = true;
    case 'd'
      fprintf('%s: discarding selection, going to next subject...\n',...
        mfilename);
      discard_flag = true;
    case 'q'
      fprintf('%s: discarding selection, quitting now\n',...
        mfilename);
      discard_flag = true;
      quit_flag = true;
    case '1n'
      fprintf('%s: selection confirmed, enter notes below:\n\n',...
        mfilename);
      notes = input('input notes: ','s');
    case '1nq'
      fprintf('%s: selection confirmed, enter notes below, then will quit\n\n',...
        mfilename);
      notes = input('input notes: ','s');
      quit_flag = true;
    case 'h'
      fprintf('%s: restarting subject...\n',mfilename);
      restart_flag = true;
      display_instructions;
    otherwise
      fprintf('%s: selection input problem, restarting subject...\n',...
        mfilename);
     restart_flag = true;  
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitInfo = remove_rejects(parms,SeriesSelection,VisitInfo)
  rm_flag = ...
    input('remove files for rejected series and update VisitInfo? [y/n] ','s');
  if ~strcmp(rm_flag,'y'), return; end;
  % set DCM_RejectSeries values in Visitinfo
  nvisits = length(SeriesSelection);
  for ind_select=1:nvisits
    VisitID = SeriesSelection(ind_select).VisitID;
    ind_v = find(strcmp(VisitID,{VisitInfo.VisitID}));
    ind_s = find(strcmp(VisitID,{parms.StudyInfo.VisitID}));
    if isempty(ind_v) | isempty(ind_s)
      continue;
    end;
    % delete files for rejected series
    if ~isempty(SeriesSelection(ind_select).RejectSeries) &&...
      any(~ismember(SeriesSelection(ind_select).RejectSeries,VisitInfo(ind_v).DCM_RejectSeries))

      if SeriesSelection(ind_select).RejectSeries == 0; continue; end;
      fprintf('%s: removing files for %s...\n',mfilename,VisitID);
      for j=1:length(ind_s)
        for d=1:length(parms.rmdirlist)
          dirname = parms.rmdirlist{d};
          ContainerDir = parms.StudyInfo(ind_s(j)).(dirname);
          if isempty(ContainerDir), continue; end;
          ContainerPath = sprintf('%s/%s',parms.RootDirs.(dirname),ContainerDir);
          if ~exist(ContainerPath,'dir'), continue; end;
          cmds = {};
          switch dirname
            case 'raw'
              fnamelist = {...
                sprintf('%s/ContainerInfo.mat',ContainerPath)...
                sprintf('%s/SeriesInfo.csv',ContainerPath)...
              };
            case 'proc_dti'
              fnamelist = {ContainerPath};
          end;
          for f=1:length(fnamelist)
            fname = fnamelist{f};
            dlist = dir(fname);
            if ~isempty(dlist)
              cmds{end+1} = ['rm -r ' fname];
            end;
          end;
          %% todo: is there a way to speed this up?
          for c=1:length(cmds)
            cmd = cmds{c};
            if parms.really_rm_flag
               [s,r] = unix(cmd);
              if s
                error('cmd %s failed:\n%s',cmd,r);
              end;
            else
              fprintf('%s\n',cmd);
            end;
          end;
        end;
      end;
      % update visitinfo
      VisitInfo(ind_v).DCM_RejectSeries = ...
        union(VisitInfo(ind_v).DCM_RejectSeries,SeriesSelection(ind_select).RejectSeries);
    end
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitInfo = sort_VisitInfo(parms,VisitInfo)
  VisitIDs = {VisitInfo.VisitID};
  [tmp,ind_sort] = sort(VisitIDs);
  VisitInfo = VisitInfo(ind_sort);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_VisitInfo(parms,VisitInfo)
  tmpinfo = squeeze(struct2cell(VisitInfo))';
  tmpnames = fieldnames(VisitInfo)';
  ind_SubjID = find(strcmp(tmpnames,'SubjID'));
  ind_VisitID = find(strcmp(tmpnames,'VisitID'));
  ind_other = setdiff([1:length(tmpnames)],[ind_SubjID,ind_VisitID]);
  new_order = [ind_SubjID,ind_VisitID,ind_other];
  tmpinfo = tmpinfo(:,new_order);
  tmpnames = tmpnames(new_order);
  mmil_write_csv(parms.fname_VisitInfo,tmpinfo,'col_labels',tmpnames);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subplot_tight(n,m,i)
% from: http://www.briandalessandro.com/blog/how-to-make-a-borderless-subplot-of-images-in-matlab/
  [c,r] = ind2sub([m n], i);
  subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 1/n])
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



