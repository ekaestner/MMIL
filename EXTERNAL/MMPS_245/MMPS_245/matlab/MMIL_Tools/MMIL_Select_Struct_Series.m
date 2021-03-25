function VisitInfo = MMIL_Select_Struct_Series(ProjID,varargin)
%function VisitInfo = MMIL_Select_Struct_Series(ProjID,[options])
% 
% Purpose: visually compare and identify usable structural MRI series
%   e.g. MPR, FLASHhi
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%
% Optional Input:
%  'type': series type as named in ContainerInfo.SeriesInfo
%     MPR, FLASHhi, FLASHlo, XetaT2, MEDIChi, MEDIClo, B1HC, B1BC, GEB1CAL
%     {default = 'MPR'}
%  'infix': string attached to file name indicating processing steps
%     e.g. 'nu', 'uw_nu'
%     {default = []}
%  'ext': file extension
%     e.g. '.mgh','.mgz'
%     {default = '.mgz'}
%  'multi_flag': [0|1] only review data for visits with > 1 series
%     otherwise review all data
%     {default = 1}
%  'attempt_flag': [0|1] review data for multiple attempts together
%     VisitInfo must contain SubjID
%     {default = 0}
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
%     {default = 0}
%  'imagine_flag': [0|1] use imagine to allow browsing of images
%     otherwise display static montages with multiple slices
%     {default = 0}
%
% Created:  09/25/13 by Matt Erhart
% Prev Mod: 04/27/17 by Don Hagler
% Last Mod: 09/21/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  % sort rows by VisitID
  VisitInfo = sort_VisitInfo(VisitInfo);
  % write VisitInfo csv file
  if parms.really_rm_flag
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
    'type','MPR',{'MPR','FLASHhi','FLASHlo',...
                  'XetaT2','MEDIChi','MEDIClo','B1HC','B1BC','GEB1CAL'},...
    'infix',[],[],...
    'ext','mgz',{'mgh','mgz'},...
    'multi_flag',true,[false true],...
    'attempt_flag',true,[false true],...
    'attempt_window',90,[1,1e5],...
    'add_missing_flag',true,[false true],...
    'SubjID_pattern','(?<SubjID>^[^_]+)_\w+',[],...
    'rm_rejects_flag',false,[false true],...
    'imagine_flag',false,[false true],...
  ... % hidden
    'really_rm_flag',true,[false true],...
    'img_dir','Struct_SeriesSelect',[],...
    'img_size',[6 4],[],...
    'img_dpi',200,[],...
    'img_titles',{'EXCLUDE','KEEP'},[],...
    'img_crop',[0.1,0.9],[0,1],...
    'img_fract_slice',0.55,[],...
    'img_fract_slices',[0.3:.1:0.7],[],...
    'img_title_fractx',0.05,[0,1],...
    'img_title_fracty',0.05,[0,1],...
    'RootDirs',[],[],...
    'StudyInfo',[],[],...
    'VisitIDs',[],[],...
    'SubjIDs',[],[],...
    'user',[],[],...
    'rmdirlist',{'raw','proc','proc_dti','proc_bold','fsurf','fsico'},[],...
    'struct_rmdirlist',{'proc_dti','proc_bold'},[],...
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
  parms.suffix = [parms.suffix '.' parms.ext];
  % set series selection file name
  parms.fname_SeriesSelection = ...
    sprintf('%s/ProjInfo/%s/%s_%s_Selection.csv',...
      parms.RootDirs.home,ProjID,ProjID,parms.type);
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

function [SeriesSelection,VisitInfo] = review_data(parms,SeriesSelection,VisitInfo)
  quit_flag = false;
  nvisits = length(SeriesSelection);

  % get SubjIDs from VisitInfo
  VI_SubjIDs = {VisitInfo.SubjID};

  % disable attempt_flag if all SubjIDs are empty
  if all(cellfun(@isempty,VI_SubjIDs))
    parms.attempt_flag = false;
  end;

  % get VisitIDs from VisitInfo and SeriesSelection
  SS_VisitIDs = {SeriesSelection.VisitID};
  VI_VisitIDs = {VisitInfo.VisitID};

  for i=1:nvisits
    VisitID = SeriesSelection(i).VisitID;
    
    % find matching visit in StudyInfo or skip
    ind_study = find(strcmp(VisitID,parms.VisitIDs));
    if isempty(ind_study), continue; end;
    StudyDate = parms.StudyInfo(ind_study).StudyDate;

    ind_study_attempts = ind_study;
    ind_select_attempts = i;
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
    
    % if multiple attempts, set outstem = concatenated VisitIDs
    outstem = [];
    for j=1:nattempts
      if ~isempty(outstem), outstem = [outstem '_+_']; end;
      outstem = [outstem parms.StudyInfo(ind_study_attempts(j)).VisitID];
    end;

    % skip if selection already made
    if nattempts>1
      % skip if only one study has accepted series or if all are rejected      
      naccept = zeros(nattempts,1);
      nreject = zeros(nattempts,1);
      for j=1:nattempts
        ind_select = ind_select_attempts(j);
        naccept(j) = length(SeriesSelection(ind_select).AcceptSeries);
        nreject(j) = length(SeriesSelection(ind_select).RejectSeries);
      end;
      if length(find(naccept))==1 ||...
         sum(naccept)==0 && length(find(nreject))==nattempts
        continue;
      end;
    else
      if ~isempty(SeriesSelection(i).AcceptSeries) ||...
         ~isempty(SeriesSelection(i).RejectSeries)
        continue;
      end;
    end;

    % find scans for each attempt
    ind_series = [];
    ind_scan = [];
    ind_study_series = [];
    ind_select_series = [];
    for j=1:nattempts
      ind_study = ind_study_attempts(j);
      ind_select = ind_select_attempts(j);
      % check for more than one scan
      if isempty(parms.StudyInfo(ind_study).proc), continue; end;
      ContainerPath = sprintf('%s/%s',...
        parms.RootDirs.proc,parms.StudyInfo(ind_study).proc);
      [ContainerInfo, errcode] = MMIL_Load_ContainerInfo(ContainerPath);
      if errcode, continue; end;
      nscans = ContainerInfo.([parms.type '_cntr']);
      if nscans == 0; continue; end;
      ind_series_tmp = [ContainerInfo.ScanInfo.MPR.SeriesIndex]';
      ind_scan_tmp = [1:nscans]';
      % exclude scans that have already been rejected
      ind_reject = SeriesSelection(ind_select).RejectSeries;
      if ~isempty(ind_reject)
        [ind_series_tmp,ind_keep] = setdiff(ind_series_tmp,ind_reject);
        ind_scan_tmp = ind_scan_tmp(ind_keep);
        nscans = length(ind_scan_tmp);
      end;
      ind_series = cat(1,ind_series,ind_series_tmp);
      ind_scan = cat(1,ind_scan,ind_scan_tmp);
      ind_study_series = ...
        cat(1,ind_study_series,ind_study*ones(nscans,1));
      ind_select_series = ...
        cat(1,ind_select_series,ind_select*ones(nscans,1));
    end;
    nscans = length(ind_scan);
    if nscans == 0, continue; end;
    if nscans == 1 && parms.multi_flag
      SeriesSelection(ind_select_series).AcceptSeries = ind_series;
      % update SeriesSelection csv file
      write_selection(parms.fname_SeriesSelection,SeriesSelection,outstem);
      continue;
    end;

    % load data
    fprintf('%s: loading data for %s...\n',mfilename,outstem);
    vol_PRS = {};
    all_exist_flag = 1;
    for j=1:nscans 
      ind_study_tmp = ind_study_series(j);
      ind_series_tmp = ind_series(j);
      ind_scan_tmp = ind_scan(j);
      tmp_VisitID = parms.StudyInfo(ind_study_tmp).VisitID;
      ContainerPath = sprintf('%s/%s',...
        parms.RootDirs.proc,parms.StudyInfo(ind_study_tmp).proc);
      fname_data = sprintf('%s/%s%d%s',...
        ContainerPath,parms.type,ind_scan_tmp,parms.suffix);
      if ~exist(fname_data)
        fprintf('%s: WARNING: for VisitID %s, file %s not found\n',...
          mfilename,tmp_VisitID,fname_data);
        all_exist_flag = 0;
        break;
      end;
      try
        [vol,M] = fs_load_mgh(fname_data,[],1); % first frame only
        vol_PRS(j) = {fs_reorient(vol,M,'PRI')};
        clear vol M;
      catch me
        fprintf('%s: WARNING: for VisitID %s, failed to load %s:\n%s\n',...
          mfilename,tmp_VisitID,fname_data,me.message);
        all_exist_flag = 0;
        break;
      end
    end
    if ~all_exist_flag, continue; end;

    % view and select
    tic
    if parms.imagine_flag
      imagine(vol_PRS{:});
    else
      figure(1);
      show_images(vol_PRS,parms);
    end;
    result = [];
    restart_flag = true;
    discard_flag = false;
    while restart_flag || isempty(result)
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
        % make image to summarize selection
        figure(2); clf;
        save_image(vol_PRS,result,outstem,parms);
        % confirm QC, quit, add notes
        [quit_flag,restart_flag,discard_flag,notes] = ...
          check_result(result,quit_flag,restart_flag,discard_flag);
        for j=1:nattempts
          ind_select = ind_select_attempts(j);
          SeriesSelection(ind_select).notes = notes;
        end;
        figure(2); clf;
        if discard_flag
          if quit_flag, return; end;
          result = [];
          break;
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
    clear vol_PRS;
    if discard_flag, continue; end;

    % record AcceptSeries and RejectSeries numbers
    for j=1:nattempts
      ind_select = ind_select_attempts(j);
      tmp_accept = ind_series(result==1 & ind_select_series==ind_select);
      if ~isempty(tmp_accept)
        SeriesSelection(ind_select).AcceptSeries = tmp_accept;
      else
        SeriesSelection(ind_select).AcceptSeries = [];
      end;
      tmp_reject = union(...
        SeriesSelection(ind_select).RejectSeries,...
        ind_series(result==0 & ind_select_series==ind_select));
      if ~isempty(tmp_reject)
        SeriesSelection(ind_select).RejectSeries = tmp_reject;
      else
        SeriesSelection(ind_select).RejectSeries = [];
      end;
    end;
    
    % update SeriesSelection csv file
    write_selection(parms.fname_SeriesSelection,SeriesSelection,outstem);
    if quit_flag; return; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_selection(fname_SeriesSelection,SeriesSelection,outstem)
  mmil_write_csv(fname_SeriesSelection,...
    squeeze(struct2cell(SeriesSelection))',...
    'col_labels',fieldnames(SeriesSelection)')
  fprintf('%s: wrote selection for %s\n',mfilename,outstem)
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_images(vols,parms)
  nscans = length(vols);
  nslices = length(parms.img_fract_slices);
  figure(1); clf;
  i = 1;
  for j=1:nscans
    [nx,ny,nz] = size(vols{j});
    for k=1:nslices
      img_slice = round(parms.img_fract_slices(k)*ny);
      x0 = round(parms.img_crop(1)*nx);
      x1 = round(parms.img_crop(2)*nx);
      z0 = round(parms.img_crop(1)*nz);
      z1 = round(parms.img_crop(2)*nz);
      subplot_tight(nscans,nslices,i);
      imagesc(flipud(rot90(squeeze(vols{j}(x0:x1,img_slice,z0:z1)))));
      axis off;
      colormap(gray);
      tstr = sprintf('%s %d slice %d',parms.type,j,img_slice);
      x = parms.img_title_fractx*nx;
      y = parms.img_title_fracty*nz;
      text(x,y,tstr,'color','w','fontsize',12);
      i = i + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_image(vols,result,outstem,parms)
  nscans = length(vols);  
  for j=1:nscans
    [nx,ny,nz] = size(vols{j});
    img_slice = round(parms.img_fract_slice*ny);
    subplot(1,nscans,j);
    imagesc(flipud(rot90(squeeze(vols{j}(:,img_slice,:)))));
    axis off;
    colormap(gray);
    tstr = sprintf('%s %d %s',parms.type,j,parms.img_titles{result(j)+1});
    title(tstr);
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
  for i=1:nvisits
    VisitID = SeriesSelection(i).VisitID;
    ind_v = find(strcmp(VisitID,{VisitInfo.VisitID}));
    ind_s = find(strcmp(VisitID,{parms.StudyInfo.VisitID}));
    if isfield(parms.StudyInfo,'STRUCT_VisitID')
      ind_s = union(ind_s,find(strcmp(VisitID,{parms.StudyInfo.STRUCT_VisitID})));
    end;
    if isempty(ind_v) || isempty(ind_s)
      continue;
    end;
    % delete files for rejected series
    if ~isempty(SeriesSelection(i).RejectSeries) &&...
      isempty(intersect(VisitInfo(ind_v).DCM_RejectSeries,...
                        SeriesSelection(i).RejectSeries))
      fprintf('%s: removing files for %s...\n',mfilename,VisitID);
      for j=1:length(ind_s)
        if strcmp(VisitID,parms.StudyInfo(ind_s(j)).VisitID)
          rmdirlist = parms.rmdirlist;
        else
          rmdirlist = parms.struct_rmdirlist;
        end;
        for d=1:length(rmdirlist)
          dirname = rmdirlist{d};
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
              fnamelist = {...
                sprintf('%s/*T1*',ContainerPath)...
                sprintf('%s/*register*',ContainerPath)...
                sprintf('%s/*tkmedit*',ContainerPath)...
                sprintf('%s/*_res_info.mat',ContainerPath)...
                sprintf('%s/AtlasTrack',ContainerPath)...
                sprintf('%s/*calc',ContainerPath)...
                sprintf('%s/*analysis',ContainerPath)...
                sprintf('%s/exportDTI*',ContainerPath)...
              };
            case 'proc_bold'
              fnamelist = {...
                sprintf('%s/*T1*',ContainerPath)...
                sprintf('%s/*register*',ContainerPath)...
                sprintf('%s/*tkmedit*',ContainerPath)...
                sprintf('%s/*_res_info.mat',ContainerPath)...
                sprintf('%s/*analysis',ContainerPath)...
                sprintf('%s/exportBOLD*',ContainerPath)...
              };
            otherwise
              fnamelist = {ContainerPath};
          end;
          for f=1:length(fnamelist)
            fname = fnamelist{f};
            dlist = dir(fname);
            if ~isempty(dlist)
              cmds{end+1} = ['rm -r ' fname];
            end;
          end;
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
        union(VisitInfo(ind_v).DCM_RejectSeries,SeriesSelection(i).RejectSeries);
    end
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VisitInfo = sort_VisitInfo(VisitInfo)
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



