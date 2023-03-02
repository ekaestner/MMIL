function abcd_radreview(varargin)
%function abcd_radreview([options])
%
% Purpose: identify T1- and T2-weighted images 
%            with acceptable quality and that are protocol compliant
%          and send the corresponding tgz files for review by radiologist
%
% Optional Input: ('key', value pairs)
%   'indir': input directory containing individual series dirs
%     {default = '/space/syn05/1/data/MMILDB/DAL_ABCD_QC/incoming'}
%   'fname_info': spreadsheet containing PC and QC info for each series
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_merged_pcqcinfo.csv'}
%   'fname_enrollment': spreadsheet containing pGUIDs for enrolled subjects only
%     {default = []}
%   'fname_review': name of file created to record date when files were sent
%     for radiologist review
%     {default = 'ABCD_radreview.csv'}
%   'logfile': name of file created to log messages
%     {default = []}
%   'local_user': user name on local server
%     {default = 'ABCDRadReview'}
%   'local_ip': IP address of local server
%     {default = '169.228.56.166'}
%   'local_indir': input directory on local server
%     {default = '/data/home/acquisition_sites'}
%   'local_subdir': input subdirectory directory on local server
%       relative to site directory within local_indir
%     {default = 'fiona/outbox/'}
%   'remote_user': user name on remote server
%     {default = 'abcd'}
%   'remote_ip': IP address of remote server
%     {default = 'researchradiology.com'}
%   'remote_outdir': output directory on remote server
%     {default = '/uploads'}
%   'qcflag': [0|1] only send data with good QC
%     {default = 1}
%   'allflag': [0|1] send all image types (T1 and T2) or send nothing
%     if 0, send whichever image types are available
%     {default = 1}
%   'fname_deplck': file name of lock file of dependent job.
%     {default = []}
%
% Created:  11/26/16 by Don Hagler
% Last Mod: 05/19/17 by Feng Xue 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin);

%check lock files
if ~isempty(parms.fname_deplck)
  while exist(parms.fname_deplck,'file')
    fprintf('%s\n','lock files exist, waiting for previous process to finish.');
    pause(30);
  end;
end;


% create log file and display starting message
parms = create_log(parms);

% initialize rad review file
review_info = get_review_info(parms);

% load series info
series_info = get_series_info(parms);

% exclude non-enrolled subjects
if ~isempty(parms.fname_enrollment)
  series_info = check_enrollment_info(series_info,parms);
end;

% exclude series for events that have already been reviewed
series_info = filter_series_info(series_info,review_info,parms);
if ~isempty(series_info)
  tempfname = mmil_tempfname;
  status = mkdir(tempfname);
  if ~status
    error('%s: can not create temp folder: %s',mfilename,tempfname);
    return
  end
  parms.tempfname = tempfname;
  % loop over events, sending files for radiologist review
  send_files(series_info,parms);
  status = rmdir(tempfname, 's');
  if ~status
    error('%s: can not remove temp folder: %s',mfilename,tempfname);
    return
  end
end

% close log file
if parms.flog>0, fclose(parms.flog); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','/space/syn05/1/data/MMILDB/DAL_ABCD_QC/incoming',[],...
    'fname_info','/home/abcddaic/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_combined_pcinfo.csv',[],...
    'fname_enrollment',[],[],...
    'fname_review','radreview.csv',[],...
    'logfile',[],[],...
    'local_user','ABCDRadReview',[],...
    'local_ip','169.228.56.166',[],...
    'local_indir','/data/home/acquisition_sites',[],...
    'local_subdir','fiona/outbox/',[],...
    'remote_user','abcd',[],...
    'remote_ip','researchradiology.com',[],...
    'remote_outdir','/uploads',[],...
    'qcflag',true,[false true],...
    'allflag',true,[false true],...
    'fname_deplck',[],[],...
    ...
    'series_types',{'t1','t2'},[],...
    'date_pat','yyyy-mm-dd',[],...
  });
  if ~exist(parms.indir)
    error('input directory %s not found',parms.indir);
  end;
  if ~exist(parms.fname_info,'file')
    error('info file %s not found',parms.fname_info);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_log(parms)
  % create logfile
  if ~isempty(parms.logfile)
    logdir = fileparts(parms.logfile);
    mmil_mkdir(logdir);
    parms.flog = fopen(parms.logfile,'a');
    if parms.flog==-1
      error('failed to open logfile %s for writing\n',parms.logfile);
    end;
  else
    parms.flog = -1;
  end;
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf('%s: starting abcd rad review from %s at %s...\n',...
    mfilename,parms.indir,datestr(now));
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  if parms.flog>0
    fprintf(parms.flog,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(parms.flog,'%s: starting abcd rad review from %s at %s...\n',...
      mfilename,parms.indir,datestr(now));
    fprintf(parms.flog,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function review_info = get_review_info(parms)
  % create column headers if file does not already exist
  if ~exist(parms.fname_review,'file')
    fprintf('%s: creating %s...\n',mfilename,parms.fname_review);
    if parms.flog>0
      fprintf(parms.flog,'%s: creating %s...\n',mfilename,parms.fname_review);
    end;
    fid = fopen(parms.fname_review,'w');
    if fid < 0
      error('failed to open file %s for writing',parms.fname_review);
    end;
    fprintf(fid,'"pguidevent","id_redcap","redcap_event_name","mrif_t1_sent","mrif_t2_sent","mrif_daic2rad_send","mrif_ready_review"\n');
    fclose(fid);
  end;
  fprintf('%s: reading %s...\n',mfilename,parms.fname_review);
  if parms.flog>0
    fprintf(parms.flog,'%s: reading %s...\n',mfilename,parms.fname_review);
  end;
  review_info = mmil_csv2struct(parms.fname_review);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_info = get_series_info(parms)
  fprintf('%s: reading %s...\n',mfilename,parms.fname_info);
  if parms.flog>0
    fprintf(parms.flog,'%s: reading %s...\n',mfilename,parms.fname_info);
  end;
  % read PC and QC info
  series_info = abcd_load_csv(parms.fname_info);
  
  % exclude empty pGUID
  if isfield(series_info,'pGUID')
    i_valid = find(~cellfun(@isempty,{series_info.pGUID}) |...
                   ~cellfun(@isnumeric,{series_info.pGUID}));
    series_info = series_info(i_valid);
    %% exclude non-NDAR pGUID
    %i_valid = find(~cellfun(@isempty,regexp({series_info.pGUID},'^NDAR')));
    % include non-Phantom pGUID
    i_valid = find(cellfun(@isempty,regexp({series_info.pGUID},'Phantom')));
  else
    i_valid = find(~cellfun(@isempty,{series_info.id_redcap}) |...
                   ~cellfun(@isnumeric,{series_info.id_redcap}));
    series_info = series_info(i_valid);
    %% exclude non-NDAR id_redcap
    %i_valid = find(~cellfun(@isempty,regexp({series_info.id_redcap},'^NDAR')));
    % include non-Phantom pGUID
    i_valid = find(cellfun(@isempty,regexp({series_info.id_redcap},'Phantom')));
  end
  series_info = series_info(i_valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_info = check_enrollment_info(series_info,parms)
  if isempty(parms.fname_enrollment), return; end;
  fprintf('%s: reading %s...\n',mfilename,parms.fname_enrollment);
  if parms.flog>0
    fprintf(parms.flog,'%s: reading %s...\n',mfilename,parms.fname_info);
  end;
  % read enrollment info
  enroll_info = abcd_load_csv(parms.fname_enrollment);
  % find series with pGUIDs that are members of enrollment pGUIDs
  enroll_pGUIDs = unique({enroll_info.pGUID});
  if isfield(series_info,'pGUID')
    i_valid = ismember({series_info.pGUID},enroll_pGUIDs);
  else
    i_valid = ismember({series_info.id_redcap},enroll_pGUIDs);
  end
  fprintf('%s: %d enrolled subjects...\n',...
    mfilename,length(enroll_pGUIDs));
  fprintf('%s: %d series from enrolled subjects...\n',...
    mfilename,length(i_valid));
  if parms.flog>0
    fprintf(parms.flog,'%s: %d enrolled subjects to send...\n',...
      mfilename,length(enroll_pGUIDs));
    fprintf(parms.flog,'%s: %d series from enrolled subjects to send...\n',...
      mfilename,length(i_valid));
  end;
  series_info = series_info(i_valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_info = filter_series_info(series_info,review_info,parms)
  fprintf('%s: filtering by series type...\n',mfilename);
  if parms.flog>0
    fprintf(parms.flog,'%s: filtering by series type...\n',mfilename);
  end;

  % filter by series type
  SeriesTypes = lower({series_info.SeriesType});
  ind_keep = find(ismember(SeriesTypes,lower(parms.series_types)));
  series_info = series_info(ind_keep);
  % filter series that have been sent for review
  if isempty(review_info), return; end;
  fprintf('%s: excluding events that have already been sent for review...\n',...
    mfilename);
  if parms.flog>0
    fprintf(parms.flog,'%s: excluding events that have already been sent for review...\n',...
      mfilename);
  end;
  % get combined pGUID-EventNames for each series
  if isfield(series_info,'pguidevent')
    s_pGUID_EventNames={series_info.pguidevent};
  else
    s_pGUIDs = {series_info.pGUID};
    s_EventNames = {series_info.EventName};
    s_pGUID_EventNames = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                                 s_pGUIDs,s_EventNames,...
                                 'UniformOutput',false);
  end
  % get combined pGUID-EventNames for each reviewed event
%  r_pGUIDs = {review_info.};
%  r_EventNames = {review_info.EventName};
%  r_pGUID_EventNames = cellfun(@(x,y) sprintf('%s_%s',x,y),...
%                               r_pGUIDs,r_EventNames,...
%                               'UniformOutput',false);
                               
  % exclude series from reviewed events
  %ind_exclude = find(ismember(s_pGUID_EventNames,r_pGUID_EventNames));
  ind_exclude = find(ismember(s_pGUID_EventNames,{review_info.pguidevent}));
  ind_keep = setdiff(1:length(s_pGUID_EventNames),ind_exclude);
  series_info = series_info(ind_keep);
  fprintf('%s: %d series to send\n',...
    mfilename,length(series_info));
  if parms.flog>0
    fprintf(parms.flog,'%s: %d series to send\n',...
      mfilename,length(series_info));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function send_files(series_info,parms)
  % get unique subject IDs
  if isfield(series_info,'pGUID')
    pGUIDs = unique({series_info.pGUID});
  else
    pGUIDs = unique({series_info.id_redcap});
  end
  % loop over unique subject IDs
  for i=1:length(pGUIDs)
    pGUID = pGUIDs{i};
    % identify all series for this subject
    subj_info = select_subj_series(pGUID,series_info,parms);
    if isempty(subj_info), continue; end;
    % loop over events
    if isfield(subj_info,'EventName')
      EventNames = {subj_info.EventName};
      site = subj_info(1).SiteName;
    else
      EventNames = {subj_info.redcap_event_name};
      site = subj_info(1).site;
    end
    uniq_EventNames = unique(EventNames);
    errcode = 0;
    for j=1:length(uniq_EventNames)
      EventName = uniq_EventNames{j};
      i_event = find(strcmp(EventNames,EventName));
      event_info = subj_info(i_event);
      StudyDate = event_info(1).StudyDate;
%      if isfield(event_info,'SiteName')
%        site = event_info(1).SiteName;
%      else
%        site = event_info(1).site;
%      end
      % check if have T1 and T2
      series_types = upper({event_info.SeriesType});
      i_T1 = find(strcmp(series_types,'T1'));
      if isempty(i_T1)
        fprintf('%s: WARNING: missing or bad T1 for %s %s %d from %s\n',...
          mfilename,pGUID,EventName,StudyDate,site);
        if parms.flog>0
          fprintf(parms.flog,'%s: WARNING: missing or bad T1 for %s %s %d from %s\n',...
            mfilename,pGUID,EventName,StudyDate,site);
        end;
        if parms.allflag
          continue;
        end;
      end;
      i_T2 = find(strcmp(series_types,'T2'));
      if isempty(i_T2)
        fprintf('%s: WARNING: missing or bad T2 for %s %s %d from %s\n',...
          mfilename,pGUID,EventName,StudyDate,site);
        if parms.flog>0
          fprintf(parms.flog,'%s: WARNING: missing or bad T2 for %s %s %d from %s\n',...
            mfilename,pGUID,EventName,StudyDate,site);
        end;
        if parms.allflag
          continue;
        end;
      end;
      % skip if both scan types are missing
      if isempty(i_T1) && isempty(i_T2)
        continue;
      elseif ~isempty(i_T1) && ~isempty(i_T2)
        % use the last good scan of each
        event_info = [event_info(i_T1(end)),event_info(i_T2(end))];
      elseif ~isempty(i_T1)
        event_info = [event_info(i_T1(end))];
      elseif ~isempty(i_T2)
        event_info = [event_info(i_T2(end))];
      end;
      nseries = length(event_info);
      series_flags = zeros(2,1);
      if ~isempty(i_T1), series_flags(1) = 1; end;
      if ~isempty(i_T2), series_flags(2) = 1; end;
      %if strcmp(pGUID,'NDAR_INVKM38D8JZ'), keyboard; end; 
      %if strcmp(pGUID,'NDAR_INVU8KM6L2C'), keyboard; end; 
      for k=1:nseries
        % identify tgz file in incoming (need file name and site name)
        SUID = event_info(k).StudyInstanceUID;
        SeUID = event_info(k).SeriesInstanceUID;
        if isfield(subj_info,'EventName')
          site = event_info(k).SiteName;
        else
          site = event_info(k).site;
        end

        cmd = sprintf('(ssh %s@%s find %s/%s/fiona/outbox/ -name ''NDAR_INV*%s_%s.tgz\''>/dev/tty) >& /dev/null',parms.local_user,parms.local_ip,parms.local_indir,site,SUID,SeUID);
        [s,flist] = unix(cmd);
        %flist = dir(sprintf('%s/%s/*%s_%s*.tgz',parms.indir,site,SUID,SeUID));
        if isempty(flist)
          %if strcmp(pGUID,'NDAR_NDAR_INV569JK9GU'), keyboard; end; 
          fprintf('%s: WARNING: file for %s %s %s with SUID %s and SeUID %s not found in %s/%s\n',...
            mfilename,pGUID,EventName,StudyDate,SUID,SeUID,parms.indir,site);
          if parms.flog>0
            fprintf(parms.flog,'%s: WARNING: file for %s %s %s with SUID %s and SeUID %s not found in %s/%s\n',...
              mfilename,pGUID,EventName,StudyDate,SUID,SeUID,parms.indir,site);
          end;
          errcode = 1;
          break;
        end;
        [pathstr,fname,ext] = fileparts(flist);
        
        fprintf('%s: fetching file %s to %s...\n',mfilename,fname,...
          parms.tempfname);
        if parms.flog>0
          fprintf(parms.flog,'%s: fetching file %s to %s...\n',mfilename,fname,...
            parms.tempfname);
        end;

        cmd = sprintf('scp %s@%s:%s/%s.tgz %s/',...
          parms.local_user,parms.local_ip,...
          pathstr,fname,parms.tempfname...
          );
        [s,r] = unix(cmd);

        % check for error
        %if s, error('cmd %s failed:\n%s',cmd,r); end;
        if s
          fprintf('%s: failed to fetch file %s and save in %s...\n',mfilename,fname,parms.tempfname);
          if parms.flog>0
            fprintf(parms.flog,'%s: failed to fetch file %s and save in %s...\n',mfilename,fname,...
              parms.tempfname);
          end;
          errcode = 1;
          break;
        end

        cmd = sprintf('tar xf %s/%s.tgz -C %s',parms.tempfname,fname,parms.tempfname);
        [s,r] = unix(cmd);

        % check for error
        %if s, error('cmd %s failed:\n%s',cmd,r); end;
        if s
          fprintf('%s: failed to extract file %s/%s...\n',mfilename,parms.tempfname,fname);
          if parms.flog>0
            fprintf(parms.flog,'%s: failed to extract file %s/%s...\n',mfilename,parms.tempfname,fname);
          end;

          cmd = sprintf('rm -rf %s/*',parms.tempfname);
          [s,r] = unix(cmd);

          errcode = 1;
          break;
        end

        cmd = sprintf('dcm2niix -o %s %s',parms.tempfname,parms.tempfname);
        [s,result_conversion] = unix(cmd);

        % check for error
        %if s, error('cmd %s failed:\n%s',cmd,r); end;
        if s
          fprintf('%s: failed to convert DICOM files under %s...\n',mfilename,parms.tempfname);
          if parms.flog>0
            fprintf(parms.flog,'%s: failed to convert DICOM files under %s...\n',mfilename,parms.tempfname);
          end;

          cmd = sprintf('rm -rf %s/*',parms.tempfname);
          [s,r] = unix(cmd);

          errcode = 1;
          break;
        end

        cmd = sprintf('ls %s/*.nii|wc -l',parms.tempfname);
        [s,r] = unix(cmd);

        % check for error
        %if s, error('cmd %s failed:\n%s',cmd,r); end;
        if s
          fprintf('%s: DICOM conversion failure: no nii file found in %s...\n',mfilename,parms.tempfname);
          if parms.flog>0
            fprintf(parms.flog,'%s: DICOM conversion failure: no nii file found in %s...\n',mfilename,parms.tempfname);
          end;

          cmd = sprintf('rm -rf %s/*',parms.tempfname);
          [s,r] = unix(cmd);

          errcode = 1;
          break;
        end
        if str2num(r(1:end-1)) >1
          fprintf('%s: DICOM conversion failure: more than one nii file found in %s...\n',mfilename,parms.tempfname);
          if parms.flog>0
            fprintf(parms.flog,'%s: DICOM conversion failure: more than one nii file found in %s...\n',mfilename,parms.tempfname);
          end;

          cmd = sprintf('rm -rf %s/*',parms.tempfname);
          [s,r] = unix(cmd);

          errcode = 1;
          break;
        end
        
        pos_found = findstr(result_conversion,'Found ');
        r_rest = result_conversion(pos_found:end);
        pos_dicom = findstr(r_rest,'DICOM');
        total_dicom_count = str2num(r_rest(7:pos_dicom(1)-1));
        pos_convert = findstr(r_rest,'Convert ');
        %r_rest = result_conversion(pos_found:end);
        %pos_dicom = findstr(r_rest,'DICOM');
        converted_dicom_count = str2num(r_rest(pos_convert+8:pos_dicom(2)-1));
        if total_dicom_count ~= converted_dicom_count
          fprintf('%s: DICOM conversion failure: total count of DICOM files does not equal converted DICOM files in %s...\n',mfilename,parms.tempfname);
          if parms.flog>0
            fprintf(parms.flog,'%s: DICOM conversion failure: total count of DICOM files does not equal converted DICOM files in %s...\n',mfilename,parms.tempfname);
          end;

          cmd = sprintf('rm -rf %s/*',parms.tempfname);
          [s,r] = unix(cmd);

          errcode = 1;
          break;
        end

        cmd = sprintf('rm -rf %s/*',parms.tempfname);
        [s,r] = unix(cmd);

        % ssh to an ip address as radreview (destination FIONA) and sftp to another ip address
        cmd = sprintf('ssh %s@%s ''echo put %s/%s/%s%s.tgz | sftp -b- %s@%s:%s''',...
          parms.local_user,parms.local_ip,...
          parms.local_indir,site,parms.local_subdir,fname,...
          parms.remote_user,parms.remote_ip,parms.remote_outdir);
        fprintf('%s: sending file %s to %s@%s...\n',mfilename,fname,...
          parms.remote_user,parms.remote_ip);
        if parms.flog>0
          fprintf(parms.flog,'%s: sending file %s to %s@%s...\n',mfilename,fname,...
            parms.remote_user,parms.remote_ip);
        end;
        [s,r] = unix(cmd);

        % check for error
        if s, error('cmd %s failed:\n%s',cmd,r); end;
      end;
      % save date of rad review
      if ~errcode
        save_review_info(pGUID,EventName,series_flags,parms);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subj_info = select_subj_series(pGUID,series_info,parms)
  subj_info = [];
  if isfield(subj_info,'pGUID')
    pGUIDs = {series_info.pGUID};
  else
    pGUIDs = {series_info.id_redcap};
  end;
  i_subj = find(strcmp(pGUIDs,pGUID));
  subj_info = series_info(i_subj);
  if isfield(subj_info,'ABCD_Compliant')
    PC = {subj_info.ABCD_Compliant};
  else
    PC = [];
  end;
  if isfield(subj_info,'Completed')
    CC = {subj_info.Completed};
  else
    CC = [];
  end;
  if isfield(subj_info,'QC')
    QC = {subj_info.QC};
  else
    QC = [];
  end;

  if all([isempty(PC),isempty(QC),isempty(CC)]), return; end;

  % skip subject if series have not been checked yet
  if parms.qcflag
    i_checked = find(~cellfun(@isempty,PC) & ~cellfun(@isempty,QC));
  else
    i_checked = find(~cellfun(@isempty,PC));
  end;
  if isempty(i_checked)
    fprintf('%s: WARNING: no checked series for pGUID %s\n',mfilename,pGUID);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: no checked series for pGUID %s\n',mfilename,pGUID);
    end;
    return;
  end;
  % skip subject if no good series
  if parms.qcflag
    i_good = find(~cellfun(@isempty,regexp(PC(i_checked),'Yes')) &...
                  [CC{i_checked}] &...
                  [QC{i_checked}]);
  else
    i_good = find(cellfun(@isempty,regexp(PC(i_checked),'NA')) &...
                  [CC{i_checked}]);
    if ~isempty(QC)
      i_qc = find(~cellfun(@isempty,QC));
      if ~isempty(i_qc)
        qc_vals = cell2mat(QC(i_qc));
        if any(qc_vals==0)
          fprintf('%s: WARNING: sending data with QC = 0 for pGUID %s\n',mfilename,pGUID);
        end;
      end;
    end;
  end;
  if isempty(i_good)
    fprintf('%s: WARNING: no good series for pGUID %s\n',mfilename,pGUID);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: no good series for pGUID %s\n',mfilename,pGUID);
    end;
    return;
  end;
  subj_info = subj_info(i_checked(i_good));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_review_info(pGUID,EventName,series_flags,parms)
  fid = fopen(parms.fname_review,'a');
  if fid < 0
    error('failed to open file %s for writing',parms.fname_review);
  end;
  ready_flag = all(series_flags);
  %fprintf(fid,'"%s_%s","%s","%s",%d,%d,"%s",%d\n',...
  fprintf(fid,'%s_%s,%s,%s,%d,%d,%s,%d\n',...
    pGUID,EventName,...
    pGUID,EventName,...
    series_flags(1),series_flags(2),...
    datestr(now,parms.date_pat),...
    ready_flag);
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
