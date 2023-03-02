function abcd_check_incoming2radreview(varargin)
%function abcd_check_incoming2radreview([options])
%
% Purpose: check three gaps:
%          1. mrif_unpack_gap: between scan date and unpack date
%          2. mrif_daic2rad_gap: between scan date and send-to-radreview date
%          3. mrif_radreview_gap: between send-to-radreview date and radreviewed date
%
% Optional Input: ('key', value pairs)
%   'indir': input directory
%     {default = 'MetaData/DAL_ABCD_QC'}
%   'outdir': output directory
%     {default = 'MetaData/DAL_ABCD_QC'}
%   'instem': instem
%     {default = 'DAL_ABCD_QC'}
%   'outstem': outstem
%     {default = 'DAL_ABCD_QC'}
%   'fname_info': spreadsheet containing classified incoming info
%     {default = ''}
%   'fname_radreview_info': spreadsheet containing radreview info
%     {default = ''}
%   'fname_enrollment': spreadsheet containing pGUIDs for enrolled subjects only
%     {default = []}
%   'fname_incoming2radreview': name of output file
%     {default = ''}
%   'fname_projinfo': file name of whole project info
%     {default = []}
%   'logfile': name of file created to log messages
%     {default = []}
%   'infix': file suffix of input file
%     containing info about classified incoming info
%     {default = 'combined_incoming_info_classified'}
%   'rad_infix': file suffix of input file
%     containing info about radreview
%     {default = 'radreview_dates'}
%   'outfix : file suffix for combined info file
%     {default = 'incoming2radreview_info'}
%   'qcflag': [0|1] only send data with good QC
%     {default = 1}
%   'allflag': [0|1] send all image types (T1 and T2) or send nothing
%     if 0, send whichever image types are available
%     {default = 1}
%   'forceflag': ignore existing results and run
%     {default = 1}
%   'date_pat': pattern of date time
%     {default = 'yyyy-mm-dd'}
%    'series_types' all valid SeriesTypes
%      {default = {'T1','T2'}
%
% Created:  04/21/17 by Feng Xue
% Last Mod: 08/08/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin);

% create log file and display starting message
parms = create_log(parms);

% initialize rad review file
review_info = abcd_load_csv(parms.fname_radreview_info);

% load series info
series_info = abcd_load_csv(parms.fname_info);


% exclude non-enrolled subjects
if ~isempty(parms.fname_enrollment)
  series_info = check_enrollment_info(series_info,parms);
end;

series_info = filter_series_info(series_info,parms);


if length(series_info)
  incoming2radreview_info=[];
  if ~isfield(series_info,'pguidevent'), series_info = set_pguidevents(series_info); end;
  if ~isfield(review_info,'pguidevent'), review_info = set_pguidevents(review_info); end;
  incoming2radreview_info = check_incoming2radreview(series_info,parms,review_info);
  incoming2radreview_info = mmil_sortstruct(incoming2radreview_info,{'site','id_redcap','redcap_event_name'});
  mmil_struct2csv(incoming2radreview_info,parms.fname_incoming2radreview);
end

% close log file
if parms.flog>0, fclose(parms.flog); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'fname_info','',[],...
    'fname_radreview_info',[],[],...
    'fname_incoming2radreview',[],[],...
    'fname_enrollment',[],[],...
    'fname_projinfo',[],[],...
    'logfile',[],[],...
    'infix','combined_incoming_info_classified',[],...
    'rad_infix','radreview_dates',[],...
    'outfix','incoming2radreview_info',[],...
    'forceflag',true,[false true],...
    'series_types',{'T1','T2'},[],...
    'qcflag',true,[false true],...
    'allflag',true,[false true],...
    'date_pat','yyyy-mm-dd',[],...
  });
  if parms.outdir(1) ~= '/'
     parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir);
  end
  if parms.indir(1) ~= '/', parms.indir=parms.outdir; end;

  if ~exist(parms.indir)
    error('input directory %s not found',parms.indir);
  end;
  if isempty(parms.fname_info)
    parms.fname_info = sprintf('%s/%s_%s.csv',...
      parms.indir,parms.instem,parms.infix);
  end;
  if ~exist(parms.fname_info,'file')
    error('info file %s not found',parms.fname_info);
  end;
  if isempty(parms.fname_radreview_info)
    parms.fname_radreview_info = sprintf('%s/%s_%s.csv',...
      parms.indir,parms.instem,parms.rad_infix);
  end;
  if ~exist(parms.fname_radreview_info,'file')
    error('info file %s not found',parms.fname_radreview_info);
  end;
  if isempty(parms.fname_incoming2radreview)
    parms.fname_incoming2radreview = sprintf('%s/%s_%s.csv',...
      parms.indir,parms.instem,parms.outfix);
  end;
  if isempty(parms.fname_projinfo), parms.fname_projinfo = sprintf('%s/ProjInfo/MMIL_ProjInfo_all.csv',getenv('HOME')); end;
  if ~exist(parms.fname_projinfo,'file'), error('info file %s not found',parms.fname_projinfo); end;
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

function series_info = filter_series_info(series_info,parms)

  % exclude empty pGUID
  if isfield(series_info,'pGUID')
    i_valid = find(~cellfun(@isempty,{series_info.pGUID}) |...
                   ~cellfun(@isnumeric,{series_info.pGUID}));
  else
    i_valid = find(~cellfun(@isempty,{series_info.id_redcap}) |...
                   ~cellfun(@isnumeric,{series_info.id_redcap}));
  end
  series_info = series_info(i_valid);
  % exclude non-NDAR pGUID
  if isfield(series_info,'pGUID')
    i_valid = find(~cellfun(@isempty,regexp({series_info.pGUID},'^NDAR')));
  else
    i_valid = find(~cellfun(@isempty,regexp({series_info.id_redcap},'^NDAR')));
  end
  series_info = series_info(i_valid);

  fprintf('%s: filtering by series type...\n',mfilename);
  if parms.flog>0
    fprintf(parms.flog,'%s: filtering by series type...\n',mfilename);
  end;

  % filter by series type
  SeriesTypes = lower({series_info.SeriesType});
  ind_keep = find(ismember(SeriesTypes,lower(parms.series_types)));
  series_info = series_info(ind_keep);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function incoming2radreview_info=check_incoming2radreview(series_info,parms,review_info)
  % get unique subject IDs
  pguidevents = unique({series_info.pguidevent});
  
  incoming2radreview_info = [];
  incoming2radreview_info_pos=0;

  % load project info
  projinfo = abcd_load_projinfo_all(parms);

  % loop over unique pguidevents
  for i=1:length(pguidevents)
    pguidevent = pguidevents{i};
    % identify all series for this subject
    [subj_info_struct,subj_info_review] = select_subj_series(pguidevent,series_info,review_info,parms);
    if (isempty(subj_info_struct) || isempty(subj_info_review)), continue; end;

    % check if have T1 and T2
    series_types = upper({subj_info_struct.SeriesType});
    i_T1 = find(strcmp(series_types,'T1'));
    i_T2 = find(strcmp(series_types,'T2'));
    % skip if both scan types are missing
    if isempty(i_T1) && isempty(i_T2), continue; end;

    incoming2radreview_info_pos = incoming2radreview_info_pos+1;
    incoming2radreview_info(incoming2radreview_info_pos).pguidevent = pguidevent;
    if isfield(subj_info_struct,'pGUID')
      incoming2radreview_info(incoming2radreview_info_pos).id_redcap = subj_info_struct(1).pGUID;
      incoming2radreview_info(incoming2radreview_info_pos).redcap_event_name = subj_info_struct(1).EventName;
      incoming2radreview_info(incoming2radreview_info_pos).site = subj_info_struct(1).SiteName;
    else
      incoming2radreview_info(incoming2radreview_info_pos).id_redcap = subj_info_struct(1).id_redcap;
      incoming2radreview_info(incoming2radreview_info_pos).redcap_event_name = subj_info_struct(1).redcap_event_name;
      incoming2radreview_info(incoming2radreview_info_pos).site = subj_info_struct(1).site;
    end

    if isempty(subj_info_review(end).mrif_scan_read_dte)
      scan_read_date=[];
      incoming2radreview_info(incoming2radreview_info_pos).mrif_radreview_gap=[];
    else
      scan_read_date=datenum(regexprep(num2str(subj_info_review(end).mrif_scan_read_dte),'_',''),'yyyymmdd');
      incoming2radreview_info(incoming2radreview_info_pos).mrif_radreview_gap=scan_read_date - datenum(regexprep(num2str(subj_info_review(end).mrif_daic2rad_send),'_',''),'yyyymmdd');
    end

    %incoming2radreview_info(incoming2radreview_info_pos).mrif_radreviewdate = datestr(datenum(regexprep(num2str(subj_info_review(end).mrif_daic2rad_send),'_',''),'yyyymmdd'),parms.date_pat);
    if ~isempty(i_T1) && ~isempty(i_T2)
      % use the last good scan of each
      incoming2radreview_info(incoming2radreview_info_pos).mrif_t1_studydate = datestr(datenum(num2str(subj_info_struct(i_T1(end)).StudyDate),'yyyymmdd'),parms.date_pat);
      incoming2radreview_info(incoming2radreview_info_pos).mrif_t2_studydate = datestr(datenum(num2str(subj_info_struct(i_T2(end)).StudyDate),'yyyymmdd'),parms.date_pat);

      if isfield(subj_info_struct,'pGUID')
        incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,subj_info_struct(i_T1(end)).fname_pc_json,...
                                                  subj_info_struct(i_T1(end)).pGUID,subj_info_struct(i_T1(end)).EventName,subj_info_struct(i_T1(end)).StudyInstanceUID,...
                                                  subj_info_struct(i_T1(end)).SeriesInstanceUID,subj_info_struct(i_T1(end)).SiteName,'T1',parms,projinfo);

        incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,subj_info_struct(i_T2(end)).fname_pc_json,...
                                                  subj_info_struct(i_T2(end)).pGUID,subj_info_struct(i_T2(end)).EventName,subj_info_struct(i_T2(end)).StudyInstanceUID,...
                                                  subj_info_struct(i_T2(end)).SeriesInstanceUID,subj_info_struct(i_T2(end)).SiteName,'T2',parms,projinfo);
      else
        incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,subj_info_struct(i_T1(end)).json,...
                                                  subj_info_struct(i_T1(end)).id_redcap,subj_info_struct(i_T1(end)).redcap_event_name,subj_info_struct(i_T1(end)).StudyInstanceUID,...
                                                  subj_info_struct(i_T1(end)).SeriesInstanceUID,subj_info_struct(i_T1(end)).site,'T1',parms,projinfo);

        incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,subj_info_struct(i_T2(end)).json,...
                                                  subj_info_struct(i_T2(end)).id_redcap,subj_info_struct(i_T2(end)).redcap_event_name,subj_info_struct(i_T2(end)).StudyInstanceUID,...
                                                  subj_info_struct(i_T2(end)).SeriesInstanceUID,subj_info_struct(i_T2(end)).site,'T2',parms,projinfo);
      end
      gap_T1=datenum(regexprep(num2str(subj_info_review(end).mrif_daic2rad_send),'_',''),'yyyymmdd') - datenum(num2str(subj_info_struct(i_T1(end)).StudyDate),'yyyymmdd');
      gap_T2=datenum(regexprep(num2str(subj_info_review(end).mrif_daic2rad_send),'_',''),'yyyymmdd') - datenum(num2str(subj_info_struct(i_T2(end)).StudyDate),'yyyymmdd');
      if gap_T1>=gap_T2
        incoming2radreview_info(incoming2radreview_info_pos).mrif_daic2rad_gap=gap_T1;
      else
        incoming2radreview_info(incoming2radreview_info_pos).mrif_daic2rad_gap=gap_T2;
      end

      if ~isempty(incoming2radreview_info(incoming2radreview_info_pos).mrif_t1_unpackdate)
        gap_T1=datenum(incoming2radreview_info(incoming2radreview_info_pos).mrif_t1_unpackdate,parms.date_pat) - datenum(num2str(subj_info_struct(i_T1(end)).StudyDate),'yyyymmdd');
      else
        gap_T1=[];
      end
      if ~isempty(incoming2radreview_info(incoming2radreview_info_pos).mrif_t2_unpackdate)
        gap_T2=datenum(incoming2radreview_info(incoming2radreview_info_pos).mrif_t2_unpackdate,parms.date_pat) - datenum(num2str(subj_info_struct(i_T2(end)).StudyDate),'yyyymmdd');
      else
        gap_T2=[];
      end
      if any([isempty(gap_T1),isempty(gap_T2)])
        incoming2radreview_info(incoming2radreview_info_pos).mrif_unpack_gap=[];
      else
        if gap_T1>=gap_T2
          incoming2radreview_info(incoming2radreview_info_pos).mrif_unpack_gap=gap_T1;
        else
          incoming2radreview_info(incoming2radreview_info_pos).mrif_unpack_gap=gap_T2;
        end
      end

    elseif ~isempty(i_T1)
      incoming2radreview_info(incoming2radreview_info_pos).mrif_t1_studydate = datestr(datenum(num2str(subj_info_struct(i_T1(end)).StudyDate),'yyyymmdd'),parms.date_pat);
      if isfield(series_info,'pGUID')
        incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,subj_info_struct(i_T1(end)).fname_pc_json,...
                                                  subj_info_struct(i_T1(end)).pGUID,subj_info_struct(i_T1(end)).EventName,subj_info_struct(i_T1(end)).StudyInstanceUID,...
                                                  subj_info_struct(i_T1(end)).SeriesInstanceUID,subj_info_struct(i_T1(end)).SiteName,'T1',parms,projinfo);
      else
        incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,subj_info_struct(i_T1(end)).json,...
                                                  subj_info_struct(i_T1(end)).id_redcap,subj_info_struct(i_T1(end)).redcap_event_name,subj_info_struct(i_T1(end)).StudyInstanceUID,...
                                                  subj_info_struct(i_T1(end)).SeriesInstanceUID,subj_info_struct(i_T1(end)).site,'T1',parms,projinfo);
      end
      gap_T1=datenum(regexprep(num2str(subj_info_review(end).mrif_daic2rad_send),'_',''),'yyyymmdd') - datenum(num2str(subj_info_struct(i_T1(end)).StudyDate),'yyyymmdd');
      incoming2radreview_info(incoming2radreview_info_pos).mrif_daic2rad_gap=gap_T1;
      
      if ~isempty(incoming2radreview_info(incoming2radreview_info_pos).mrif_t1_unpackdate)
        gap_T1=datenum(incoming2radreview_info(incoming2radreview_info_pos).mrif_t1_unpackdate,parms.date_pat) - datenum(num2str(subj_info_struct(i_T1(end)).StudyDate),'yyyymmdd'); 
      else
        gap_T1=[];
      end
      incoming2radreview_info(incoming2radreview_info_pos).mrif_unpack_gap=gap_T1;
    elseif ~isempty(i_T2)
      incoming2radreview_info(incoming2radreview_info_pos).mrif_t2_studydate = datestr(datenum(num2str(subj_info_struct(i_T2(end)).StudyDate),'yyyymmdd'),parms.date_pat);
      if isfield(subj_info_struct,'pGUID')

        incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,subj_info_struct(i_T2(end)).fname_pc_json,...
                                                  subj_info_struct(i_T2(end)).pGUID,subj_info_struct(i_T2(end)).EventName,subj_info_struct(i_T2(end)).StudyInstanceUID,...
                                                  subj_info_struct(i_T2(end)).SeriesInstanceUID,subj_info_struct(i_T2(end)).SiteName,'T2',parms,projinfo);
      else

        incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,subj_info_struct(i_T2(end)).json,...
                                                  subj_info_struct(i_T2(end)).id_redcap,subj_info_struct(i_T2(end)).redcap_event_name,subj_info_struct(i_T2(end)).StudyInstanceUID,...
                                                  subj_info_struct(i_T2(end)).SeriesInstanceUID,subj_info_struct(i_T2(end)).site,'T2',parms,projinfo);
      end
      gap_T2=datenum(regexprep(num2str(subj_info_review(end).mrif_daic2rad_send),'_',''),'yyyymmdd') - datenum(num2str(subj_info_struct(i_T2(end)).StudyDate),'yyyymmdd');
      incoming2radreview_info(incoming2radreview_info_pos).mrif_daic2rad_gap=gap_T2;
      if ~isempty(incoming2radreview_info(incoming2radreview_info_pos).mrif_t2_unpackdate)
        gap_T2=datenum(incoming2radreview_info(incoming2radreview_info_pos).mrif_t2_unpackdate,parms.date_pat) - datenum(num2str(subj_info_struct(i_T2(end)).StudyDate),'yyyymmdd');
      else
        gap_T2=[];
      end
      incoming2radreview_info(incoming2radreview_info_pos).mrif_unpack_gap=gap_T2;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function incoming2radreview_info = set_unpack_date(incoming2radreview_info,incoming2radreview_info_pos,fname_pc_json,pGUID,EventName,StudyInstanceUID,SeriesInstanceUID,SiteName,ScanType,parms,projinfo)
  if regexp(fname_pc_json,'/pc/')
    unpack_dir=regexprep(fname_pc_json,'/pc/.*','/unpack');
    unpackdir=dir(sprintf('%s/%s/%s_%s_*_%s_%s',unpack_dir,SiteName,pGUID,EventName,StudyInstanceUID,SeriesInstanceUID));
  else
    for i=1:length(projinfo)
      unpack_dir = projinfo(i).unpack;
      unpackdir=dir(sprintf('%s/%s/%s_%s_*_%s_%s',unpack_dir,SiteName,pGUID,EventName,StudyInstanceUID,SeriesInstanceUID));
      if ~isempty(unpackdir),break; end;
    end
  end
  if ~isempty(unpackdir)
    unpackinfo=dir(sprintf('%s/%s/%s/%s.unpacked',unpack_dir,SiteName,unpackdir(end).name,unpackdir(end).name));
    if isempty(unpackinfo)
      eval(sprintf('incoming2radreview_info(incoming2radreview_info_pos).mrif_%s_unpackdate = [];',lower(ScanType)));
    else
      eval(sprintf('incoming2radreview_info(incoming2radreview_info_pos).mrif_%s_unpackdate = datestr(unpackinfo.datenum,parms.date_pat);',lower(ScanType)));
    end
  else
    eval(sprintf('incoming2radreview_info(incoming2radreview_info_pos).mrif_%s_unpackdate = [];',lower(ScanType)));
  end
return

function [subj_info_struct,subj_info_review] = select_subj_series(pguidevent,series_info,review_info,parms)
  subj_info_struct = [];
  subj_info_review = [];
  pguidevents = {series_info.pguidevent};
  i_subj = find(strcmp(pguidevents,pguidevent));
  subj_info_struct = series_info(i_subj);

  idx = ~cellfun('isempty',strfind({review_info.pguidevent},pguidevent));
  subj_info_review  = review_info(idx);

  if isfield(subj_info_struct,'ABCD_Compliant')
    PC = {subj_info_struct.ABCD_Compliant};
  else
    PC = [];
  end
  if isfield(subj_info_struct,'Completed')
    CC = {subj_info_struct.Completed};
  else
    CC = [];
  end
  if isfield(subj_info_struct,'QC')
    QC = {subj_info_struct.QC};
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
    fprintf('%s: WARNING: no checked series for event %s\n',mfilename,pguidevent);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: no checked series for event %s\n',mfilename,pguidevent);
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
          fprintf('%s: WARNING: sending data with QC = 0 for event %s\n',mfilename,pguidevent);
        end;
      end;
    end;
  end;
  if isempty(i_good)
    fprintf('%s: WARNING: no good series for event %s\n',mfilename,pguidevent);
    if parms.flog>0
      fprintf(parms.flog,'%s: WARNING: no good series for event %s\n',mfilename,pguidevent);
    end;
    return;
  end;
  subj_info_struct = subj_info_struct(i_checked(i_good));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function series_info = set_pguidevents(series_info)
  if isfield(series_info,'pGUID')
    pguidevents = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                    {series_info.pGUID},{series_info.EventName},...
                    'UniformOutput',false);
    [series_info.pguidevent] = deal(pguidevents{:});
  else
    pguidevents = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                    {series_info.id_redcap},{series_info.redcap_event_name},...
                    'UniformOutput',false);
    [series_info.pguidevent] = deal(pguidevents{:});
  end
return;
