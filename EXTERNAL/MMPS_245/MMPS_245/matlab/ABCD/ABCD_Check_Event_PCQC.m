function abcd_check_event_pcqc(varargin)
%function abcd_check_event_pcqc(varargin)
%
% Optional input:
%   'indir': incoming directory
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC'}
%   'outdir': output directory
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC'}
%   'instem': instem
%     {default = 'DAL_ABCD_QC'}
%   'outstem': outstem
%     {default = 'DAL_ABCD_QC'}
%   'infix : filename of input info file
%     {default = 'merged_pcqcinfo'}
%   'outfix':  file suffix for output file
%     {default = 'event_pcqc_info'}
%   'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  03/09/17 by Don Hagler
% Last Mod: 03/12/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin);
while exist(sprintf('%s/.merged_pcqcinfo.lck',parms.outdir),'file')
  fprintf('%s\n','lock files exist, waiting for previous process to finish.');
  pause(30);
end;


if ~exist(parms.fname_out,'file') || parms.forceflag
  % load series info
  series_info = abcd_load_csv(parms.fname_in);
  % select ABCD-relevant series
  series_info = get_relevant_series(series_info);
  % set pGUID_EventName for each series, get list of unique events
  [series_info,pguidevents] = set_pguidevents(series_info);
  % find qc notes fields
  parms = find_qc_notes(series_info,parms);
  % check each event
  nevents = length(pguidevents);
  fprintf('%s: checking %d events...\n',mfilename,nevents);
  tic;
  event_info = [];
  series_pguidevents = {series_info.pguidevent};
  for i=1:nevents
%    fprintf('%s: checking event %d of %d...\n',mfilename,i,nevents);
    ind_series = find(strcmp(pguidevents{i},series_pguidevents));
    new_info = check_event(series_info(ind_series),parms);
    event_info = concat_info(event_info,new_info);
  end;    
  toc;
  % write output
  fprintf('%s: writing results to %s...\n',mfilename,nevents,parms.fname_out);
  tic;
  event_info = mmil_sortstruct(event_info,{'SiteName','id_redcap','redcap_event_name'});
  mmil_struct2csv(event_info,parms.fname_out);
  toc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'infix','merged_pcqcinfo',[],...
    'outfix','event_pcqcinfo',[],...
    'forceflag',true,[false true],...
    ...
    'stypes',{'t1','t2','dmri','rsfmri','mid','sst','nback'},[],...
    'min_nscans',[1,1,1,3,2,2,2],[],...
    'max_nscans',[3,3,3,12,6,6,6],[],...
    'SeriesTypes',{'T1','T2','dMRI','rsfMRI','fMRI_MID_task','fMRI_SST_task','fMRI_nBack_task'},[],...
    'fmap_stypes',{'dmri','rsfmri','mid','sst','nback'},[],...
    'qc_issues',{'dis','dco','fa','gh','ht','hb','mo','rf','sd',...
                 'si','sus','vco','wr','SLICE','line','dark','moire',...
                 'RECON','MISS','OTHER'},[],...
  });
  if parms.outdir(1) ~= '/', parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir); end;
  if parms.indir(1) ~= '/', parms.indir= parms.outdir; end;
  parms.fname_in = sprintf('%s/%s_%s.csv',...
    parms.indir,parms.instem,parms.infix);
  parms.fname_out = sprintf('%s/%s_%s.csv',...
    parms.outdir,parms.outstem,parms.outfix);
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  parms.ntypes = length(parms.stypes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_info = get_relevant_series(series_info)
  % find series with required fields
  ind_valid = find(~cellfun(@isempty,{series_info.ABCD_Compliant}));
  series_info = series_info(ind_valid);
  % select series that are non-compliant or compliant series (exclude NA)
  %   and exclude Undefined or Undefined_fMRI series types
  ind_valid = find(ismember({series_info.ABCD_Compliant},{'Yes','No'}) &...
                   cellfun(@isempty,regexp({series_info.SeriesType},'Undefined')));
  series_info = series_info(ind_valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [series_info,pguidevents] = set_pguidevents(series_info)
  pguidevents = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                  {series_info.pGUID},{series_info.EventName},...
                  'UniformOutput',false);
  [series_info.pguidevent] = deal(pguidevents{:});  
  pguidevents = unique(pguidevents);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = find_qc_notes(series_info,parms)
  fnames = fieldnames(series_info);
  ind = find(~cellfun(@isempty,regexp(fnames,'notes_\w+')));
  parms.tags_qc_notes = fnames(ind);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function event_info = check_event(series_info,parms)
  % set event-specific variables
  event_info = init_event_info(series_info);
  % check protocol compliance and quality control for each series type
  event_info = check_series_pcqc(event_info,series_info,parms);
  % check whether event includes all series types
  event_info = check_event_completeness(event_info,parms);
  % check whether qc was performed for all valid series of this event
  event_info = check_event_qc(event_info,series_info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function event_info = init_event_info(series_info)
  event_info.pguidevent = series_info(1).pguidevent;
  event_info.id_redcap = series_info(1).pGUID;
  event_info.redcap_event_name = series_info(1).EventName;
  event_info.SiteName = series_info(1).SiteName;
  dates = sort([series_info.StudyDate]);
  event_info.LatestDate = dates(end);
  %% todo: initialize all other variables to [] so order is fixed?
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function event_info = check_series_pcqc(event_info,series_info,parms)
  SeriesTypes = {series_info.SeriesType};
  for i=1:parms.ntypes
    stype = parms.stypes{i};
    % identify scans for current series type
    ind_series = find(strcmp(parms.SeriesTypes{i},SeriesTypes));
    % get pc/qc info for individual scans
    [event_info,scan_info] = ...
      check_series_scans(event_info,series_info(ind_series),stype,parms);
    % set summary variables for event
    event_info = check_series_event(event_info,scan_info,stype,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,scan_info] = check_series_scans(event_info,...
                                                       series_info,stype,parms)
  % set index for stype
  i = find(strcmp(stype,parms.stypes));
  % initialize scan counters, etc.
  scan_info = init_scan_info(parms);
  % loop over possible scans for current series type
  for j=1:parms.max_nscans(i);
    if j>length(series_info), continue; end;
    % qc score
    qc_score = series_info(j).QC;
    tag = sprintf('iqc_%s_%d_qc_score',stype,j);
    event_info.(tag) = qc_score;
    if isempty(qc_score), qc_score = 0; end;
    % qc notes
    qc_notes = get_qc_notes(series_info(j),parms);
    tag = sprintf('iqc_%s_%d_qc_notes',stype,j);
    event_info.(tag) = qc_notes;
    % pc score
    if strcmpi(series_info(j).ABCD_Compliant,'Yes')
      pc_score = 1;
    else
      pc_score = 0;
    end;
    tag = sprintf('iqc_%s_%d_pc_score',stype,j);
    event_info.(tag) = pc_score;
    % pc notes
    pc_notes = series_info(j).AdditionalInfo;
    tag = sprintf('iqc_%s_%d_pc_notes',stype,j);
    event_info.(tag) = pc_notes;
    % completed (no missing dicoms)
    completed = series_info(j).Completed;
    tag = sprintf('iqc_%s_%d_complete',stype,j);
    event_info.(tag) = completed;
    % count number of scans
    scan_info.total_ser = scan_info.total_ser + 1;
    if pc_score
      scan_info.total_passpc = scan_info.total_passpc + 1;
    else
      scan_info.ser_pc_issues = scan_info.ser_pc_issues + 1;
    end;
    if pc_score && qc_score && completed
      scan_info.good_ser = scan_info.good_ser + 1;
    end;
    if ~qc_score
      scan_info.ser_qcs = scan_info.ser_qcs + 1;
    end;
    if ~completed
      scan_info.ser_incomp = scan_info.ser_incomp + 1;
    end;
    % check for QC issues (based on QC notes)
    scan_info = check_qc_issues(scan_info,qc_notes,parms);

    %% todo: check stype-specific variables for each scan (e.g fMRI motion sub_02)

    % check field map variables for each scan for some series types
    if ismember(stype,parms.fmap_stypes)

      %% todo: set variables about field maps

    end;
  end;
  % whether qc was done
  qc_vals = {series_info.QC};
  if any(~cellfun(@isempty,qc_vals))
    scan_info.qcd = 1;
  else
    scan_info.qcd = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scan_info = init_scan_info(parms)
  scan_info = [];
  scan_info.total_ser = 0;
  scan_info.total_passpc = 0;
  scan_info.good_ser = 0;
  scan_info.ser_pc_issues = 0;
  scan_info.ser_qcs = 0;
  scan_info.ser_incomp = 0;
  scan_info.qcd = 0;
  for i=1:length(parms.qc_issues)
    scan_info.(lower(parms.qc_issues{i})) = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scan_info = check_qc_issues(scan_info,qc_notes,parms)
  if isempty(qc_notes), return; end;
  for i=1:length(parms.qc_issues)
    qc_issue = parms.qc_issues{i};
    if ~isempty(regexp(qc_notes,qc_issue))
      scan_info.(lower(qc_issue)) = scan_info.(lower(qc_issue)) + 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function event_info = check_series_event(event_info,scan_info,stype,parms)
  % number of scans of current series type
  tag = sprintf('iqc_%s_total_ser',stype);
  event_info.(tag) = scan_info.total_ser;
  % number of protocol compliant scans
  tag = sprintf('iqc_%s_total_passpc',stype);
  event_info.(tag) = scan_info.total_passpc;
  % number of compliant, complete, and good QC scans
  if scan_info.qcd
    % only set if at least one series was qc'd
    tag = sprintf('iqc_%s_good_ser',stype);
    event_info.(tag) = scan_info.good_ser;
    tag = sprintf('iqc_%s_bad_ser',stype);
    event_info.(tag) = scan_info.total_ser - scan_info.good_ser;
    % set variables for QC issues
    for i=1:length(parms.qc_issues)
      qc_issue = lower(parms.qc_issues{i});
      tag = sprintf('iqc_%s_%s',stype,qc_issue);
      event_info.(tag) = scan_info.(qc_issue);
    end;
  else
    tag = sprintf('iqc_%s_good_ser',stype);
    event_info.(tag) = [];
    tag = sprintf('iqc_%s_bad_ser',stype);
    event_info.(tag) = [];
    for i=1:length(parms.qc_issues)
      tag = sprintf('iqc_%s_%s',stype,lower(parms.qc_issues{i}));
      event_info.(tag) = [];
    end;
  end;
  % number of scans that failed pc
  tag = sprintf('iqc_%s_ser_pc_issues',stype);
  event_info.(tag) = scan_info.ser_pc_issues;
  % number of scans that failed qc
  tag = sprintf('iqc_%s_ser_qcs',stype);
  event_info.(tag) = scan_info.ser_qcs;
  % number of incomplete scans (missing dicoms)
  tag = sprintf('iqc_%s_ser_incomp',stype);
  event_info.(tag) = scan_info.ser_incomp;

  %% todo: check stype-specific variables for each event (e.g fMRI motion)

  % check field map variables for each event for some series types
  if ismember(stype,parms.fmap_stypes)

    %% todo: set variables about field maps

  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function old_info = concat_info(old_info,new_info)
  if isempty(new_info), return; end;
  if isempty(old_info)
    old_info = new_info;
  else
    % reconcile fieldnames (add empty fields to match)
    [old_info,new_info] = reconcile_fields(old_info,new_info);
    old_info = cat(2,old_info,new_info);    
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,B] = reconcile_fields(A,B)
  fnamesA = fieldnames(A);
  fnamesB = fieldnames(B);
  fnames = setdiff(fnamesB,fnamesA);
  for f=1:length(fnames)
    A(1).(fnames{f}) = [];
  end;
  fnames = setdiff(fnamesA,fnamesB);
  for f=1:length(fnames)
    B(1).(fnames{f}) = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function notes = get_qc_notes(series_info,parms)
  notes = [];
  for i=1:length(parms.tags_qc_notes)
    note = series_info.(parms.tags_qc_notes{i});
    if ~isempty(note)
      notes = sprintf('%sR%d: %s ',notes,i,note);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function event_info = check_event_completeness(event_info,parms)
  % check that all series types have enough scans
  event_info.iqc_eventcomplete = 1;
  for i=1:parms.ntypes
    tag = sprintf('iqc_%s_total_ser',parms.stypes{i});
    if event_info.(tag) < parms.min_nscans(i)
      event_info.iqc_eventcomplete = 0;
      break;
    end;
  end;
  % check that all series types have enough scans with good pc/qc
  event_info.iqc_eventcomplete_c_pc_qc = 1;
  for i=1:parms.ntypes
    tag = sprintf('iqc_%s_good_ser',parms.stypes{i});
    if event_info.(tag) < parms.min_nscans(i)
      event_info.iqc_eventcomplete_c_pc_qc = 0;
      break;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function event_info = check_event_qc(event_info,series_info)
  series_info = get_valid_series(series_info);
  qc_vals = {series_info.QC};
  if any(cellfun(@isempty,qc_vals))
    event_info.iqc_event_qcd = 0;
  else
    event_info.iqc_event_qcd = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_info = get_valid_series(series_info)
  % find classified series
  ind_valid = find(strcmpi({series_info.ABCD_Compliant},'Yes') &...
                [series_info.Completed]==1);
  series_info = series_info(ind_valid);
  %% todo: exclude series for which processing failed (manual QC not possible)
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

