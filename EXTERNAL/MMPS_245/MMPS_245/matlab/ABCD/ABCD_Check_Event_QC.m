function abcd_check_event_qc(varargin)
%function abcd_check_event_qc(varargin)
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
%     {default = 'event_qc_info'}
%   'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  03/07/17 by Feng Xue
% Last Mod: 03/09/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin);

if ~exist(parms.fname_out,'file') || parms.forceflag
  % load input info file
  series_info = abcd_load_csv(parms.fname_in);
  % filter series by protocol compliance and completeness
  series_info = filter_series(series_info);
  % check qc for each event
  event_info = check_events(series_info);
  % write output
  mmil_struct2csv(event_info,parms.fname_out);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','/home/mmilrec14/MetaData/DAL_ABCD_QC',[],...
    'outdir','/home/mmilrec14/MetaData/DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'infix','merged_pcqcinfo',[],...
    'outfix','event_qc_info',[],...
    'forceflag',true,[false true],...
  });
  parms.fname_in = sprintf('%s/%s_%s.csv',...
    parms.indir,parms.instem,parms.infix);
  parms.fname_out = sprintf('%s/%s_%s.csv',...
    parms.outdir,parms.outstem,parms.outfix);
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_info = filter_series(series_info)
  % find series with required fields
  fprintf('%s: identifying valid series...\n',mfilename);
  ind_valid = find(~cellfun(@isempty,{series_info.ABCD_Compliant}));
  series_info = series_info(ind_valid);

  % find compliant, complete series, with orig and raw
  fprintf('%s: identifying compliant, complete series...\n',mfilename);
  ind_valid = find(strcmpi({series_info.ABCD_Compliant},'Yes') &...
                [series_info.Completed]==1);
  series_info = series_info(ind_valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function event_info = check_events(series_info)
  fprintf('%s: checking QC for each event...\n',mfilename);
  event_info = [];
  % get combined pGUID-EventNames for each series
  pguidevents = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                         {series_info.pGUID},{series_info.EventName},...
                         'UniformOutput',false);
  [tmp,ind_uniq] = unique(pguidevents,'first');
  for i=1:length(ind_uniq)
    j = ind_uniq(i);
    % copy event-specific info
    event_info(i).pguidevent = pguidevents{j};
    event_info(i).id_redcap = series_info(j).pGUID;
    event_info(i).redcap_event_name = series_info(j).EventName;
    event_info(i).SiteName = series_info(j).SiteName;
    event_info(i).StudyInstanceUID = series_info(j).StudyInstanceUID;
    event_info(i).StudyDate = series_info(j).StudyDate;
    % check QC
    ind_series = find(strcmp(event_info(i).pguidevent,pguidevents));
    qc_vals = {series_info(ind_series).QC};
    if any(cellfun(@isempty,qc_vals))
      event_info(i).iqc_event_qcd = 0;
    else
      event_info(i).iqc_event_qcd = 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

