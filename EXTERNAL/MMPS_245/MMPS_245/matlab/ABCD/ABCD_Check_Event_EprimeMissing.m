function abcd_check_event_eprimemissing(varargin)
%function abcd_check_event_eprimemissing(varargin)
%
% Optional input:
%  'outdir': output metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'outstem': output file stem
%     {default = 'DAL_ABCD_QC'}
%  'instem': input file stem
%     {default = 'DAL_ABCD_QC'}
%  'infix': file suffix of eprime csv
%       {default = 'combined_eprime'}
%  'ra_infix': file suffix of ra_checklist
%       {default = 'ra_checklist'}
%  'nbackrec_infix': file suffix of nback memo from redcap
%       {default = 'nback_rec'}
%  'indir': input metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'outfix': file suffix of output file
%     {default = 'combined_eprime_missing'}
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Required input:   Octavio Ruiz, 2017oct10
%  'sers_eprime_time_diff_max': Eprime_T_Series_Max (minutes)
% 
% Created:  03/30/17 by Feng Xue
% Last Mod: 05/25/17 by Feng Xue,  and  Octavio Ruiz, 2017aug18-25, oct10
%
% todo: make sure pGUID are from ra_checklist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check input parameters


% % % % % % Octavio: ----------
parms = check_input(varargin);
Eprime_T_Series_Max = parms.sers_eprime_time_diff_max;
% ---------- :Octavio % % % % %


%check lock files
fname_lck=sprintf('%s/.check_event_eprimemissing.lck',parms.outdir);
if exist(fname_lck,'file')
  fprintf('%s\n','lock files exist!.');
  return;
end

%Place lock file
fclose(fopen(fname_lck, 'w'));

fname_combined_info=sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,parms.outfix);

fname_eprime_info = sprintf('%s/%s_%s.csv',parms.indir,parms.instem,parms.infix);
fname_ra_checklist = sprintf('%s/%s_%s.csv',parms.indir,parms.instem,parms.ra_infix);
fname_nbackrec = sprintf('%s/%s_%s.csv',parms.indir,parms.instem,parms.nbackrec_infix);
if ~exist(fname_eprime_info,'file'), error('file %s not found',fname_eprime_info); end;
if ~exist(fname_ra_checklist,'file'), error('file %s not found',fname_ra_checklist); end;
if ~exist(fname_nbackrec,'file'), error('file %s not found',fname_nbackrec); end;

eprime_info = abcd_load_csv(fname_eprime_info);
ra_checklist_data = abcd_load_csv(fname_ra_checklist);

% data reduction
ra_checklist_data = ra_checklist_data(cellfun(@(x) isequal(x,1),{ra_checklist_data.enroll_total_1}));
nbackrec_data = abcd_load_csv(fname_nbackrec);
nbackrec_data = nbackrec_data(cellfun(@(x) isequal(x,1),{nbackrec_data.enroll_total_1}));

[eprime_info_tmp(1:numel(eprime_info)).pGUID] = deal(eprime_info.pGUID);
[eprime_info_tmp(1:numel(eprime_info)).EventName] = deal(eprime_info.EventName);
[eprime_info_tmp(1:numel(eprime_info)).site] = deal(eprime_info.SiteName);
[eprime_info_tmp(1:numel(eprime_info)).taskname] = deal(eprime_info.taskname);
[eprime_info_tmp(1:numel(eprime_info)).eprime_file_found] = deal(eprime_info.eprime_file_found);
[eprime_info_tmp(1:numel(eprime_info)).eprime_pguid_match] = deal(eprime_info.eprime_pguid_match);
[eprime_info_tmp(1:numel(eprime_info)).eprime_extra_file_found] = deal(eprime_info.eprime_extra_file_found);
[eprime_info_tmp(1:numel(eprime_info)).eprime_extra_pguid_match] = deal(eprime_info.eprime_extra_pguid_match);


% % % % % % Octavio: ----------

[eprime_info_tmp(1:numel(eprime_info)).eprime_datetime_found]  = deal(eprime_info.eprime_datetime_found);
[eprime_info_tmp(1:numel(eprime_info)).StudyDate]  = deal(eprime_info.StudyDate);
[eprime_info_tmp(1:numel(eprime_info)).SeriesTime] = deal(eprime_info.SeriesTime);
[eprime_info_tmp(1:numel(eprime_info)).eprime_file_res]        = deal(eprime_info.eprime_file_res);
[eprime_info_tmp(1:numel(eprime_info)).eprime_date]            = deal(eprime_info.eprime_date);
[eprime_info_tmp(1:numel(eprime_info)).eprime_time]            = deal(eprime_info.eprime_time);
[eprime_info_tmp(1:numel(eprime_info)).eprime_t_series_mode]   = deal(eprime_info.eprime_t_series_mode);
[eprime_info_tmp(1:numel(eprime_info)).eprime_t_rel_to_series] = deal(eprime_info.eprime_t_rel_to_series);
[eprime_info_tmp(1:numel(eprime_info)).eprime_t_series_match]  = deal(eprime_info.eprime_t_series_match);

% ---------- :Octavio % % % % %


eprime_info=eprime_info_tmp';
clear eprime_info_tmp;

combined_info = eprime_stats( eprime_info, ra_checklist_data, nbackrec_data, Eprime_T_Series_Max);

% write output file
mmil_struct2csv(combined_info,fname_combined_info);

%delete lock file
delete(fname_lck);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function combined_info = eprime_stats( eprime_info, ra_info, nbackrec_info, Eprime_T_Series_Max)
  eprime_info = set_pguidevents(eprime_info);
  ra_info = set_pguidevents(ra_info);
  nbackrec_info = set_pguidevents(nbackrec_info);
  unique_pguidevent=unique({eprime_info.pguidevent});
  unique_pguidevent_ra=unique({ra_info.pguidevent});
  idx=ismember(unique_pguidevent,unique_pguidevent_ra);
  unique_pguidevent=unique_pguidevent(idx);
  %combined_info=cell2struct(unique_pguidevent,'pguidevent');
  combined_info=[];
  tasklist = {'sst','mid','nback'};
  tasklist_ra = {'sstrc','rcom','vemorc'};
  tasklist_lap = {'sst','mid','nbac'};

  for i=1:length(unique_pguidevent)
    idx = ~cellfun('isempty',strfind({eprime_info.pguidevent},unique_pguidevent{i}));
    event_data = eprime_info(idx);
    idx = ~cellfun('isempty',strfind({ra_info.pguidevent},unique_pguidevent{i}));
    ra_data = ra_info(idx);
    idx = ~cellfun('isempty',strfind({nbackrec_info.pguidevent},unique_pguidevent{i}));
    nbackrec_data = nbackrec_info(idx);

    combined_info(i).pguidevent=unique_pguidevent{i};
    combined_info(i).id_redcap=event_data(1).pGUID;
    combined_info(i).redcap_event_name=event_data(1).EventName;
    combined_info(i).site=event_data(1).site;
    combined_info(i).iqc_eprime_lap_missing=0;
    combined_info(i).iqc_eprime_scan_missing=0;


    % % % % % % Octavio: ----------

    % Set per-sesion variables (as opposed to per-task vars.)
    combined_info(i).iqc_eprime_scan_t_match = [];
    
    scan_and_eprime_num = 0;   % Number of tasks for wich scan and eprime files were found and readable
    t_match_num = 0;           % Number of tasks with both scan and eprime files, and where t_diff < Threshold

    % ---------- :Octavio % % % % %


    for t=1:length(tasklist)
      idx = ~cellfun('isempty',strfind(lower({event_data.taskname}),tasklist{t}));
      if isempty(find(idx==1))
        record_cnt = 0;
      else
        record_cnt = 1;
        taskdata = event_data(idx);
        idx = cellfun(@(x) isequal(x,0), {taskdata.eprime_file_found});
        if ~isempty(find(idx==1)), record_cnt = 0; end;
      end

      if record_cnt < 1
        eval(sprintf('combined_info(i).iqc_%s_ep_file_found=0;',tasklist{t}));
        eval(sprintf('isCollected = ra_data.ra_scan_check_list_%s;',tasklist_ra{t}));
        if ~isempty(isCollected) && isCollected >0
          eval(sprintf('collection_method = ra_data.ra_scan_cl_%s_scan_lap;',tasklist_lap{t}));
          if collection_method == 2
            eval(sprintf('combined_info(i).iqc_%s_ep_lap_missing=1;',tasklist{t}));
            eval(sprintf('combined_info(i).iqc_%s_ep_scan_missing=0;',tasklist{t}));
            combined_info(i).iqc_eprime_lap_missing=1;
          else
            eval(sprintf('combined_info(i).iqc_%s_ep_lap_missing=0;',tasklist{t}));
            eval(sprintf('combined_info(i).iqc_%s_ep_scan_missing=1;',tasklist{t}));
            combined_info(i).iqc_eprime_scan_missing=1;
          end
        else
          eval(sprintf('combined_info(i).iqc_%s_ep_lap_missing=0;',tasklist{t}));
          eval(sprintf('combined_info(i).iqc_%s_ep_scan_missing=0;',tasklist{t}));
        end

      else
        eval(sprintf('combined_info(i).iqc_%s_ep_file_found=1;',tasklist{t}));
        idx = cellfun(@(x) isequal(x,0), {taskdata.eprime_pguid_match});
        if ~isempty(find(idx==1)), eval(sprintf('combined_info(i).iqc_%s_ep_pguid_mis=1;',tasklist{t})); end;
        
        eval(sprintf('combined_info(i).iqc_%s_ep_lap_missing=0;',tasklist{t}));
        eval(sprintf('combined_info(i).iqc_%s_ep_scan_missing=0;',tasklist{t}));
        

        % % % % % % Octavio: ----------
        [tmin,tmax,event] = find_time_diffs( event_data, tasklist{t} );
        
        tag = sprintf('iqc_%s_ep_datetime_found', tasklist{t});        combined_info(i).(tag) = event.eprime_datetime_found;
        tag = sprintf('iqc_%s_ep_file_res', tasklist{t});              combined_info(i).(tag) = event.eprime_file_res;
        tag = sprintf('iqc_%s_study_date',  tasklist{t});              combined_info(i).(tag) = dateconvert(event.StudyDate);
        tag = sprintf('iqc_%s_series_time', tasklist{t});              combined_info(i).(tag) = timeconvert(event.SeriesTime);

        tag = sprintf('iqc_%s_ep_date', tasklist{t});                  combined_info(i).(tag) = dateconvert(event.eprime_date);
        tag = sprintf('iqc_%s_ep_time', tasklist{t});                  combined_info(i).(tag) = timeconvert(event.eprime_time);
        tag = sprintf('iqc_%s_ep_t_rel_to_series_min', tasklist{t});   combined_info(i).(tag) = round( tmin, 2);
        tag = sprintf('iqc_%s_ep_t_rel_to_series_max', tasklist{t});   combined_info(i).(tag) = round( tmax, 2);

        tag = sprintf('iqc_%s_ep_t_series_match', tasklist{t});
        if event.eprime_datetime_found
            scan_and_eprime_num = scan_and_eprime_num + 1;
            if abs(tmin) <= Eprime_T_Series_Max
                combined_info(i).(tag) = true;
                t_match_num = t_match_num +1;
            else
                combined_info(i).(tag) = false;
            end
        else
% % %             combined_info(i).(tag) = [];
            combined_info(i).(tag) = -1;
        end
        % ---------- :Octavio % % % % %
        
      end;
      
      if strcmp(tasklist{t},'nback')
        combined_info(i).iqc_nback_epx_missing=0;
        idx = cellfun(@(x)isequal(x,1),{event_data.eprime_extra_file_found});
        if isempty(find(idx==1))
          record_cnt = 0;
        else
          record_cnt = 1;
          taskdata = event_data(idx);
          idx = cellfun(@(x)isequal(x,0),{taskdata.eprime_extra_pguid_match});
          if ~isempty(find(idx==1)), eval(sprintf('combined_info(i).iqc_%s_epx_pguid_mis=1;',tasklist{t})); end;
        end
            
        eval(sprintf('combined_info(i).iqc_%s_epx_file_found=%d;',tasklist{t},record_cnt));
        eval(sprintf('combined_info(i).iqc_%s_epx_missing=0;',tasklist{t}));
        if record_cnt == 0 && ~isempty(nbackrec_data)  && ~isempty(nbackrec_data.rec_mem) && nbackrec_data.rec_mem >0
          eval(sprintf('combined_info(i).iqc_%s_epx_missing=1;',tasklist{t}));
          %eval(sprintf('combined_info(i).iqc_%s_eprime_missing=1;',tasklist{t}));
          %combined_info(i).iqc_eprime_missing=1;
          combined_info(i).iqc_eprime_lap_missing=1;
        end;
      end
    end


    % % % % % % Octavio: ----------
    
% %     % If all task eprime files match an imaging series...
% %     if scan_and_eprime_num > 0
% %         if t_match_num >= scan_and_eprime_num
% %             combined_info(i).iqc_eprime_scan_t_match = true;
% %         else
% %             combined_info(i).iqc_eprime_scan_t_match = false;
% %         end
% %     else
% %         combined_info(i).iqc_eprime_scan_t_match = -1;
% %     end

    % If, for all tasks, eprime files were found and time extracted:
    if 0 < scan_and_eprime_num  &&  scan_and_eprime_num == length(tasklist)
        % If all time_differnces less than threshold:
        if t_match_num == scan_and_eprime_num
            combined_info(i).iqc_eprime_scan_t_match = true;
        else
            combined_info(i).iqc_eprime_scan_t_match = false;
        end
    else
        combined_info(i).iqc_eprime_scan_t_match = -1;
    end

    % ---------- :Octavio % % % % %


  end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % % % Octavio: ----------

function y = dateconvert(x)
% Convert single-number encoded date to string
    try
        y = datetime( num2str(x), 'InputFormat', 'yyyyMMdd' );
        y = datestr( y, 'yyyy-mm-dd' );
    catch
        y = '';
    end
return

function y = timeconvert(x)
% Convert single-number encoded date to string
    try
        y = datetime( num2str(floor(x)), 'InputFormat', 'HHmmss' );
        y = datestr( y, 'HH:MM' );
    catch
        y = '';
    end
return

function [tmin,tmax,event] = find_time_diffs( event_data, task )
    aaa  =  regexp( lower({event_data.taskname}), lower(task) );
    inds = find( ~cellfun(@isempty, aaa) );
    switch length(inds)
        case 0
            tmin  = NaN;
            tmax  = NaN;
            fprintf('\n+ abcd_check_event_eprimemissing: task %s absent in eprime_info +\n', task);
            % Instead of crashing, I should return here a structure with same, but empy-fields 
        case 1
            tmin  = event_data(inds).eprime_t_rel_to_series;
            tmax  = tmin;
            event = event_data(inds);
        otherwise
            tt = [event_data(inds).eprime_t_rel_to_series];
            [~,imin] = min(abs(tt));
            [~,imax] = max(abs(tt));
            tmin  = tt(imin);
            tmax  = tt(imax);
            event = event_data( inds(imin) );
    end
return

% ---------- :Octavio % % % % %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'outfix','eprime_missing',[],...
    'ra_infix','ra_checklist',[],...
    'nbackrec_infix','nback_rec',[],...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'infix','eprime',[],...
    'forceflag',true,[false true], ...
    'sers_eprime_time_diff_max',[],[], ...   % Octavio Ruiz (2017oct09)
  });
  if parms.outdir(1) ~= '/', parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir); end
  if parms.indir(1) ~= '/', parms.indir = parms.outdir; end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



% % % No,no,no:
% % % ddd = {eprime_info.StudyDate};
% % % [eprime_info_tmp(1:numel(eprime_info)).StudyDate]  = deal(cellfun(@dateconvert, ddd, 'UniformOutput', false));
% % % ttt = {eprime_info.SeriesTime};
% % % [eprime_info_tmp(1:numel(eprime_info)).SeriesTime] = deal(cellfun(@timeconvert, ttt, 'UniformOutput', false));
