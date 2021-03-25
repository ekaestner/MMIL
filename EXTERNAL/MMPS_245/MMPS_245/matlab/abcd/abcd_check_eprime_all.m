function fname_out = abcd_check_eprime_all(varargin)
%function fname_out = abcd_check_eprime_all([options])
%
% Purpose: identify eprime files for each task fMRI series
%
% Optional Input: ('key', value pairs)
%   'indir': root input directory containing eprime files
%     {default = '/space/syn07/1/data/ABCD/aux_incoming'}
%   'fname_info': spreadsheet containing PC info for each series
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_pcinfo.csv'}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'ABCD_eprime'}
%   'scanner_flag': [0|1] find eprime files only for tasks done in scanner
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
%   'sers_eprime_time_diff_max': maximum time differene allowed between image-series and eprime files
%     {default = NaN}
% 
% Created:  11/21/16 by Jose Teruel} )
% Prev Mod: 08/03/17 by Don Hagler
% Last Mod: 08/14/17 by Feng Xue  and  Octavio Ruiz, 2017aug04-oct02
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: based on abcd_pc_eprime_all, created by Jose Teruel 11/21/16
%       last mod 11/28/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin);

fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,parms.outfix);
if ~exist(fname_out,'file') || parms.forceflag

  fprintf('%s: reading %s...\n',mfilename,parms.fname_info);

  pc_info = abcd_load_csv(parms.fname_info);
  pc_info = set_pguidevents(pc_info);
  % exclude non-task-fMRI series
  pc_info_fmri = select_task_fMRI(pc_info);

  if ~parms.scanner_flag
    filesys_info = get_filesys_info(parms);
    filesys_info = set_pguidevents(filesys_info);
  end;

  eprime_info = [];

  if parms.scanner_flag
    % for pguidevents exist in pc_info_fmri
    pguidevent_uniq = unique({pc_info_fmri.pguidevent});
    for i=1:length(pguidevent_uniq)
      output_info = check_pguidevent(pguidevent_uniq{i},pc_info_fmri,[],parms,1);
      eprime_info = cat(1,eprime_info,output_info);
    end;
  else
    % for pguidevents exist in pc_info_fmri
    idx = ismember({filesys_info.pguidevent},{pc_info_fmri.pguidevent});
    filesys_info_fmri = filesys_info(idx);
    for i=1:length(filesys_info_fmri)
      output_info = check_pguidevent(filesys_info_fmri(i).pguidevent,pc_info_fmri,filesys_info_fmri,parms,1);
      eprime_info = cat(1,eprime_info,output_info);
    end;

    % for pguidevents exist in filesys_info only
    filesys_info_rest = filesys_info(~idx);
    for i=1:length(filesys_info_rest)
      output_info = check_pguidevent(filesys_info_rest(i).pguidevent,pc_info_fmri,filesys_info_rest,parms,2);
      eprime_info = cat(1,eprime_info,output_info);
    end;

    % in case there are events that exist in pc_info_fmri but not in filesys_info
    idx = ismember({pc_info_fmri.pguidevent},{filesys_info.pguidevent});
    pc_info_rest = pc_info_fmri(~idx);
    pguidevent_uniq = unique({pc_info_rest.pguidevent});
    for i=1:length(pguidevent_uniq)
      output_info = check_pguidevent(pguidevent_uniq{i},pc_info_rest,[],parms,3);
      eprime_info = cat(1,eprime_info,output_info);
    end;
  end;

  eprime_info = mmil_sortstruct(eprime_info,{'SiteName','pGUID','EventName'}); 
  % save info to csv file
  mmil_struct2csv(eprime_info,fname_out);
  fname_cache = sprintf('%s/cache/%s_%s_%s.mat',...
    parms.outdir,parms.outstem,parms.outfix,datestr(now,'yyyymmdd'));
  abcd_info = eprime_info;
  save(fname_cache,'abcd_info');

end;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','/space/syn07/1/data/ABCD/aux_incoming',[],...
    'fname_info','/home/abcddaic/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_pcinfo.csv',[],...
    'tasknames',{'MID','SST','nBack'},{'MID','SST','nBack'},...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'outfix','eprime',[],...
    'scanner_flag',false,[false true],...
    'forceflag',false,[false true],...
    'sers_eprime_time_diff_max',NaN,[], ...   % Octavio Ruiz (2017oct02)
  ...
    'pcinfo_col_names',{'pGUID','VisitID','EventName','SessionType',...
                        'StudyDate','SeriesTime','SiteName','SeriesType',...
                        'ABCD_Compliant','SeriesDescription','Completed',...
                        'AdditionalInfo','SeriesInstanceUID','StudyInstanceUID_SeriesTime'},[],...
  });
  if ~exist(parms.indir)
    error('input directory %s not found',parms.indir);
  end;
  if ~exist(parms.fname_info,'file')
    error('info file %s not found',parms.fname_info);
  end;
  if parms.outdir(1) ~= '/'
     parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir);
  end
  if ~iscell(parms.tasknames)
    parms.tasknames = {parms.tasknames};
  end;
  parms.series_types = cell(size(parms.tasknames));
  for i=1:length(parms.tasknames)
    parms.series_types{i} = sprintf('fMRI_%s_task',parms.tasknames{i});
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output_info = select_task_fMRI(pc_info,parms)
  k = regexp({pc_info.SeriesType},'fMRI_\w+_task');
  ind = find(~cellfun(@isempty,k));
  output_info = pc_info(ind);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% Octavio (2017jul31-aug04, sep21): --------------------

function sitename = check_sitename(sitename)
  sitename = upper(sitename);
  switch sitename
    case 'DAIC'
      sitename = 'UCSD';
    case 'WASHU'
      sitename = 'WUSTL';
    case 'OAHU'
      sitename = 'UMB';
  end;
return;

function  tdiff = eprime_t_rel_to_series( StudyDate, SeriesTime, eprime_date, eprime_time )
    % Calculate time difference between: eprime_file date&time and Imag.series date&time
    %
    % Date Time are not strings but numbers that represent literally date and time.
    % The need to be translated to strings and then to datetime info.
    %
    tdiff = NaN;
    % Get integer part of SeriesTime, discarding fractions of a second
    ser_datt_str = datetime([num2str(StudyDate) ',' num2str(floor(SeriesTime))], ...
                            'InputFormat', 'yyyyMMdd,HHmmss');
    ser_datt_n  =  datenum( ser_datt_str );
    
    epr_datt_str = datetime([num2str(eprime_date) ',' num2str(floor(eprime_time))], ...
                            'InputFormat', 'yyyyMMdd,HHmmss');
    epr_datt_n  =  datenum( epr_datt_str );

    tdiff = epr_datt_n - ser_datt_n;   % days
    tdiff = tdiff * 24 * 60;           % minutes
return;
% ------------ :Octavio (2017jul31-aug04) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% Octavio (2017aug04-14): ---------------------
% % function output_info = check_series(pc_info,parms)

function output_info = check_series(pc_info,parms, mode)

% ---------------- :Octavio (2017aug04-14) %%%%%%%%%%%%%%%%%%%%%%
  output_info = [];
  k = regexp(pc_info.SeriesType,'fMRI_(?<taskname>\w+)_task','names');
  if isempty(k), return; end;

  eprime_info = abcd_check_eprime(pc_info.pGUID,...
                                  pc_info.EventName,pc_info.SessionType,...
                                  k.taskname,pc_info.SiteName, parms.indir);
  if ~isempty(eprime_info)
    tags = parms.pcinfo_col_names;
    for f=1:length(tags)
      output_info.(tags{f}) = pc_info.(tags{f});
    end;
    tags = fieldnames(eprime_info);
    for f=1:length(tags)
      output_info.(tags{f}) = eprime_info.(tags{f});
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Octavio: -----------------------------
    % Calculate time difference between: eprime-file date&time and Imag.series date&time
    % Append result to output record
    output_info.eprime_t_series_mode   = mode;
    output_info.eprime_t_rel_to_series = NaN;
    output_info.eprime_t_series_match  =  0;

    if output_info.eprime_datetime_found
        try
            tdiff = eprime_t_rel_to_series(pc_info.StudyDate, pc_info.SeriesTime, ...
                                           eprime_info.eprime_date, eprime_info.eprime_time);
            output_info.eprime_t_rel_to_series = round(tdiff,2);
            
            if abs(tdiff) <= parms.sers_eprime_time_diff_max
                output_info.eprime_t_series_match = 1;
            end
        catch err
            fprintf('Unable to calculate output_info.eprime_t_rel_to_series\n');
            fprintf('%s\n', err.message);
            fprintf('parms.fname_info, eprime_info.eprime_file_name:\n');
            fprintf('%s\n', parms.fname_info, eprime_info.eprime_file_name);
            fprintf('pc_info.StudyDate & SeriesTime = %.0f, %.0f;   ', pc_info.StudyDate, pc_info.SeriesTime);
            fprintf('eprime_info.eprime_date & .eprime_time = %.0f, %.0f\n\n', ...
                     eprime_info.eprime_date, eprime_info.eprime_time);
        end
    end
    % ---------- :Octavio (2017jul31-sep07) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%% Octavio (2017aug04-14): ---------------------

% % function output_info = check_series_file(pGUID,EventName,taskname,SiteName,parms)
function output_info = check_series_file(pGUID,EventName,taskname,SiteName, parms, StudyDate,SeriesTime, mode)

% ---------------- :Octavio (2017aug04-14) %%%%%%%%%%%%%%%%%%%%%%

  output_info = [];

% Octavio (2017sep21): ---------------------------------------------------
% % %   eprime_info = abcd_check_eprime(pGUID,EventName,[],taskname, SiteName, parms.indir);
  eprime_info = abcd_check_eprime(pGUID,EventName,[],taskname, check_sitename(SiteName), parms.indir);
%---------------------------------------------------- :Octavio (2017sep21)    

  if ~isempty(eprime_info)
    output_info.pGUID = pGUID;
    output_info.VisitID = '';
    output_info.EventName = EventName;
    output_info.SessionType = '';
    output_info.StudyDate = '';
    output_info.SeriesTime = '';

% Octavio (2017sep21): ---------------------------------------------------
% % %     output_info.SiteName = SiteName;
    output_info.SiteName = check_sitename(SiteName);
%---------------------------------------------------- :Octavio (2017sep21)    

    output_info.SeriesType = sprintf('fMRI_%s_task',taskname);
    output_info.ABCD_Compliant = '';
    output_info.SeriesDescription = '';
    output_info.Completed = '';
    output_info.AdditionalInfo = '';
    output_info.SeriesInstanceUID = '';
    output_info.StudyInstanceUID_SeriesTime = '';

    tags = fieldnames(eprime_info);
    for f=1:length(tags)
      output_info.(tags{f}) = eprime_info.(tags{f});
    end;
    
    %%%%%%%%%%%%%%%%%% Octavio: ----------------------------------------
    output_info.eprime_t_series_mode   = mode;
    output_info.eprime_t_rel_to_series = NaN;
    output_info.eprime_t_series_match  =  0;
    
    if output_info.eprime_datetime_found
        try
            tdiff = eprime_t_rel_to_series(StudyDate, SeriesTime, ...
                                           eprime_info.eprime_date, eprime_info.eprime_time);
            output_info.eprime_t_rel_to_series = round(tdiff,2);

            if abs(tdiff) <= parms.sers_eprime_time_diff_max
                output_info.eprime_t_series_match = 1;
            else
                output_info.eprime_t_series_match = 0;
            end
        catch err
            fprintf('Unable to calculate output_info.eprime_t_rel_to_series\n');
            fprintf('%s\n', err.message);
            fprintf('parms.fname_info, eprime_info.eprime_file_name:\n');
            fprintf('%s\n', parms.fname_info, eprime_info.eprime_file_name);
            fprintf('StudyDate = %.0f, SeriesTime = %.0f', StudyDate, SeriesTime);
            fprintf(';  eprime_info.eprime_date & .eprime_time = %.0f, %.0f\n\n', ...
                        eprime_info.eprime_date, eprime_info.eprime_time);       
        end
    end
    % ------------------- :Octavio (2017jul25-sep07) %%%%%%%%%%%%%%%%%%
    
  end;
return;


function output_info = check_pguidevent(pguidevent,pc_info,filesys_info,parms,mode)
%mode
%  1. matched fmri pguidevents (exist in both filesys_info and pc_info_fmri)
%  2. non-matched pguidevents from filesys_info (exist only in filesys_info)
%  3. non-matched pguidevents from pc_info (exist only in pc_info_fmri)
  output_info = [];
  [pGUID, EventName] = decode_pguidevent(pguidevent);
  switch mode
    case 1
      idx = ~cellfun('isempty',strfind({pc_info.pguidevent},pguidevent));
      pc_info = pc_info(idx);
      for i=1:length(parms.tasknames)
        taskcount.(parms.tasknames{i}) = 0;
      end
      for i=1:length(pc_info)
        k = regexp(pc_info(i).SeriesType,'fMRI_(?<taskname>\w+)_task','names');
        if ~isempty(k)

          %%%%%%%%%%%%%%%%%%%% Octavio (2017aug03-14, sep21): ------------------
% %           output_info = [output_info; check_series(pc_info(i),parms)];
% %           output_info = [output_info; check_series(pc_info(i),parms, mode)];

          pc_info_crrctd = pc_info(i);
          pc_info_crrctd.SiteName = check_sitename(pc_info_crrctd.SiteName);
          output_info = [output_info; check_series(pc_info_crrctd, parms, mode)];
          %------------------ :Octavio (2017aug03-14, sep21) %%%%%%%%%%%%%%%%%%%

          taskcount.(k.taskname) = taskcount.(k.taskname) +1;
        end
      end
      for i=1:length(parms.tasknames)
        if taskcount.(parms.tasknames{i}) == 0

          %%%%%%%%%%%%%%%%%%%% Octavio (2017aug03-14): ------------------

% %           output_info = [output_info; check_series_file(pGUID,EventName,parms.tasknames{i},pc_info(1).SiteName,parms)];
          output_info = [output_info; ...
                         check_series_file(pGUID,EventName,parms.tasknames{i}, pc_info(1).SiteName, parms, pc_info(1).StudyDate, pc_info(1).SeriesTime,mode)];
          %------------------ :Octavio (2017aug03-14) %%%%%%%%%%%%%%%%%%%

        end
      end

    case 2
      idx = ~cellfun('isempty',strfind({filesys_info.pguidevent},pguidevent));
      filesys_info = filesys_info(idx);
      for i=1:length(parms.tasknames)

        %%%%%%%%%%%%%%%%%%%% Octavio (2017aug03-14): -------------------
% %         output_info = [output_info; check_series_file(pGUID,EventName,parms.tasknames{i},filesys_info(1).site,parms)];
        output_info = [output_info; ...
                       check_series_file(pGUID,EventName,parms.tasknames{i}, filesys_info(1).site, parms, [], [], mode)];
        %------------------ :Octavio (2017aug03-14) %%%%%%%%%%%%%%%%%%%%

      end

    case 3
      idx = ~cellfun('isempty',strfind({pc_info.pguidevent},pguidevent));
      pc_info = pc_info(idx);
      for i=1:length(pc_info)
        output_info_tmp = [];
        k = regexp(pc_info(i).SeriesType,'fMRI_(?<taskname>\w+)_task','names');
        if ~isempty(k)
          output_info_tmp.pGUID = pGUID;
          output_info_tmp.VisitID = pc_info(i).VisitID;
          output_info_tmp.EventName = EventName;
          output_info_tmp.SessionType = pc_info(i).SessionType;
          output_info_tmp.StudyDate = pc_info(i).StudyDate;
          output_info_tmp.SeriesTime = pc_info(i).SeriesTime;

% Octavio (2017sep21): ---------------------------------------------------
% % %           output_info_tmp.SiteName = pc_info(i).SiteName;
          output_info_tmp.SiteName = check_sitename(pc_info(i).SiteName);
%---------------------------------------------------- :Octavio (2017sep21)    

          output_info_tmp.SeriesType = pc_info(i).SeriesType;
          output_info_tmp.ABCD_Compliant = pc_info(i).ABCD_Compliant;
          output_info_tmp.SeriesDescription = pc_info(i).SeriesDescription;
          output_info_tmp.Completed = pc_info(i).Completed;
          output_info_tmp.AdditionalInfo = pc_info(i).AdditionalInfo;
          output_info_tmp.SeriesInstanceUID = pc_info(i).SeriesInstanceUID;
          output_info_tmp.StudyInstanceUID_SeriesTime = pc_info(i).StudyInstanceUID_SeriesTime;
          output_info_tmp.taskname = k.taskname;
          output_info_tmp.eprime_folder_found = 0;
          output_info_tmp.eprime_file_found = 0;
          output_info_tmp.eprime_pguid_match = 0;
          output_info_tmp.eprime_file_naming = 0;
          output_info_tmp.eprime_file_name = '';
          output_info_tmp.eprime_extra_folder_found = 0;
          output_info_tmp.eprime_extra_file_found = 0;
          output_info_tmp.eprime_extra_pguid_match = 0;
          output_info_tmp.eprime_extra_file_naming = 0;
          output_info_tmp.eprime_extra_file_name = '';

          %%%%%%%%%%%%%%%%%% Octavio: -----------------------------------
          output_info_tmp.eprime_datetime_found  =  0;
          output_info_tmp.eprime_file_res = 0;
          output_info_tmp.eprime_date = NaN;
          output_info_tmp.eprime_time = NaN;
          output_info_tmp.eprime_t_series_mode   = mode;
          output_info_tmp.eprime_t_rel_to_series = NaN;
          output_info_tmp.eprime_t_series_match  =  0;
          % --------------- :Octavio (2017jul25-sep07) %%%%%%%%%%%%%%%%%%

          output_info = [output_info; output_info_tmp];
        end
      end
    otherwise
      error('wrong mode');
  end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pGUID, EventName] = decode_pguidevent(pguidevent)
  k = regexp(pguidevent,'_','split');
  if length(k)>=2
    pGUID = sprintf('%s_%s',k{1},k{2});
  else
    pGUID = 'invalid';
  end;
  if length(k)>=5
    %EventName = sprintf('%s_%s_%s',k{3},k{4},k{5});
    for i = 3:length(k)
      switch i
        case 3
          EventName = k{i};
        otherwise
          EventName = sprintf('%s_%s',EventName,k{i});
      end
    end
    %EventName = strjoin(k(3:end),'_');
  else
    EventName = 'invalid';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filesys_info = get_filesys_info(parms)
  fprintf('%s: getting filesys info...\n',mfilename);
  tic;
  if 1
    filesys_info = [];
    k = 1;
    slist = dir(sprintf('%s/*',parms.indir));
    slist = setdiff({slist.name},{'.','..'});
    for s=1:length(slist)
      site = slist{s};
      plist = dir(sprintf('%s/%s/*',parms.indir,site));
      plist = setdiff({plist.name},{'.','..'});
      for p=1:length(plist)
        pGUID = plist{p};
        if regexp(pGUID, '^NDAR_INV[^_]*$')
          elist = dir(sprintf('%s/%s/%s/*',parms.indir,site,pGUID));
          elist = setdiff({elist.name},{'.','..'});
          for e=1:length(elist)
            EventName = elist{e};
            if isempty(strfind(EventName,'backup'))
              filesys_info(k).site = lower(site);
              filesys_info(k).pGUID = pGUID;
              filesys_info(k).EventName = EventName;
              k = k + 1;
            end;
          end;
        end;
      end;
    end;
  else
    filesys_info = [];
    filesys_info_pos = 0;
    if parms.indir(end) == '/', parms.indir=parms.indir(1:end-1); end;
    cmd = sprintf('echo `ls -d %s/*/NDAR_INV*/*`',parms.indir);
    [errcode,dirlist] = unix(cmd);
    if errcode
      error('cmd %s failed:\n%s\n',cmd,dirlist);
    end;
%    [~,dirlist] = unix(sprintf('ls -d %s/*/NDAR_INV*/*',parms.indir));
    dirlist = regexprep(dirlist,'\s/',';/');
    dirlist = regexprep(dirlist,'\s','');
    dirlist = regexprep(dirlist,parms.indir,'');
    index = strfind(dirlist,';');
    for i=1:length(index)+1
      switch i
        case 1
          sitestr = dirlist(1:index(i)-1);
        case length(index)+1
          sitestr = dirlist(index(i-1)+1:end);
        otherwise
          sitestr = dirlist(index(i-1)+1:index(i)-1);
      end
      sitestr = regexprep(sitestr,'^/','');
      if regexp(sitestr,'^[^/]*/NDAR_INV[^_/]*/.*') & isempty(strfind(sitestr,'backup'))
        slashindex = strfind(sitestr,'/');
        filesys_info_pos = filesys_info_pos +1;
        if slashindex(1) ==1
          filesys_info(filesys_info_pos).site = lower(sitestr(2:slashindex(2)-1));
          filesys_info(filesys_info_pos).pGUID = sitestr(slashindex(2)+1:slashindex(3)-1);
          filesys_info(filesys_info_pos).EventName = sitestr(slashindex(3)+1:end);
        else
          filesys_info(filesys_info_pos).site = lower(sitestr(1:slashindex(1)-1));
          filesys_info(filesys_info_pos).pGUID = sitestr(slashindex(1)+1:slashindex(2)-1);
          filesys_info(filesys_info_pos).EventName = sitestr(slashindex(2)+1:end);
        end
      end
    end
  end;
  toc;
return
