function eprime_info = abcd_check_eprime(pGUID,eventname,sesstype,taskname,sitename,indir)
%function eprime_info = abcd_check_eprime(pGUID,eventname,sesstype,taskname,sitename,[indir])
%
% Purpose: identify eprime file for given visit and task
%
% Required Input:
%   pGUID: globally unique ID
%   eventname: name of data acquisition "event"
%     e.g. 'baseline_year_1_arm_1'
%   sesstype: Session Type
%     i.e. 'SessionC'
%   taskname: name of fMRI task
%     i.e. 'MID', 'SST', or 'nBack'
%   sitename: name of data acquisition sitename
%     e.g. 'UCSD',' WUSTL', etc.
%
% Optional Input: ('key', value pairs)
%   indir: root input directory containing eprime files
%     {default = '/space/syn05/1/data/MMILDB/DAL_ABCD/aux_incoming'}
%
% Output:
%   eprime_info: struct containing fields:
%     'eprime_folder_found'
%     'eprime_file_found'
%     'eprime_file_naming'
%     'eprime_pguid_match'
%     'eprime_file_name'
%     'eprime_extra_file_found'
%     'eprime_extra_file_naming'
%     'eprime_extra_file_name'
%
%   % Octavio: -----------------------------------------------------------
%      eprime_datetime_found
%      eprime_file_res
%      eprime_date
%      eprime_time
%   %------------------------------------------ :Octavio (2017jul25-sep07)
% 
% Created:  12/20/16 by Don Hagler
% Prev Mod: 07/19/17 by Feng Xue
% Last Mod: 2017aug01-sep07 by Octavio Ruiz
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eprime_info = init_eprime_info(taskname);
if ~mmil_check_nargs(nargin,5), return; end;
if ~exist('indir','var') || isempty(indir)
  indir = '/space/syn05/1/data/MMILDB/DAL_ABCD/aux_incoming';
end;
if ~exist('extlist','var') || isempty(extlist)
  extlist = {'.csv','.txt'};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch taskname
  case 'nBack'
    dir_infix = 'nback-wm';
    task_infix = 'WM';
    extra_task_infix = 'REC';
    extra_dir_infix = 'nback-rec';
  case 'MID'
    dir_infix = 'mid';
    task_infix = 'MID';
    extra_task_infix = [];
    extra_dir_infix = [];
  case 'SST'
    dir_infix = 'sst';
    task_infix = 'SST';
    extra_task_infix = [];
    extra_dir_infix = [];
  otherwise
    error('invalid task name %s',taskname);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find event directory

% % % Octavio (2017sep21) removed:  sitename = check_sitename(sitename);

eventdir = fullfile(indir,sitename,pGUID,eventname);
if ~exist(eventdir,'dir'), return; end;

% find task directory
dpat = sprintf('%s/%s-exported_%s_%s*',eventdir,dir_infix,pGUID,eventname);
eprimedir_list=dir(dpat);
eprimedir=[];
if isempty(eprimedir_list)
  if isempty(extra_task_infix), return; end;
else
  eprime_info.eprime_folder_found = 1;
  if ~isempty(sesstype)
    for folderid=1:length(eprimedir_list)
      if regexp(eprimedir_list(folderid).name,sesstype)
        eprimedir = sprintf('%s/%s',eventdir,eprimedir_list(folderid).name);
        break;
      end
    end
  end
  if isempty(eprimedir)
    eprimedir = sprintf('%s/%s',eventdir,eprimedir_list(1).name);
    folderid=1;
  end
end

% find extra task directory
if ~isempty(extra_task_infix)
  dpat = sprintf('%s/%s-exported_%s_%s*',eventdir,extra_dir_infix,pGUID,eventname);
  extra_eprimedir_list=dir(dpat);
  extra_eprimedir=[];
  if ~isempty(extra_eprimedir_list)
    eprime_info.eprime_extra_folder_found = 1;
    if ~isempty(sesstype)
      for folderid_extra=1:length(extra_eprimedir_list)
        if regexp(extra_eprimedir_list(folderid_extra).name,sesstype)
          extra_eprimedir = sprintf('%s/%s',eventdir,extra_eprimedir_list(folderid_extra).name);
          break;
        end
      end
    end
    if isempty(extra_eprimedir)
      extra_eprimedir = sprintf('%s/%s',eventdir,extra_eprimedir_list(1).name);
      folderid_extra=1;
    end
  end;
end;

for i=1:length(extlist)
  ext = extlist{i};
  
  if eprime_info.eprime_folder_found && ~eprime_info.eprime_file_found
    fname_corr = sprintf('%s/%s_%s%s',eprimedir,pGUID,task_infix,ext);
    dlist_pguid_match = dir(sprintf('%s/*%s*%s',...
                            eprimedir,regexprep(pGUID,'^NDAR_INV',''),ext));
    dlist_ncorr = dir(sprintf('%s/*%s',eprimedir,ext));
%    dlist_pguid_match = dir(sprintf('%s/*%s*%s*%s',...
%                            eprimedir,regexprep(pGUID,'^NDAR_INV',''),task_infix,ext));
    %dlist_ncorr = dir(sprintf('%s/*%s*%s',eprimedir,task_infix,ext));
    % check for file with correct naming convention
    if exist(fname_corr,'file')
      eprime_info.eprime_file_found = 1;
      eprime_info.eprime_pguid_match = 1;
      eprime_info.eprime_file_naming  = 1;
      eprime_info.eprime_file_name = fname_corr;  
    % check if pguid is part of the filename
    elseif ~isempty(dlist_pguid_match)
      eprime_info.eprime_file_found = 1;
      eprime_info.eprime_pguid_match = 1;
      
      % Octavio (2017aug10): ---------------------------------------------
      % Pick "best" file name
% %       eprime_info.eprime_file_name = fullfile(eprimedir,dlist_pguid_match(1).name);
      
      file_sel = file_name_filter(dlist_pguid_match);
      eprime_info.eprime_file_name = fullfile(eprimedir,file_sel.name);
      % ---------------------------------------------: Octavio (2017aug10)
      
      
    % check for file with incorrect naming convention
    elseif ~isempty(dlist_ncorr)
      eprime_info.eprime_file_found = 1;
      
      % Octavio (2017aug10): ---------------------------------------------
      % Pick "best" file name
% %       eprime_info.eprime_file_name = fullfile(eprimedir,dlist_ncorr(1).name);
      
      file_sel = file_name_filter(dlist_ncorr);
      eprime_info.eprime_file_name = fullfile(eprimedir,file_sel.name);
      % ---------------------------------------------: Octavio (2017aug10)
      
    end

    %try other folders
    if ~eprime_info.eprime_file_found
      for j=1:length(eprimedir_list)
        %skip sesstype matched 
        if j==folderid, continue; end;

        eprimedir = sprintf('%s/%s',eventdir,eprimedir_list(j).name);

        fname_corr = sprintf('%s/%s_%s%s',eprimedir,pGUID,task_infix,ext);
        dlist_pguid_match = dir(sprintf('%s/*%s*%s',...
                                eprimedir,regexprep(pGUID,'^NDAR_INV',''),ext));
        dlist_ncorr = dir(sprintf('%s/*%s',eprimedir,ext));
    %    dlist_pguid_match = dir(sprintf('%s/*%s*%s*%s',...
    %                            eprimedir,regexprep(pGUID,'^NDAR_INV',''),task_infix,ext));
        %dlist_ncorr = dir(sprintf('%s/*%s*%s',eprimedir,task_infix,ext));
        % check for file with correct naming convention
        if exist(fname_corr,'file')
          eprime_info.eprime_file_found = 1;
          eprime_info.eprime_pguid_match = 1;
          eprime_info.eprime_file_naming  = 1;
          eprime_info.eprime_file_name = fname_corr;
        % check if pguid is part of the filename
        elseif ~isempty(dlist_pguid_match)
          eprime_info.eprime_file_found = 1;
          eprime_info.eprime_pguid_match = 1;
          
          % Octavio (2017aug10): ---------------------------------------------
          % Pick "best" file name
% %           eprime_info.eprime_file_name = fullfile(eprimedir,dlist_pguid_match(1).name);

          file_sel = file_name_filter(dlist_pguid_match);
          eprime_info.eprime_file_name = fullfile(eprimedir,file_sel.name);
          % ---------------------------------------------: Octavio (2017aug10)
          
          
        % check for file with incorrect naming convention
        elseif ~isempty(dlist_ncorr)
          eprime_info.eprime_file_found = 1;
          
          % Octavio (2017aug10): ---------------------------------------------
          % Pick "best" file name
% %           eprime_info.eprime_file_name = fullfile(eprimedir,dlist_ncorr(1).name);
          
          file_sel = file_name_filter(dlist_ncorr);
          eprime_info.eprime_file_name = fullfile(eprimedir,file_sel.name);
          % ---------------------------------------------: Octavio (2017aug10)
          
      
        end
        if eprime_info.eprime_file_found, break; end;
      end
    end
  end

  if ~isempty(extra_task_infix) && eprime_info.eprime_extra_folder_found && ~eprime_info.eprime_extra_file_found
    fname_corr = sprintf('%s/%s_%s%s',extra_eprimedir,pGUID,extra_task_infix,ext);
    dlist_pguid_match = dir(sprintf('%s/*%s*%s',...
                            extra_eprimedir,regexprep(pGUID,'^NDAR_INV',''),ext));
    dlist_ncorr = dir(sprintf('%s/*%s',extra_eprimedir,ext));
%    dlist_pguid_match = dir(sprintf('%s/*%s*%s*%s',...
%                            extra_eprimedir,regexprep(pGUID,'^NDAR_INV',''),extra_task_infix,ext));
    %dlist_ncorr = dir(sprintf('%s/*%s*%s',extra_eprimedir,extra_task_infix,ext));
%    if ~exist(fname_corr,'file') && isempty(dlist_pguid_match) && isempty(dlist_ncorr)
%      fname_corr = sprintf('%s/%s_%s%s',extra_eprimedir,pGUID,extra_task_infix_2,ext);
%      dlist_pguid_match = dir(sprintf('%s/*%s*%s*%s',...
%                              extra_eprimedir,regexprep(pGUID,'^NDAR_INV',''),extra_task_infix_2,ext));
%      %dlist_ncorr = dir(sprintf('%s/*%s*%s',extra_eprimedir,extra_task_infix_2,ext));
%      dlist_ncorr = dir(sprintf('%s/*%s',extra_eprimedir,ext));
%    end
    % check for file with correct naming convention
    if exist(fname_corr,'file')
      eprime_info.eprime_extra_file_found = 1;
      eprime_info.eprime_extra_pguid_match = 1;
      eprime_info.eprime_extra_file_naming  = 1;
      eprime_info.eprime_extra_file_name = fname_corr;  
    % check if pguid is part of the filename
    elseif ~isempty(dlist_pguid_match)
      eprime_info.eprime_extra_file_found = 1;
      eprime_info.eprime_extra_pguid_match = 1;
      eprime_info.eprime_extra_file_name = fullfile(extra_eprimedir,dlist_pguid_match(1).name);
    % check for file with incorrect naming convention
    elseif ~isempty(dlist_ncorr)
      eprime_info.eprime_extra_file_found = 1;
      eprime_info.eprime_extra_file_name = fullfile(extra_eprimedir,dlist_ncorr(1).name);
    end
   
    if ~eprime_info.eprime_extra_file_found
      for j=1:length(extra_eprimedir_list)
        %skip sesstype matched 
        if j==folderid_extra, continue; end;

        extra_eprimedir = sprintf('%s/%s',eventdir,extra_eprimedir_list(j).name);
        fname_corr = sprintf('%s/%s_%s%s',extra_eprimedir,pGUID,extra_task_infix,ext);
        dlist_pguid_match = dir(sprintf('%s/*%s*%s',...
                                extra_eprimedir,regexprep(pGUID,'^NDAR_INV',''),ext));
        dlist_ncorr = dir(sprintf('%s/*%s',extra_eprimedir,ext));
    %    dlist_pguid_match = dir(sprintf('%s/*%s*%s*%s',...
    %                            extra_eprimedir,regexprep(pGUID,'^NDAR_INV',''),extra_task_infix,ext));
        %dlist_ncorr = dir(sprintf('%s/*%s*%s',extra_eprimedir,extra_task_infix,ext));
    %    if ~exist(fname_corr,'file') && isempty(dlist_pguid_match) && isempty(dlist_ncorr)
    %      fname_corr = sprintf('%s/%s_%s%s',extra_eprimedir,pGUID,extra_task_infix_2,ext);
    %      dlist_pguid_match = dir(sprintf('%s/*%s*%s*%s',...
    %                              extra_eprimedir,regexprep(pGUID,'^NDAR_INV',''),extra_task_infix_2,ext));
    %      %dlist_ncorr = dir(sprintf('%s/*%s*%s',extra_eprimedir,extra_task_infix_2,ext));
    %      dlist_ncorr = dir(sprintf('%s/*%s',extra_eprimedir,ext));
    %    end
        % check for file with correct naming convention
        if exist(fname_corr,'file')
          eprime_info.eprime_extra_file_found = 1;
          eprime_info.eprime_extra_pguid_match = 1;
          eprime_info.eprime_extra_file_naming  = 1;
          eprime_info.eprime_extra_file_name = fname_corr;
        % check if pguid is part of the filename
        elseif ~isempty(dlist_pguid_match)
          eprime_info.eprime_extra_file_found = 1;
          eprime_info.eprime_extra_pguid_match = 1;
          eprime_info.eprime_extra_file_name = fullfile(extra_eprimedir,dlist_pguid_match(1).name);
        % check for file with incorrect naming convention
        elseif ~isempty(dlist_ncorr)
          eprime_info.eprime_extra_file_found = 1;
          eprime_info.eprime_extra_file_name = fullfile(extra_eprimedir,dlist_ncorr(1).name);
        end
      end
    end
  end;
  %if regexp(pGUID,'NDAR_INV4LUBR4AE'), keyboard; end;
  
  %if eprime_info.eprime_file_found, break; end;
end;


% Octavio (2017jul24-sep07): ---------------------------------------------
% Append behavioral-file date and time to output record
%
eprime_info.eprime_datetime_found = 0;
eprime_info.eprime_file_res = 0;
eprime_info.eprime_date = NaN;
eprime_info.eprime_time = NaN;

if eprime_info.eprime_file_found
    try
        [date,time,dtres,fres] = eprime_datetime_get( eprime_info.eprime_file_name );

        if dtres  % We got valid date and time, so:
            eprime_info.eprime_datetime_found = dtres;
            eprime_info.eprime_file_res = fres;
            eprime_info.eprime_date = str2num(date);
            eprime_info.eprime_time = str2num(time);
        end
    catch err
        fprintf('\ndpat = %s\n', dpat );
        eprime_info
        if isfield(eprime_info,'eprime_file_name')
            fprintf('eprime_info.eprime_file_name =\n%s\n', eprime_info.eprime_file_name );
        end
        fprintf('%s\n', err.message);
        fprintf('Unable to get eprime date & time\n');
    end
end
%---------------------------------------------- :Octavio (2017jul24-sep07)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Octavio (2017sep21) ----------------------------------------------------
% removed from this script and put it in abcd_check_eprime_all:
% % % function sitename = check_sitename(sitename)
% % %   sitename = upper(sitename);
% % %   switch sitename
% % %     case 'DAIC'
% % %       sitename = 'UCSD';
% % %     case 'WASHU'
% % %       sitename = 'WUSTL';
% % %     case 'OAHU'
% % %       sitename = 'UMB';
% % %   end;
% % % return;
%---------------------------------------------------- :Octavio (2017sep21)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eprime_info = init_eprime_info(taskname)
  eprime_info.taskname = taskname;
  eprime_info.eprime_folder_found = 0;
  eprime_info.eprime_file_found = 0;
  eprime_info.eprime_pguid_match = 0;
  eprime_info.eprime_file_naming = 0;
  eprime_info.eprime_file_name = '';
  eprime_info.eprime_extra_folder_found = 0;
  eprime_info.eprime_extra_file_found = 0;
  eprime_info.eprime_extra_pguid_match = 0;
  eprime_info.eprime_extra_file_naming = 0;
  eprime_info.eprime_extra_file_name = '';

  % Octavio (2017jul25-aug10): -------------------------------------------
  eprime_info.eprime_datetime_found = 0;
  eprime_info.eprime_file_res = 0;
  eprime_info.eprime_date = NaN;
  eprime_info.eprime_time = NaN;
  %-------------------------------------------- :Octavio (2017jul25-aug10)

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio (2017jul25-aug16):

function sel = file_name_filter(list0)
% Select "best" ABCD file name from list.  Octavio Ruiz, 2017aug10-16
% 
% Favor *.txt. over *.csv.  Remove problematic file names like:
%   containing a space and a number between parentheses.
%   containing questions marks
% 
% Input is a struct array with fields:
%     name
%     date
%     bytes
%     isdir
%     datenum

if ~isempty(list0)
    if length(list0) > 1
        listn = list0;

        % Check if .txt files are listed, if so remove the others
        ss = regexp( {listn.name}, '.txt' );
        ind = cellfun( @isempty, ss );
        listn(ind) = [];
        
        if length(listn) > 1
            % Check for spaces in filename, if found, remove corresponding entries
            ss = regexp( {listn.name}, '\s' );
            ind = cellfun( @isempty, ss );
            listn(ind) = [];
            
            if length(listn) > 1
                % Check for ? in filename, if found, remove corresponding entries
                ss = regexp( {listn.name}, '?' );
                ind = cellfun( @isempty, ss );
                listn(ind) = [];
                fprintf('\n- abcd_check_eprime: file_name_filter: "?" found in file-name list -');
                list0;
                listn;
            end
        end
        if length(listn) >= 1
            sel = listn(1);
        else
            % All fnames filtered out, return first entry in original list
            sel = list0(1);
        end
    else
        % There is only one file name in list, return it
        sel = list0;
    end
else
    fprintf('- Empty file-name list -');
    sel = list0;
end
%                                                         :Octavio (2017aug16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
