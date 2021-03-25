function [subject, date, task, condition, data_info, data_proc] = ts_iEEG_decon_fname (file_name)

% [subject, date, task, condition, data_info,data_procs] = ts_iEEG_decon_fname(file_name)
%
% Function to deconstruct file names of the agreed upon naming structure
% of iEEG data files.
%
% Standard information about the data will be separated by '_' (eg. subject
% information, data type, ect. and that information about preprocessing
% steps applied will be separated by '.' (eg. filters applied ect.).
%
% The first two parts of the filename must be the subject-date_task.
% The end of the filename must be condition number preceded by '.'.
%
% The function behaves that as soon as the first '.' is econuntered in the file
% name the rest of the information will automatically be interpretted as
% data preprocessing information.  And subsequent '_' will be treated as
% being a part of any previous '.' string section.
%
% The file extension will not be returned.
%
% In the case of file_names being a cell array of strings, if there are
% file names of differing length then the shorter file names will have
% empty cells in the returned cell array of strings.
%
% Input:
%   file_name: can be a single string or cell array of strings
%
% Outputs:
%   Subject 
%   Date
%   Task
%   Condition
%   data_info: Cell array of strings containing the further data info.
%   data_proc: Cell array of strings containing the preprocessing steps.
%
%     both arrays will have multiple rows (dimension 1) if there are
%     multiple file_names is a cell array of strings.
%
% Created: 08/22/2007 Rajan Patel
%
%


if iscell(file_name), file_name=char(file_name); end;
[path,file_name] = fileparts(file_name);
clear path
data_info = {};
data_proc = {};
subject   = {};
date      = {};
task      = {};
condition = {};

for j = 1:size(file_name,1)
    file_name(j,end+1) = '.';
    curr_char = 1;
    curr_info = 1;
    curr_info_str = [];
    while (file_name(j,curr_char) ~= '.') && (curr_char<size(file_name,2))
        if file_name(j,curr_char)=='_'
            if curr_info==1
                subject(j) = {curr_info_str(1:strfind(curr_info_str,'-')-1)};
                date(j)    = {curr_info_str(strfind(curr_info_str,'-')+1:length(curr_info_str))};
            elseif curr_info==2
                task(j) = {curr_info_str};
            else
                data_info(j,curr_info-2) = {curr_info_str};
            end
            curr_info_str=[];
            curr_info=curr_info+1;
            curr_char=curr_char+1;
        else
            curr_info_str= [curr_info_str file_name(j,curr_char)];
            curr_char=curr_char+1;
        end
    end
    if curr_info==1
        subject(j) = {curr_info_str(1:strfind(curr_info_str,'-')-1)};
        date(j)    = {curr_info_str(strfind(curr_info_str,'-')+1:length(curr_info_str))};
    elseif curr_info==2, task(j) = {curr_info_str};
    elseif curr_info>2, data_info(j,curr_info-2) = {curr_info_str};
    else data_info(j,1) = {''}; end

    %% Extract data processing options
    curr_char = curr_char + 1;
    curr_proc = 1;
    curr_proc_str=[];
    data_proc(j,curr_proc)={curr_proc_str};
    while curr_char <= size(file_name,2)
        if file_name(j,curr_char)=='.'
            data_proc(j,curr_proc) = {curr_proc_str};
            curr_proc_str=[];
            curr_proc=curr_proc+1;
            curr_char=curr_char+1;
        else
            curr_proc_str= [curr_proc_str file_name(j,curr_char)];
            curr_char=curr_char+1;
        end
    end
    if curr_proc>1
        if findstr(char(data_proc(j,curr_proc-1)),'cond')
            condition(j) = data_proc(j,curr_proc-1);
        end
    end;
    %     data_proc(j,curr_proc-1) = '';
end
