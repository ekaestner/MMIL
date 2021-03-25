function [date,time,dtres,fres] = eprime_datetime_get(fname,varargin)
% function [date,time,dtres,fileres] = eprime_datetime_get(fname,varargin)
%
% Purpose: extract timestamp from eprime data file
% Based partly on abcd_extract_eprime_mid.m
% 
% Major modifications by Octavio Ruiz (2017sep-) to use his python script
% to extract data from E-Prime spreadsheet files
%
% Input (required):
%   fname: name of input eprime file
% Output:
%   date
%   time
%   dtres   Result of date&time conversion: 0,1
%   fres    Issues found while reading and interpreting file
%
% Created:  10/07/16  by Don Hagler
% Last Mod: 2017aug04-oct20 by Octavio Ruiz
%

date = '';   time = '';   dtres = 0;   fres = 10000;  % Unable to read
fname_contns_parenthes = false;

if ~mmil_check_nargs(nargin,1), return; end;

% fname sometimes is a copy of the original file, signaled by ' (1)' etc. in its name
% Here I remove those ' (*)' in order to access the original file.
i1 = strfind(fname,' (');
if ~isempty(i1)
    fname_contns_parenthes = true;
    fprintf('\nRequested file name contains parentheses:\n%s\n', fname);
    i2 = strfind(fname,')');
    if ~isempty(i2) && i2 > i1
        % Remove one space, parentheses, and anything between them
        fname(i1:i2) = '';   % Attempted to read file without ' (*)' in file name
        fprintf('Name replaced with\n%s\n', fname);
    else
        fprintf('Unable to fix\n');
    end
end


% Get date and time from eprime file
date = '';
time = '';

% Use python script to read file and get task-start date and time
% % % cmnd0 = './eprime_sprdsht_get.py';
cmnd0 = 'python $MMPS_DIR/python/eprime_sprdsht_get.py';
cmdline = [cmnd0 ' ' fname ' DateTimeDiagnos'];
[status, cmdout] = system(cmdline);

% Parse returned information into date, time, and file diagnostics
if ~status
    [datime,fres_s] = strtok(cmdout, ',');
    fres = str2num(fres_s);
    try
        date = datestr( datime, 'yyyymmdd');
        time = datestr( datime, 'HHMMSS');
        % If arrived here, everything is right
        dtres = 1;
    catch err
        fprintf('%s\n- Unable to get datetime -\n', fname);
        fprintf('%s\n', err.message);
        date, time, fres
    end
else
    fres = 1000;
    fprintf('%s\n- Unable to read or interpret file -\n', fname);
end

if fname_contns_parenthes
    % Update fres to reflect presence of parentheses in file name
    fres = fres + 100000;  % diagnos.number added here must be different from those returned by eprime_sprdsht_get.py
end
