function ts_iEEG_recodestats (filename,varargin)
% Usage: ts_iEEG_recodestats (filename,'option',value...
%
% To recode a between statistics with a new event code and saving it out.
%
% Required input:
%
% filename       - the name of the original statistics file
% 'neweventcode' - the event code to recode the between stats to
% 
% Created 05/29/2008 by Rajan Patel
% Last Modified 05/28/2008 by Rajan Patel

%% Check Inputs

if nargin < 1, help(mfilename); end

opt = mmil_args2parms(varargin,...
                     {'neweventcode',[],[],...
                      'verbose',true,sort([false true]),...
                      'logfile',[],[],...
                      'logfid',1,[],...
                      },...
                      false);      
   
if ~exist(filename,'file')
    mmil_error(opt,'Could not find file: %s.',filename);
else
    load(filename);
end

if ~exist('stat_data','var')
    mmil_error(opt,'Could not find stat_data in file: %s.',filename);
end
                  
if isempty(opt.neweventcode)
    mmil_error(opt,'You must provide a new event code.');
end

[filepath,filename,fileext] = fileparts(filename);

filename = filename(1:strfind(filename,'event')+4);

filename = [filename num2str(opt.neweventcode) fileext];

filename = fullfile(filepath,filename);

if exist(filename,'file')
    mmil_error(opt,'The target file already exists, check your neweventcode: %s.',filename);
end

%% 

indx = find([stat_data.stats.event_code] == -1);

if ~isempty(indx)
    stat_data.stats(indx).event_code = opt.neweventcode;
    stat_data.stats                   = stat_data.stats(indx);
else
    mmil_error(opt,'Could not find statistics between conditions in the given data.');
end

mmil_logstr(opt,'Saving new stat_data: %s.',filename);
save(filename,'stat_data');
