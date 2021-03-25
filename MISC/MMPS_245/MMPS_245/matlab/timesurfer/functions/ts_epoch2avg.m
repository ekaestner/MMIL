function avg_data = ts_epoch2avg (epoch_data,varargin)
% function avg_data = ts_epoch2avg (epoch_data)
%
% Use: avg_data = ts_epoch2avg (epoch_data)
%
% Averages epochs across trials and return the appropriate
% data structure.
%
% Required Input:
%
% epoch_data - valid TimeSurfer epoch_data structure
%
% Output:
%
% avg_data - TimeSurfer avg_data sructure
%
% Created by:       Rajan Patel 03/31/2008
% Last Modified by: Rajan Patel 03/31/2008

%% todo: calculate ncov

%% Check Inputs

if nargin < 1, help(mfilename); end

errors = ts_checkdata(epoch_data);
if ~isempty(errors)
    mmil_error(parms,'Errors in provided data structure: %s.',errors);
end                   

%% Create avg_data

avg_data = [];

avg_data          = epoch_data;
avg_data.averages = avg_data.epochs;
avg_data          = rmfield(avg_data,'epochs');

for c = 1:length(avg_data.averages)
 avg_data.averages(c).stdev    = squeeze( std(avg_data.averages(c).data,0,3));
 avg_data.averages(c).data     = squeeze(mean(avg_data.averages(c).data,3));
end

%% Check output

errors = ts_checkdata(avg_data);
if ~isempty(errors)
    mmil_error(parms,'Errors in created avg_data structure: %s.',errors);
end 
