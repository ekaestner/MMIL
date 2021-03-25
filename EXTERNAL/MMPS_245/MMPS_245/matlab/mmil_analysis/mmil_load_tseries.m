function [data,roinames,time] = mmil_load_tseries(fname)
%function [data,roinames,time] = mmil_load_tseries(fname)
%
% Purpose: load time series from csv file (output of mmil_extract_tseries)
%
% Usage:
%  output = mmil_load_tseries(fname, 'key1', value1,...);
%
% Required Input:
%  fname: full path name of csv file containing time vector
%     and time series values for ROIs in additional columns
%     First column should have time, additional columns should contain
%       extracted values
%     First row should be column headers, number of additional rows
%       is number of time points
%
% Output:
%   data: matrix of values; size = [ntpoints,nrois]
%   roinames: cell array of roi names (column headers)
%   time: vector of time points (first column values)
%
% Created:  02/29/12 by Don Hagler
% Last Mod: 10/29/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = []; data = []; roinames = [];
if ~mmil_check_nargs(nargin, 1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read tseries file

if ~exist(fname,'file'), error('file %s not found',fname); end;
try
  input_info = mmil_readtext(fname);
  time = cell2mat({input_info{2:end,1}});
catch
  error('failed to get time from tseries file %s -- check format',fname);
end;
roinames = input_info(1,2:end);
try
  data = cell2mat(input_info(2:end,2:end));
catch
  error('list of regressors in file %s must all be integers',fname);
end

