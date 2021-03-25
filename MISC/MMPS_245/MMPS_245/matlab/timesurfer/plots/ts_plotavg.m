function ts_plotavg(filename,varargin)
% Purpose: plot data in a Neuroscan *.avg file
% Inputs: *.avg filename or 2D average matrix ['signal' returned by loadavg()]
%
% Created by JSS on 21-May-2009

if ischar(filename)
  % load *.avg data
  [signal,variance,chan_names,pnts,rate,xmin,xmax] = loadavg(filename);
  % convert average matrix to avg_data structure
  dat = ts_matrix2avg(signal,'time',[xmin:1/rate:xmax]);  
elseif isnumeric(filename)
  dat = ts_matrix2avg(filename);  
end

% plot avg_data
ts_ezplot(dat,varargin{:});