function evnts = ts_import_events(fname)
%function evnts = ts_import_events(fname)
%
% Purpose: import events from tab-delimited text file to matlab structure
%
% Required Input:
%   fname - input file name
%     This should be a tab-delimited text file with four columns:
%       type, latency, condition, and duration
%     The first row should contain labels and will be ignored
%     Acceptable "types" are:
%       trigger, manual, reject
%     "latency" is in samples relative to start of file
%     "condition" should be integer event codes
%     "duration" is the duration of the event
%        (duration should be 1 if this is not relevant to the analysis)
%
% Output:
%   evnts - structure containing events -- see ts_read_fif_events
%
% created:  09/12/06 Don Hagler
% Last Mod: 08/05/09 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

fid = fopen(fname,'rt');
if fid==-1
  error('failed to open input file %s',fname);
end;
C = textscan(fid,'%s\t%s\t%s\t%s\n');
fclose(fid);
nlines = length(C{1});

% skip first line
j = 1;
evnts = [];
for i=2:nlines
  evnts(j).type = C{1}{i};
  evnts(j).latency = str2num(C{2}{i});
  evnts(j).condition = str2num(C{3}{i});
  evnts(j).duration = str2num(C{4}{i});
  evnts(j).fname = fname;
  j=j+1;
end;

return;
