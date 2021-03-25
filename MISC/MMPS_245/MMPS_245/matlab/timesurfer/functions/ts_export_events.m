function ts_export_events(evnts,fname)
%function ts_export_events(evnts,[fname])
%
% Purpose: export events in matlab structure to tab-delimited text file
%
% Required Input:
%   evnts - structure containing events -- see ts_read_fif_events
%
% Optional Parameters:
%   fname - output file name
%     {default: 'exported_events.txt'}
%
% created: 09/07/06 Don Hagler
% last modified: 09/07/06 Don Hagler

if nargin < 1
   help(mfilename);
   return;
end

if isempty(evnts)
  fprintf('%s: error: evnts structure is empty\n',mfilename);
  return;
end;

if ~exist('fname','var')
  fname = 'exported_events.txt';
end;

fid = fopen(fname,'wt');
if fid==-1
  fprintf('%s: error creating output file %s\n',mfilename,fname);
  return;
end;
fprintf(fid,'type\tlatency\tcondition\tduration\n');
for j=1:length(evnts)
  if isempty(evnts(j).condition)
    evnts(j).condition=0;
  end;
  fprintf(fid,'%s\t%d\t%d\t%d\n',...
    evnts(j).type,evnts(j).latency,evnts(j).condition,evnts(j).duration);
end;
fclose(fid);

return;
