function badchan_i=ts_read_txt_badchans(fname,labels)
% ts_read_txt_badchans: reads text file containing bad channel labels
%                    compares to list of labels
%                    returns indices of bad channels
%
% Usage: badchan_i=ts_read_txt_badchans(fname,labels)
%
% Input:
%   fname: name of text file containing bad channel labels
%   labels: cell array of all channel labels
%           labels are stored in hdr.sensors.label, created by read_fif_header
%
% Output:
%   badchan_i: vector of channel indices corresponding to bad channels
%
% See also read_fif_header, browseraw, avg_fif_data
%
% Created:    05/10/06 by Don Hagler
% Last Mod:   08/05/09 by Don Hagler
%

badchan_i = [];

if (~mmil_check_nargs(nargin,2)) return; end;

if ~exist(fname,'file')
  error('file %s not found',fname);
end

% make sure labels is a cell array
if ~iscell(labels) | isempty(labels)
  fprintf('%s: labels must be a cell array of strings\n',mfilename);
  return;
end;
numlabels = length(labels);

% open bad chans file
fid=fopen(fname,'rt');
j = 1;
while (~feof(fid))
  badchans{j} = fgetl(fid);
  j=j+1;
end
fclose(fid);

if isempty(badchans) || (isnumeric(badchans{1}) && badchans{1}==1)
  fprintf('%s: no bad chans found in %s\n',mfilename,fname);
  return;
end;

[tmp_labels,tmp_i,badchan_i]=intersect(badchans,labels);
badchan_i = sort(badchan_i);

return;
