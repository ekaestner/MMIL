function first_stats = mmil_read_first_stats(fname)
%function first_stats = mmil_read_first_stats(fname)
%
% Purpose: read text file containing volumes derived from FSL's FIRST
%
% Required Input:
%   fname: full path file name of text file containing volumes
%
% Output:
%   first_stats is a struct array containing:
%     roiname
%     nvox
%     volume
%
% created:  01/09/07 by Don Hagler
% last mod: 10/11/07 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
first_stats = [];

if ~exist(fname,'file'), error('file %s not found',fname); end;
tmp_stats = mmil_readtext(fname, ' ');
nroi = size(tmp_stats,1);
nstats = size(tmp_stats,2);

if nstats~=4
  error('incorrect number of stats in %s (%d)',fname,nstats);
end;

for i=1:nroi
  first_stats(i).roiname = char(tmp_stats{i,2});
  first_stats(i).nvox = double(tmp_stats{i,3});
  first_stats(i).volume  = double(tmp_stats{i,4});
end;

