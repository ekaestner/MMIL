function [data,eofstatus] = ts_read_fif_data(datafile, hdr, begsample, endsample, chanindx)
%function [data,eofstatus] = ts_read_fif_data(datafile, hdr, begsample, endsample, chanindx)
%
% Pupose: read raw data from Neuromag fif files
%
% Usage:
%   [data,eofstatus] = ts_read_fif_data(datafile, hdr, begsample, endsample, chanindx);
%     hdr is obtained using ts_read_fif_header
%
% The beginsample and endsample are counted from the beginning of the raw
% data, starting at 1.
%
% The argument "chanindx" is optional, the default is to read all channels.
%
% See also ts_read_fif_header
%
% based on FieldTrip's read_fcdc_data
% Copyright (C) 2003-2005, F.C. Donders Centre
%
% created:  03/10/06 by Don Hagler
% last mod: 10/24/11 by Don Hagler
%

%% todo: use MNE matlab toolbox instead of fiff access

data = [];
eofstatus = 0;
if (~mmil_check_nargs(nargin,4)) return; end;
if ~exist('chanindx','var') | isempty(chanindx)
  chanindx = 1:hdr.nChans;
end

SKIP_LENGTH = 100;

% read untill the end of the file if the endsample is "inf"
if isinf(endsample) && endsample>0
  endsample = hdr.nSamples*hdr.nBuffs;
end

% test whether the file exists
if ~exist(datafile)
  error('file %s not found',datafile);
end

% test whether the requested channels can be accomodataed
if min(chanindx)<1 || max(chanindx)>hdr.nChans
  error('selected channels are not present in the data');
end

% test whether the requested data segment is not outside the file
if begsample<1
  error('cannot read data before the start of the file');
end

% calculate time corresponding to begsample
begtime = hdr.tfirst + (begsample-1)/hdr.sfreq;

% check that begtime is less than tlast, otherwise return zeros
if(begtime >= hdr.tlast && hdr.tlast > 0)
  data = zeros(length(chanindx),endsample-begsample+1);
  return;
end

% if fif data is > 2GB, it gets broken into multiple files
% in that case, tfirst does not = 0
first_samp = floor(hdr.tfirst*hdr.sfreq) + 1;
try_begsample = floor(hdr.tfirst*hdr.sfreq) + begsample;
try_endsample = try_begsample - begsample + endsample;

% open file
rawdata('any',datafile);
real_begtime = rawdata('goto', begtime);
real_begsample = floor(real_begtime*hdr.sfreq) + 1;
samp_offset = try_begsample - real_begsample;

% must treat skips and data buffers differently to figure
% out how many buffers must be read
num_samples = try_endsample - real_begsample + 1;
skip_i = find(hdr.skips>=real_begsample-first_samp+1 & ...
              hdr.skips<=try_endsample-first_samp+1);
num_samples = num_samples - length(skip_i)*SKIP_LENGTH;
num_buffs = ceil(num_samples/hdr.nSamples) + length(skip_i);
num_samples = try_endsample - real_begsample + 1;

s1 = 0;
s2 = 0;
for i=1:num_buffs
  [M,status]=rawdata('next');
  if strcmp(status, 'eof')
    s1=s2+1;
    s2=num_samples;
    eofstatus = 1;
    fprintf('%s: reached end of file\n',mfilename);
    data(:,s1:s2) = ...
      zeros(length(chanindx),s2-s1+1);
    break;
  elseif strcmp(status,'error'),
    error('error reading selected data from fif-file');
  elseif strcmp(status,'skip'),
    s1=s2+1;
    s2=s1+SKIP_LENGTH-1;
    data(:,s1:s2) = ...
      zeros(length(chanindx),SKIP_LENGTH);
  elseif strcmp(status,'ok'),
    s1=s2+1;
    s2=s1+hdr.nSamples-1;
    data(:,s1:s2) = M(chanindx,:);
  end;
end;
rawdata('close');

% return only the data requested
endsample = samp_offset + endsample - begsample + 1;
begsample = samp_offset + 1;

data = data(:,begsample:endsample);
