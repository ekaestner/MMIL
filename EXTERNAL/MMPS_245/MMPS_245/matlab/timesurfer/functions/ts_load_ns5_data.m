function [data,sfreq] = ts_load_ns5_data(fname,downsampling_factor)
%function [data,sfreq] = ts_load_ns5_data(fname,downsampling_factor)
%
% Purpose: load and downsample (channel-by-channel), position
%   and save an NS5 file as Utah array datacube
%
% Required Input:
%   fname: full path name of input file
%
% Optional Input:
%   downsampling_factor: downsampling factor (integer)
%     {default = 30}
%
% Output:
%   data: data matrix with size = [nchans,ntpoints]
%   sfreq: sampling frequency (after downsampling)
%
% Notes: 1. uses NPMK toolbox function "openNSx"
%        2. anti-aliasing factor set to 3
%           e.g. 30 kHz signal with downsampling factor of 10 will result
%                in a 3 kHz signal with a cutoff frequency at 1 kHz
%
% Created:  10/09/14 by Lyle Muller (as process_ns5_file)
% Last Mod: 01/05/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = []; sfreq = [];

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('downsampling_factor','var') || isempty(downsampling_factor)
  downsampling_factor = 30;
end;

% parameters
anti_alias_factor = 3;
filter_order = 4;

assert(ceil(downsampling_factor) == floor(downsampling_factor),...
       'downsampling_factor must be integer' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load header
nsx_header = openNSx(fname);

% initialize data array
ntpoints = nsx_header.MetaTags.DataPoints;
nchans = nsx_header.MetaTags.ChannelCount;
sfreq = nsx_header.MetaTags.SamplingFreq;
data = nan(nchans,ceil(ntpoints/downsampling_factor));

% calculate filter coefficients
[b,a] = butter(filter_order,...
    ((sfreq/(downsampling_factor*anti_alias_factor))/(sfreq/2)),'low');

% load and process NSx file channel-by-channel
for ii = 1:nchans
  fprintf('%s: loading channel %d...\n',mfilename,ii);
  channel_string = sprintf('c:%d:%d',ii,ii);
  nsx_channel = openNSx(fname,...
    'read','report','precision','double',channel_string);
  % data verification check; assert all equal length
  assert(length(nsx_channel.Data) == ntpoints)
  % antialiasing filter and downsampling
  xf = filtfilt(b,a,nsx_channel.Data);
  channel_out = downsample(xf,downsampling_factor);
  % channel positioning
  data(ii,:) = channel_out;
end

sfreq = sfreq/downsampling_factor;

