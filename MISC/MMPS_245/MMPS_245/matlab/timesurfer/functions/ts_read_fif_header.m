function hdr = ts_read_fif_header(headerfile,get_eof_flag)
%function hdr = ts_read_fif_header(headerfile,get_eof_flag)
%
% This returns a header structure with the following elements
%   hdr.sfreq           sampling frequency
%   hdr.nChans          number of channels
%   hdr.nSamples        number of samples per buffer
%   hdr.skips           vector of sample latencies for "skips"
%   hdr.tfirst          first time point in file
%   hdr.tlast           last time point in file
%   hdr.nBuffs          number of buffers in file (not including skips)
%   hdr.sensors         structure containing sensor information
%
% The following elements contained in hdr.sensor come directly from fif file:
%  hdr.sensors.label     array with labels of each channel
%  hdr.sensors.kind      channel kinds
%    1 = MEG, 2 = EEG, 202 = EOG, 3 = STI (stim triggers)
%  hdr.sensors.lognum    logical numbers
%    contains same numbers found in channel labels
%  hdr.sensors.type      coil types
%    3012 = gradiometer, 3022 = magnetometer, 1 = EEG, 5 = EOG, 0 = STI
%  hdr.sensors.loc       cell-array containing 4x4 position matrices
%
% The following element is generated based on the previous elements
%  hdr.sensors.typestring   channel type
%    possible values: 'mag', 'grad1', 'grad2', 'eeg', 'eog', 'sti', 'other'
%
%
% get_eof_flag is optional (default value is 0)
%   if get_eof_flag = 1, will read entire file to get nBuffs
%   if get_eof_flag = 0, will set nBuffs = -1
%
% See also ts_read_fif_data, ts_read_fif_events
%
% based on FieldTrip's read_fcdc_header
% Copyright (C) 2003-2005, F.C. Donders Centre
%
% created:  03/10/06 by Don Hagler
% last mod: 10/24/11 by Don Hagler
%

%% todo: reformat help
%% todo: use MNE matlab toolbox instead of fiff access

hdr = [];
if (~mmil_check_nargs(nargin,1)) return; end;
if ~exist('get_eof_flag','var') | isempty(get_eof_flag), get_eof_flag = 0; end;

% test whether the file exists
if ~exist(headerfile,'file')
  error('file %s not found', headerfile);
end

% read sensor info
[hdr.sensors.label,hdr.sensors.kind,hdr.sensors.lognum,hdr.sensors.type,...
  hdr.sensors.loc]=channames(headerfile);

if isempty(hdr.sensors.label)
    fprintf('%s: WARNING: no sensor info found in fif file.', mfilename);
    hdr.sensors.typestrings = [];
else
  % generate type strings
  hdr.nChans = length(hdr.sensors.label);
  grad = find(hdr.sensors.type == 3012);
  grad1 = grad(1:2:length(grad));
  grad2 = grad(2:2:length(grad));
  for i=1:hdr.nChans
    switch hdr.sensors.type(i)
      case 3012,
        if ismember(i,grad1)
          hdr.sensors.typestring{i} = 'grad1';
        else
          hdr.sensors.typestring{i} = 'grad2';
        end;
      case {3013},
        hdr.sensors.typestring{i} = 'grad1';
      case {3022 3023 3024 4001},
        hdr.sensors.typestring{i} = 'mag';
      case 1,
        hdr.sensors.typestring{i} = 'eeg';
      case 5,
        if hdr.sensors.kind(i) == 202 | ...
           hdr.sensors.kind(i) == 2
          hdr.sensors.typestring{i} = 'eog';
        elseif hdr.sensors.kind(i) == 402
          hdr.sensors.typestring{i} = 'ecg';
        else
          hdr.sensors.typestring{i} = 'emisc';
        end;        
      case 0,
        hdr.sensors.typestring{i} = 'sti';
      otherwise,
        hdr.sensors.typestring{i} = 'other';
      end;
  end
end;
hdr.sensors.typestring = hdr.sensors.typestring';

% open file
rawdata('any',headerfile);
% get calibration values for pre-scaling channels
[hdr.sensors.range, hdr.sensors.cal] = rawdata('range');
hdr.sfreq = rawdata('sf');

% read first buffer
hdr.tfirst = rawdata('t');
[buf, status] = rawdata('next');
hdr.nChans   = size(buf,1);
hdr.nSamples = size(buf,2);
hdr.skips = [];

% read to end of file to determine file length
if(get_eof_flag)
  hdr.nBuffs = 0;
  while ~strcmp(status, 'eof')
    if strcmp(status,'error'),
      lasterr('File error');
      error('File error');
    elseif strcmp(status,'ok'),
      hdr.nBuffs = hdr.nBuffs + 1;
    elseif strcmp(status, 'skip')
      hdr.skips = [hdr.skips rawdata('t')*hdr.sfreq];
    end;
    [buf, status] = rawdata('next');
  end
  hdr.tlast = rawdata('t')-0.1;
else
  hdr.tlast = -1;
  hdr.nBuffs = -1;
end
rawdata('close');

