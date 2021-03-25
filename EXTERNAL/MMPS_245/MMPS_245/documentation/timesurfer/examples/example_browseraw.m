rawfile = '/space/monkeys/1/home/dhagler/data/meg-raw/060325LP/060325LP_raw.fif';
matfile = 'matfiles/edited_events_1.mat';
badchanfile = '060325LP_badchans.txt';
dsfact = 4;
winlength = 5;
minchan = 1;
maxchan = 60;
report_progress = 1;

load(matfile);

ts_browseraw(rawfile,...
  'events',evnts,'hdr',hdr,...
  'read_mag',0,'read_grad1',1,'read_grad2',0,...
  'read_eeg',0,'read_eog',0,'read_other',0,...
  'autoreject',0,'filter',0,...
  'winlength',winlength,...
  'badchanfile',badchanfile,...
  'dsfact',dsfact,...
  'minchan',minchan,'maxchan',maxchan,...
  'report_progress',report_progress);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% tips on how to make ts_browseraw faster:
%
% 1. do not run remotely
% 2. shorter window length (e.g. 1-5 sec)
% 3. filter = 0
% 4. pre-filter and downsample raw data before using ts_browseraw (neuromag?)
% 5. use read_mag, read_grad1, etc. to select a type of channels
% 6. specify minchan and maxchan to select subset of channels
% 7. dsfact > 1
% 8. events = []
%
%
% if report_progress=1, you can see how long steps like reading,
%   filtering, and drawing take with different settings
%
%

