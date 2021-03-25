function data = ts_MNE_loadfif(datafile,varargin)
% ts_MNE_loadfif(file,'begsample',0,'endsample',10000,'chantype','eeg');

parms = mmil_args2parms( varargin, ...
                         { 'channels',[],[],...
                           'chantype',[],[],...
                           'begsample' ,0,[],...
                           'endsample',[],[],...
                         }, ...
                         false );
                       
%                            'continuous',0,[],...
%                            'stim_delay',0,[],...
%                            'valid_event_codes',[],[],...
%                            'trigchan',[],[],...
%                            'trigthresh',[],[],...
%                            'trig_minduration',[],[],...
%                            'trig_minseparation',.01,[],...
%                            'oldchanlabels',[],[],...
%                            'newchanlabels',[],[],...
%                            'evntfile',[],[],...
                           
[FIFPATH,FIFFILE,EXT] = fileparts(datafile);
  % check EXT = FIF
FIFFILE = [FIFFILE EXT];
FILE    = datafile;%fullfile(FIFPATH,FIFFILE);
chans   = parms.channels;
typestr = parms.chantype;
begsamp = parms.begsample;
endsamp = parms.endsample;

if ~exist(datafile,'file')
  fprintf('File not found: %s\n',datafile);
  return;
end

cwd = pwd;
if exist(FIFPATH,'dir')
  cd(FIFPATH);
end

try
  hdr   = ts_read_fif_header(FILE,0);
  for k = 1:length(hdr.sensors.label)
    sens(k).label         = hdr.sensors.label{k};
    sens(k).typestring    = hdr.sensors.typestring{k};
    sens(k).type          = hdr.sensors.type(k);
    sens(k).kind          = hdr.sensors.kind(k);
    sens(k).badchan       = 0;
    sens(k).lognum        = hdr.sensors.lognum(k);
    sens(k).loc           = hdr.sensors.loc{k};
  end
  if isempty(chans) && ~isempty(typestr)
    chans = strmatch(typestr,{sens.typestring});
    sens  = sens(chans);
  end

  hdr   = fiff_setup_read_raw(FIFFILE);
  if isempty(endsamp)
    % set endsamp to total number of samples
    endsamp = hdr.last_samp;
  end
  if isempty(chans)
    [x,t] = fiff_read_raw_segment(hdr,begsamp,endsamp);
  else
    [x,t] = fiff_read_raw_segment(hdr,begsamp,endsamp,chans);
  end
  try
    trans = loadtrans(datafile);
  catch
    trans = []; % there is no trans in empty room fif data
  end
  data    = ts_matrix2epoch(x,'continuous',1,'sens',sens,'time',t,'device2head',trans);
  clear hdr x t sens chans
  cd(cwd);
catch
  rethrow(lasterror);
  cd(cwd);
end