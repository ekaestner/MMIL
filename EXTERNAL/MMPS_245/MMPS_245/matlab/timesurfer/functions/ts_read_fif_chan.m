function [B,nBuffs,tlast,skipevents,sf]=ts_read_fif_chan(filename,channels);
% ts_read_fif_chan   - Load some channels from raw data file
%
%    [B,NBUFFS,SKIPEVENTS,SF]=ts_read_fif_chan(FILENAME,CHANNELS);
%       CHANNELS can be either
%          order numbers: [1 2 3 10]
%          or names:      {'MEG 0111', 'MEG 0112', 'STI 001'}
%
%       if a "skip" is found, 100 samples will be inserted
%
%  Based on Neuromag meg_pd fiffaccess' rawchannels
%
% created:  03/11/06 by Don Hagler
% last mod: 10/24/11 by Don Hagler
%

%% todo: use MNE matlab toolbox instead of fiff access

SKIP_LENGTH = 100;

B = [];
nBuffs = [];
tlast = [];
skipevents = [];
sf = [];

try
  if ~isnumeric(channels), % Convert channel names
    if ischar(channels),
      channels={channels};
    end;
    cn=channames(filename);
    rcn=zeros(length(channels),1);
    for ii=1:length(channels),
      jj=min(find(strcmp(channels{ii},cn)));
      if isempty(jj),
        lasterr(['Unknown channel ' channels{ii}]);
        error(['Unknown channel ' channels{ii}]);
      end;
      rcn(ii)=jj;
    end;
    channels=rcn;
  end;

  rawdata('any',filename); % open raw file
  sf=rawdata('sf');
  status='';
  B=[];
  skipevents=[];
  s=0;
  nBuffs = 0;
  while ~strcmp(status,'eof')
    [M,status]=rawdata('next');
    if strcmp(status,'error'),
      lasterr('File error');
      error('File error');
    elseif strcmp(status,'ok'),
      M=M(channels,:);
      nSamples=size(M,2);
      nBuffs = nBuffs + 1;
      s=s+nSamples;
    elseif strcmp(status,'skip'),
      M=zeros(length(channels),SKIP_LENGTH);
      skipevents(end+1).type = 'skip';
      skipevents(end  ).latency = s;
      skipevents(end  ).condition = 0;
      skipevents(end  ).duration = SKIP_LENGTH;
      s=s+SKIP_LENGTH;
    end;
    B=[B M];
  end;
  tlast = rawdata('t')-0.1;
  rawdata close;
catch
  rawdata close;
  warning(lasterr);
end;

