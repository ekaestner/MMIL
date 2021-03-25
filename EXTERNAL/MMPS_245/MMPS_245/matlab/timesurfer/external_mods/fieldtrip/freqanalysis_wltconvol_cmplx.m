function [freq] = freqanalysis_wltconvol_cmplx(cfg,data,trials_flag)
% FREQANALYSIS_WLTCONVOL performs time-frequency analysis on any time series trial data
% using the 'wavelet method' based on Morlet wavelets.
%
% [freq] = wltanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from 
% the PREPROCESSING function. The configuration should be according to 
%
% cfg.output     = 'pow'; 
%
% If cfg.output is 'pow', then only the power spectra will be returned. 
% In that case, cfg should specify the names of channels to be analyzed in the  
% tag cfg.channel, which is a cell-array, which can contan individual channel labels
% or known groups of channels such as 'all', 'MEG', 'EEG'.
% See the CHANNELSELECTION function for more information.
% 
% Alternatively:
%
% cfg.output      = 'powandcsd';
%
% If cfg.output is 'powandcsd', then the power- and cross-spectra are returned. 
% In that case, cfg should contain the tag cfg.channelcmb, specifying the channel combinations:
%
% cfg.channelcmb      =[{'MZO02' 'MLP11'}; {'MZO02' 'MLP12'}; {'MZO02' 'MLP21'}];
%
% If you want to select more than just a few combinations of channels,
% e.g. all, you can use the CHANNELCOMBINATIONS helper function. Look
% at the help of that function for more options and to learn how to
% make specific channel combinations. Selecting all combinations is
% done with
%
% cfg.channelcmb   = channelcombination({'all' 'all}, data.label);
%
% The spectrum of the data is estimated at frequencies and time points
% defined by the voctors:
%
% cfg.foi (e.g. cfg.foi = 5:5:150; ) units: Hertz
% cfg.toi (e.g. cfg.toi = 0:.1:2; ) units: s 
%
% The temporal and spectral resolution of the analysis is determined
% by the 'width' of the wavelet which is set by cfg.width (default: 7).
%
% The standard deviation in the frequency domain (sf) at frequency f0 is 
% defined as: sf = f0/width
% The standard deviation in the temporal domain (st) at frequency f0 is 
% defined as: st = width/f0 = 1/(2*pi*sf)
%
% cfg.width can be specified as a constant for a 'classical constant-Q' 
% wavelet analysis or can be specified as a vector defining a variable width 
% for each frequency (e.g. cfg.width = linspace(5,10,100); in case of 100 
% frequencies).
%
% cfg.gwidth determines the length of the used wavelets in standard deviations
% of the implicit gaussian kernel and should be choosen >= 3;
%
% Thus, a possible configuration is:
%
% cfg.method         = 'wltconvol';
% cfg.foi            = 1:1:100;
% cfg.width          = 7;
%  or
% cfg.width          = linspace(5,10,length(cfg.foi));
% cfg.toi            = 0:0.1:2; 
% cfg.gwidth         = 3;
% 
% Furthermore, you can specify, whether you want to obtain the average
% spectra: 
% cfg.keeptrials  = 'no';
%
% or whether you want to keep the spectra of the individual trials:
% cfg.keeptrials  = 'yes';

% Copyright (c) 2003,2004 F.C. Donders Centre, Markus Siegel
%
% $Log: freqanalysis_wltconvol.m,v $
% Revision 1.12  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.11  2006/02/07 20:08:22  roboos
% changed all occurences of a dimord with chancmb (was previous sgncmb) into chan
%
% Revision 1.10  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.9  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.8  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.7  2005/01/19 08:42:42  jansch
% removed obsolete code for generating and checking channelcombinations
% cleaned up handling of sgnindx/sgncmbindx
% ensured that power is computed for all channels that are in channelcmb
%
% Revision 1.6  2005/01/18 15:11:39  roboos
% Cleaned up configuration for sgn/sgncmb, now exclusively using channel/channelcmb which is consistent with rest of fieldtrip and freqanalysis documentation. Also the output now only contains freq.label/labelcmb and not any more the old sgn/sgncmb.
%
% Revision 1.5  2004/12/20 14:52:56  roboos
% changed rounding off of timboi
%
% Revision 1.4  2004/11/18 16:51:58  roboos
% updated the help: pointed to CHANNELCOMBINATION for cfg.channelcmb
%
% Revision 1.3  2004/10/28 09:00:33  roboos
% fixed bug in caller function detection (for matlab 6.1 and 6.5)
%
% Revision 1.2  2004/10/28 07:21:46  roboos
% added check to ensure that the FREQANALYSIS wrapper is started instead of the
% FERQANALYSIS_xxx subfunctions
%
% Revision 1.1  2004/09/21 12:04:53  marsie
% the wltconvol method of wltanalysis.m has been moved to this seperate function
%
% Revision 1.2  2004/08/25 20:47:18  marsie
% fixed scaling problems
%
% Revision 1.1  2004/08/25 19:18:43  marsie
% initial release
%

% make sure we have double precision data
for k = 1:length(data.time),  data.time{k}  = double(data.time{k}); end;
for k = 1:length(data.trial), data.trial{k} = double(data.trial{k}); end
try data.dofvec = double(data.dofvec); end

% expand cfg.width to array if constant width
if ~isfield(cfg, 'width') || isempty(cfg.width)
    if isfield(cfg,'sf') && ~isempty(cfg.sf), 
        % if sf is provided, then width = foi/sf
        cfg.width = cfg.foi./cfg.sf; 
    elseif isfield(cfg,'st') && ~isempty(cfg.st)
        % if st is provided, then width = foi * st * 2 * pi
        cfg.width = cfg.foi.*cfg.st.*(2*pi);
    else
	cfg.width = cfg.foi;
    end
    % cfg.width = single(cfg.foi)';
    % cfg.width = ones(length(cfg.foi),1).*4;
    % cfg.width = single(cfg.foi.*4)'; % is there any reason for this??
else
    if numel(cfg.width ) == 1, cfg.width  = ones(1,length(cfg.foi)) * cfg.width;  end
end
if ~isfield(cfg, 'gwidth') || isempty(cfg.gwidth)
    % cfg.gwidth = ones(length(cfg.foi),1).*(0.49*(double(numsmp))*pi/data.fsample);
    cfg.gwidth = ones(length(cfg.foi),1).*pi;
else
    if numel(cfg.gwidth) == 1, cfg.gwidth = ones(1,length(cfg.foi)) * cfg.gwidth; end
end

% padding factor
if size(cfg.width,2) == 1, cfg.width = cfg.width'; end
try 
    wlength = (cfg.gwidth.*cfg.width)./(pi.*cfg.foi); % this is the wavelet length in seconds
catch
    keyboard
end
padfactor = max(wlength./2); % this is the required padding

% pad
for k = 1:length(data.trial)
    pad{k}   = zeros(size(data.trial{k},1),floor(size(data.trial{k},2)*padfactor));
    torig{k} = data.time{k}; tstep = torig{k}(2) - torig{k}(1);
    tnew1    = torig{k}(1)-tstep*size(pad{k},2):tstep:torig{k}(1)-tstep;
    tnew2    = torig{k}(length(torig{k}))+tstep:tstep:torig{k}(length(torig{k}))+tstep*size(pad{k},2);
    tnew     = [tnew1 torig{k} tnew2];
    origoff  = data.offset;
    try origdof  = data.dofvec; end
    origtoi  = cfg.toi;
    % assign
    data.trial{k}    = [pad{k} data.trial{k} pad{k}];
    data.time{k}     = tnew;
    data.numsamples(k)  = length(tnew);
    data.cfg.latency = [tnew(1) tnew(length(tnew))];
    data.offset(k)      = round(data.time{k}(1)*data.fsample);
    try data.dofvec      = ones(1,length(origdof)*2); end
    cfg.toi          = data.time{k};
end
    
% ensure that this function is started as a subfunction of the FREQANALYSIS wrapper
[s, i] = dbstack;
if length(s)>1
    [caller_path, caller_name, caller_ext] = fileparts(s(2).name);
else
    caller_path = '';
    caller_name = '';
    caller_ext  = '';
end
% evalin('caller', 'mfilename') does not work for Matlab 6.1 and 6.5
if ~strcmp(caller_name, 'freqanalysis')
    error(['you should call FREQANALYSIS, instead of the ' upper(mfilename) ' subfunction']);
end

% set all the defaults
if ~isfield(cfg, 'method'),        cfg.method        = 'wltconvol';  end
if ~isfield(cfg, 'keeptrials'),    cfg.keeptrials    = 'no';         end
if ~isfield(cfg, 'output'),        cfg.output        = 'powandcsd';  end
if ~isfield(cfg, 'pad'),           cfg.pad           = 'maxperlen';  end

% setting a flag (csdflg) that determines whether this routine outputs
% only power-spectra or power-spectra and cross-spectra?
if strcmp(cfg.output,'pow')
    csdflg = 0;
elseif strcmp(cfg.output,'powandcsd')
    csdflg = 1;
end

% determine the corresponding indices of all channels
sgnindx     = match_str(data.label, cfg.channel);
numsgn      = size(sgnindx,1);
if csdflg
    % determine the corresponding indices of all channel combinations
    for k=1:size(cfg.channelcmb,1)
        sgncmbindx(k,1) = strmatch(cfg.channelcmb(k,1), data.label, 'exact');
        sgncmbindx(k,2) = strmatch(cfg.channelcmb(k,2), data.label, 'exact');
    end
    numsgncmb   = size(sgncmbindx,1);
    sgnindx     = unique([sgnindx(:); sgncmbindx(:)]);
    numsgn      = length(sgnindx);
    cutdatindcmb = zeros(size(sgncmbindx));
    for sgnlop = 1:numsgn
        cutdatindcmb(find(sgncmbindx == sgnindx(sgnlop))) = sgnlop;
    end
end

% if rectan is 1 it means that trials are of equal lengths
numper = size(data.trial,2);
rectan = 1;
for perlop = 1:numper
    numdatbnsarr(perlop,1) = size(data.trial{perlop},2);
    if numdatbnsarr(perlop,1) ~= numdatbnsarr(1,1), rectan = 0; end
end

%if cfg.pad is 'maxperlen', this is realized here:
if ischar(cfg.pad) && strcmp(cfg.pad,'maxperlen'), cfg.pad = max(numdatbnsarr,[],1) ./ data.fsample; end
numsmp = int32(cfg.pad .* data.fsample);

% keeping trials and/or tapers?
if strcmp(cfg.keeptrials,'no' ), keep = 1; end
if strcmp(cfg.keeptrials,'yes'), keep = 2; end
      
% do the computation for WLTCONVOL
if strcmp(cfg.method,'wltconvol')
    numfoi = length(cfg.foi);
    if keep == 1, dimord = 'chan_freq_time';     end
    if keep == 2, dimord = 'rpt_chan_freq_time'; end
    if 0,
        % Min Xuang's method (very very slow)
        for perlop = 1:numper % trials
            for sgnlop = 1:numsgn % signals
                S = data.trial{perlop}(sgnindx(sgnlop),:);
                for foilop = 1:numfoi
                    width = ceil(foilop/2);
                    a     = traces2TF_complex(S,cfg.foi,data.fsample,width);
                    freq.spctrcmpx(perlop,sgnlop,foilop,:) =     a(foilop,1:size(a,2)/2);
                    freq.powspctrm(perlop,sgnlop,foilop,:) = abs(a(foilop,1:size(a,2)/2));
                end
            end
        end
    else
        minoffset    = min   (data.offset);
        timboi       = round (cfg.toi .* data.fsample - minoffset);
        toi          = round (cfg.toi .* data.fsample) ./ data.fsample;
        numtoi       = length(toi);
        knlspctrmstr = cell  (numfoi,1);
        % calculation of the kernel for each frequency of interest
        for foilop = 1:numfoi 
            dt                  = 1/data.fsample;                       % time differential
            sf                  = cfg.foi(foilop)/cfg.width(foilop);    % sigma_f, controls the frequency resolution
            st                  = 1/(2*pi*sf);                          % sigma_t, which controls the temporal resolution
            toi2                = -cfg.gwidth(foilop)*st:dt:cfg.gwidth(foilop)*st;
            A                   = 1/sqrt(st*sqrt(pi));                  % normalization constant for the wavelet
            tap                 = (A*exp(-toi2.^2/(2*st^2)))';          % wavelet transform
            acttapnumsmp        = size(tap,1);                          % size of the taper
            taplen(foilop)      = acttapnumsmp;                         % number of samples included in the taper
            ins                 = ceil(numsmp./2) - floor(acttapnumsmp./2);
            prezer              = zeros(ins,1);
            pstzer              = zeros(numsmp - ((ins-1) + acttapnumsmp)-1,1);%
            ind                 = (0:acttapnumsmp-1)'.*(2.*pi./data.fsample).*cfg.foi(foilop);
            knlspctrmstr{foilop}= complex(zeros(1,numsmp));%
            knlspctrmstr{foilop}= fft(complex(vertcat(prezer,tap.*cos(ind),pstzer),vertcat(prezer,tap.*sin(ind),pstzer)),[],1)';
        end
        % create output arrays
        if keep == 1
            powspctrm = zeros  (numsgn,numfoi,numtoi);
            spctrcmpx = complex(zeros(numsgn,numfoi,numtoi));
            if csdflg, crsspctrm = complex(zeros(numsgncmb,numfoi,numtoi)); end
            cntpertoi = zeros(numfoi,numtoi);
            dimord    = 'chan_freq_time';
        elseif keep == 2
            powspctrm = zeros (numper,numsgn,numfoi,numtoi);
            spctrcmpx = complex(zeros(numper,numsgn,numfoi,numtoi));
            if csdflg, crsspctrm = complex(zeros(numper,numsgncmb,numfoi,numtoi)); end
            dimord    = 'rpt_chan_freq_time';
        end
        for perlop = 1:numper
            if keep == 2, cnt = perlop; end
            numdatbns = numdatbnsarr(perlop,1);
            datspctra = complex(zeros(numsgn,numsmp));
            prepad    = zeros(1,(data.offset(perlop) - minoffset));
            pstpad    = zeros(1,minoffset + numsmp - (data.offset(perlop) + numdatbns));
            sz        = size(data.trial{perlop});
            for sgnlop = 1:numsgn
                datspctra(sgnlop,:) = fft([prepad,data.trial{perlop}(sgnindx(sgnlop),:),pstpad],[],2);
            end
            for foilop = 1:numfoi
                actfoinumsmp    = taplen(foilop); % taper length
                acttimboiind    = find(timboi >= (-minoffset + data.offset(perlop) +             (actfoinumsmp ./ 2)) & ...
                                       timboi <  (-minoffset + data.offset(perlop) + numdatbns - (actfoinumsmp ./2)));
                nonacttimboiind = find(timboi <  (-minoffset + data.offset(perlop) +             (actfoinumsmp ./ 2)) | ...
                                       timboi >= (-minoffset + data.offset(perlop) + numdatbns - (actfoinumsmp ./2)));
                acttimboi       = timboi(acttimboiind);
                numacttimboi    = length(acttimboi);
                if keep == 1, cntpertoi(foilop,acttimboiind) = cntpertoi(foilop,acttimboiind) + 1; end
                if keep == 3, cnt = 1; end
                autspctrmacttap = complex(zeros(numsgn,numacttimboi));
                if numacttimboi > 0
                    for sgnlop = 1:numsgn
                        dum = fftshift(ifft(datspctra(sgnlop,:).*knlspctrmstr{foilop},[],2));
                        autspctrmacttap(sgnlop,:) = dum(acttimboi);
                    end
                    warning off MATLAB:divideByZero;
                    powdum =  2.* (abs(autspctrmacttap).^2) ./ data.fsample;
                    powcpx =  2.*     autspctrmacttap ./data.fsample;
                    warning on MATLAB:divideByZero;
                end
                ndx2 = numacttimboi;
                vect = acttimboiind;
                if keep == 1 && numacttimboi >0
                    powspctrm(:,foilop,vect)     = powspctrm(:,foilop,vect) + reshape(powdum(:,1:ndx2),[numsgn,1,ndx2]);
                    spctrcmpx(:,foilop,vect)     = spctrcmpx(:,foilop,vect) + reshape(powcpx(:,1:ndx2),[numsgn,1,ndx2]);
                elseif keep == 2 && numacttimboi >0
                    powspctrm(cnt,:,foilop,vect) = powspctrm(cnt,:,foilop,vect) + reshape(powdum(:,1:ndx2),[1,numsgn,1,ndx2]);
                    spctrcmpx(cnt,:,foilop,vect) = spctrcmpx(cnt,:,foilop,vect) + reshape(powcpx(:,1:ndx2),[1,numsgn,1,ndx2]);
                    powspctrm(cnt,:,foilop,nonacttimboiind) = nan;
                    spctrcmpx(cnt,:,foilop,nonacttimboiind) = nan;
                elseif keep == 2 && numacttimboi == 0
                    powspctrm(cnt,:,foilop,nonacttimboiind) = nan;
                    spctrcmpx(cnt,:,foilop,nonacttimboiind) = nan;
                end
                if csdflg
                    csddum = 2.* (autspctrmacttap(cutdatindcmb(:,1),:) .*conj(autspctrmacttap(cutdatindcmb(:,2),:)))./data.fsample;
                    if keep == 1 && numacttimboi > 0
                        crsspctrm(:,foilop,acttimboiind) = crsspctrm(:,foilop,acttimboiind) + reshape(csddum,[numsgncmb,1,numacttimboi]);
                    elseif keep == 2 && numacttimboi > 0
                        crsspctrm(cnt,:,foilop,acttimboiind) = crsspctrm(cnt,:,foilop,acttimboiind) + reshape(csddum,[1,numsgncmb,1,numacttimboi]);
                        crsspctrm(cnt,:,foilop,nonacttimboiind) = nan;
                    elseif keep == 2 && numacttimboi == 0
                        crsspctrm(cnt,:,foilop,nonacttimboiind) = nan;
                    end
                end
            end % of foilop
        end % of perlop
        if trials_flag
            freq.powspctrm  = single(powspctrm);
            freq.spctrcmpx  = single(spctrcmpx);
        elseif keep == 1
            warning off matlab:dividebyzero;
            freq.powspctrm(:,:,:) = single(powspctrm(:,:,:) ./ repmat(permute(cntpertoi,[3,1,2]),[numsgn,1,1]));
            freq.spctrcmpx(:,:,:) = single(spctrcmpx(:,:,:) ./ repmat(permute(cntpertoi,[3,1,2]),[numsgn,1,1]));
            if csdflg, freq.crsspctrm(:,:,:) = single(crsspctrm(:,:,:) ./ repmat(permute(cntpertoi,[3,1,2]),[numsgncmb,1,1])); end
            warning on matlab:dividebyzero;
        end
    end
end

% collect the results
freq.label      = data.label(sgnindx);
freq.dimord     = dimord;
freq.freq       = cfg.foi;
freq.time       = cfg.toi;

if csdflg
    freq.labelcmb   = cfg.channelcmb;
    freq.crsspctrm  = crsspctrm;
end

try freq.grad = data.grad; end   % remember the gradiometer array
try freq.elec = data.elec; end   % remember the electrode array

% add information about the version of this function to the configuration
try
    % get the full name of the function
    cfg.version.name = mfilename('fullpath');
catch
    % required for compatibility with Matlab versions prior to release 13 (6.5)
    [st, i1] = dbstack;
    cfg.version.name = st(i1);
end
cfg.version.id = '$Id: freqanalysis_wltconvol.m,v 1.12 2006/02/23 10:28:16 roboos Exp $';
% remember the configuration details of the input data
try cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

% go back (remove padding)
if trials_flag  
  freq.powspctrm = freq.powspctrm(:,:,:,length(tnew1)+1:length(tnew1)+length(torig{1}));
  freq.spctrcmpx = freq.spctrcmpx(:,:,:,length(tnew1)+1:length(tnew1)+length(torig{1}));
  if csdflg, freq.crsspctrm = freq.crsspctrm(:,:,:,length(tnew1)+1:length(tnew1)+length(torig{1})); end
  freq.time      = torig{1};
  freq.cfg.toi   = torig{1};  
else
  freq.powspctrm = freq.powspctrm(:,:,length(tnew1)+1:length(tnew1)+length(torig{1}));
  freq.spctrcmpx = freq.spctrcmpx(:,:,length(tnew1)+1:length(tnew1)+length(torig{1}));
  if csdflg, freq.crsspctrm = freq.crsspctrm(:,:,length(tnew1)+1:length(tnew1)+length(torig{1})); end
  freq.time      = torig{1};
  freq.cfg.toi   = torig{1};
end
return


