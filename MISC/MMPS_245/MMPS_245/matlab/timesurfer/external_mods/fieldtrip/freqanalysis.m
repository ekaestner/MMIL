function [freq] = freqanalysis(cfg, data, trials_flag)
% FREQANALYSIS performs frequency and time-frequency analysis
% on time series data over multiple trials
%
% Use as
%   [freq] = freqanalysis(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the PREPROCESSING function. The configuration depends on the type
% of computation that you want to perform. If you want to compute only
% power-spectra, you should specify
%
%   cfg.output     = 'pow'; 
%   cfg.channel    = cell-array with selection of channels, default 'all'
%                    see CHANNELSELECTION for details
%
% If you want to compute both power-spectra and cross-spectra, you should
% specify
%
%   cfg.output     = 'powandcsd';
%   cfg.channelcmb = cell-array with selection of combinations between
%                    channels, see CHANNELCOMBINATION for details
%
% There are different methods of calculating the spectra, which are handled
% by separate subfunctions which are called by freqanalysis. The method
% is specified with one of the following
%
%   cfg.method = 'mtmfft'    (or alternatively cfg.method = 'fft')
%   cfg.method = 'mtmconvol' (or alternatively cfg.method = 'convol')
%   cfg.method = 'wltconvol' 
%   cfg.method = 'tfr' 
%
% Use MTMFFT if you want to analyse an entire spectrum for the entire data
% length. Use MTMCONVOL, WLTCONVOL or TFR if you want to analyse one or more
% frequencies for multiple analysis windows.  You should read the help of
% the respective subfunction for the corresponding parameter options and for
% a detailled explanaiton of each method. A short summary is given below.
%
% Method 'mtmfft' is handled by the function FREQANALYSIS_MTMFFT. This 
% implements multitaper frequency transformation. It has the additional parameters
%   cfg.foilim, cfg.tapsmofrq, cfg.pad, cfg.keeptrials. cfg.keeptapers,
%   cfg.pad
%
% Method 'mtmconvol' is handled by the function FREQANALYSIS_MTMCONVOL. This 
% implements multitaper time-frequency transformation based on multiplication in
% the frequency domain. It has the additional parameters
%   cfg.foi, cfg.toi, cfg.t_ftimwin, cfg.tapsmofrq, cfg.pad, 
%   cfg.keeptrials. cfg.keeptapers, cfg.pad
%
% Method 'wltconvol' is handled by the function FREQANALYSIS_WLTCONVOL. This 
% implements wavelet time frequency transformation (using Morlet wavelets) based on
% multiplication in the frequency domain. It has the additional parameters
%   cfg.foi, cfg.toi, cfg.width, cfg.gwidth, cfg.keeptrials 
%
% Method 'tfr' is handled by the function FREQANALYSIS_TFR. This implements
% wavelet time frequency transformation (using Morlet wavelets) based on
% convolution in the time domain. It has the additional parameters
%   cfg.foi, cfg.waveletwidth, cfg.downsample 
%
% See also FREQANALYSIS_MTMFFT, FREQANALYSIS_MTMCONVOL, 
%          FREQANALYSIS_WLTCONVOL, FREQANALYSIS_TFR

% Copyright (C) 2003, F.C. Donders Centre, Pascal Fries
% Copyright (C) 2004, F.C. Donders Centre, Markus Siegel
%
% $Log: freqanalysis.m,v $
% Revision 1.27  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.26  2005/08/01 12:41:22  roboos
% added some pragma-like comments that help the matlab compiler with the dependencies
%
% Revision 1.25  2005/07/03 20:15:30  roboos
% replaced eval() by feval() for better support by Matlab compiler
%
% Revision 1.24  2005/06/02 12:24:47  roboos
% replaced the local conversion of average into raw trials by a call to the new helper function DATA2RAW
%
% Revision 1.23  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.22  2005/05/04 07:31:55  roboos
% remove avg after converting to single trial
%
% Revision 1.21  2005/03/08 10:39:50  roboos
% changed ordering: first check for backward compatibility, then set the (new) default
%
% Revision 1.20  2005/01/19 08:39:51  jansch
% added defaults for cfg.channel and cfg.channelcmb, delete cfg.channelcmb if output is not powandcsd
%
% Revision 1.19  2005/01/19 08:05:27  jansch
% fixed small bug in channelcombinations
%
% Revision 1.18  2005/01/18 15:20:38  roboos
% this function now already takes care of channelselection and, if needed, of selection of channelcombinations
%
% Revision 1.17  2005/01/18 15:05:15  roboos
% cleaned up configuration for sgn/label, now consistently using cfg.channel and cfg.channelcmb
%
% Revision 1.16  2004/11/01 11:34:53  roboos
% added support for timelocked trials, i.e. the result of timelockanalysis with keeptrials=yes
% tthese are now converted to raw trials prior to calling the freqanalysis_xxx subfunction
%
% Revision 1.15  2004/10/01 10:23:44  roboos
% fixed error in help for mtmconvol: f_timwin should be t_ftimwin
%
% Revision 1.14  2004/09/28 14:26:43  roboos
% fixed a typo in the help, added a line of comments
%
% Revision 1.13  2004/09/22 10:20:27  roboos
% converted to use external subfunctions time2offset and offset2time
% and add offset field to data structure if it is missing
%
% Revision 1.12  2004/09/21 12:07:28  marsie
% this version of the wrapper implements the new "freqanalysis_METHOD" convention
% for subfunction naming.
%
% Revision 1.11  2004/09/02 11:59:06  roboos
% restructured the help, fixed the copyrights
%
% Revision 1.10  2004/09/01 13:33:03  marsie
% switch to a wrapper function that calls the functions
% multitaperanalysis.m, wltanalysis.m and waveletanalysis.m for the corresponding
% methods
%
% LAST MOD: SRD 2012/09/11 Fix memory waste, or most of it.  Add memory use  $$$$$$$$$$$$$$$$$$$$$$$$$
%                      monitor options (but commented out)

% aid the Matlab Compiler to find some subfunctions that are called by feval
%#function freqanalysis_mtmconvol
%#function freqanalysis_mtmwelch
%#function freqanalysis_wltconvol
%#function freqanalysis_mtmfft
%#function freqanalysis_tfr

% make sure we have double precision, if not convert to it
data.avg   = double(data.avg); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ToDo: Do this at a higher level to save memory
data.var   = double(data.var);
data.trial = double(data.trial);

% for backward compatibility with old data structures
data = fixdimord(data);

% for backward compatibility: the old configuration used sgn/sgncomb or
% label/labelcmb, the new configuration uses channel/channelcmb
if (isfield(cfg, 'sgn') || isfield(cfg, 'label')) && ~isfield(cfg, 'channel')
  % clean up the configuration
  try, cfg.channel = cfg.label;     end
  try, cfg = rmfield(cfg, 'label'); end
  try, cfg.channel = cfg.sgn;       end
  try, cfg = rmfield(cfg, 'sgn');   end
end
if (isfield(cfg, 'sgncmb') || isfield(cfg, 'labelcmb')) && ~isfield(cfg, 'channelcmb')
  % clean up the configuration
  try, cfg.channelcmb = cfg.labelcmb;  end
  try, cfg = rmfield(cfg, 'labelcmb'); end
  try, cfg.channelcmb = cfg.sgncmb;    end
  try, cfg = rmfield(cfg, 'sgncmb');   end
end

% set the defaults for the channel selection
if ~isfield(cfg, 'channelcmb'),    cfg.channelcmb = 'all'; end
if ~isfield(cfg, 'channel'),       cfg.channel    = 'all'; end

% support 'fft' and 'convol' as methods for backward compatibility
if strcmpi(cfg.method,'fft'),    cfg.method = 'mtmfft';    end
if strcmpi(cfg.method,'convol'), cfg.method = 'mtmconvol'; end

% no cross-spectrum needs to be computed, hence remove the combinations from cfg
if isfield(cfg, 'output') && ~strcmp(cfg.output, 'powandcsd'), cfg = rmfield(cfg, 'channelcmb'); end
if isfield(cfg, 'output') &&  strcmp(cfg.output, 'powandcsd'), csdflag = 1; else csdflag = 0; end
  
% ensure that channelselection and selection of channelcombinations is 
% perfomed consistently, prior to calling freqanalysis_XXX
% cfg.channel = {'G1'};
cfg.channel = channelselection(cfg.channel, data.label);
if isfield(cfg, 'channelcmb'), cfg.channelcmb = channelcombination(cfg.channelcmb, data.label); end

% convert average to raw data using simple helper function
% this allows the same function to work on multiple data types
[data, inputdimord] = data2raw(data);

% orig            = data; %%%% DOUBLES epoch data, a waste, try this instead.  SRD 12.6.14
orig.numsamples = data.numsamples; %%%%%
orig.time =  data.time; %%%%%

for trial = 1:length(data.trial)         %%%%% INCREMENTALLY GROWS DATA STRUCTURE a trial at a time, memory waste, try preallocaion of memory for 'data'. XXX
    endd              = data.time{trial}(length(data.time{trial})); % Last time in trial epoch
    begin             = data.time{trial}(2);                        % First time
    stp               = (endd-begin)/length(data.time{trial});
    sz                = [1 length(cfg.channel) length(cfg.foi) length(cfg.toi)];
    % pad data                                                      % Expandes out the time for padding.
    data.numsamples   = data.numsamples*2;
    pre               = fliplr(data.trial{trial}(:,1:ceil(size(data.trial{trial},2)/2)));
    post              = fliplr(data.trial{trial}(:,  ceil(size(data.trial{trial},2)/2+1:size(data.trial{trial},2))));
    data.trial{trial} = [pre data.trial{trial} post];
    data.trial{trial} = data.trial{trial}(:,1:orig.numsamples*2);
    pret              = begin-size(pre,2)*stp:stp:begin-stp;
    postt             = endd+stp:stp:endd+size(post,2)*stp;
    data.time{trial}  = [pret data.time{trial} postt];
end

% the older freqnalayis_xxx functions work with the offset instead of time
% if the offset in samples is not in the data, compute it from the time axes
if ~isfield(data, 'offset')
    for i=1:length(data.time), data.offset(i) = time2offset(data.time{i}, data.fsample); end
end
cfg_c = cfg; orig_foi = cfg.foi; ndx = 1; clear freq freq_;  % dta = data;  %%%%% Another big memory waster 'dta' - doubles data, and isn't used.   SRD 12.6.14 
% if length(cfg.foi) > 1, 
%     limit = length(cfg.foi); unit = 2;
% else
%     limit = 1; unit = 1;
% end
% for k = 1:unit:limit
%     try cfg_c.foi    = cfg.foi   (k:k+1); catch cfg_c.foi    = cfg.foi   (k):cfg.foi   (k)+1; end
%     try cfg_c.sf     = cfg.sf    (k:k+1); catch cfg_c.sf     = cfg.sf    (k):cfg.sf    (k)+1; end
%     try cfg_c.gwidth = cfg.gwidth(k:k+1); catch cfg_c.gwidth = cfg.gwidth(k):cfg.gwidth(k)+1; end
%     data  = dta;
%     freq_ = feval(sprintf('freqanalysis_wltconvol_cmplx'),cfg_c,data,trials_flag); ndx = ndx + 1;
%     if k == 1, 
%         freq = freq_;
%     else
%         if cfg_c.foi(2) > cfg.foi(length(cfg.foi))
%             freq.powspctrm(:,k,:) = freq_.powspctrm(:,1,:);
%             freq.spctrcmpx(:,k,:) = freq_.spctrcmpx(:,1,:);
%         else
%             freq.powspctrm(:,k:k+1,:) = freq_.powspctrm;
%             freq.spctrcmpx(:,k:k+1,:) = freq_.spctrcmpx;
%         end
%     end
% end
for k = 1:length(cfg.foi) %%%%% INCREMENTALLY BUILDS FREQ DATA one freq at a time, memory waste, try preallocaion of memory for freq
    try cfg_c.foi    = cfg.foi   (k); end 
    try cfg_c.sf     = cfg.sf    (k); end
    try cfg_c.st     = cfg.st    (k); end
    try cfg_c.gwidth = cfg.gwidth(k); end
    try cfg_c.width  = cfg.width (k); end
%     data         = dta; %%%%% REMOVED, SRD 
    freq_        = feval(sprintf('freqanalysis_wltconvol_cmplx'),cfg_c,data,trials_flag);
    if k == 1 %%%%% INITIALIZE ONCE - don't grow it incrementally
%%%%% Dimensions of freq are trials, channels, foi, & time.  
% % %         freq = freq_; 
        freq.powspctrm = zeros(length(data.trial),length(cfg.channel),length(cfg.foi),length(freq_.time),'single'); %%%%%
        freq.spctrcmpx = zeros(length(data.trial),length(cfg.channel),length(cfg.foi),length(freq_.time),'single'); %%%%%
        freq.powspctrm(:,:,1,:) = freq_.powspctrm(:,:,1,:);  %%%%% 
        freq.spctrcmpx(:,:,1,:) = freq_.spctrcmpx(:,:,1,:);  %%%%%  
        freq.label = freq_.label; %%%%%
        freq.dimord = freq_.dimord; %%%%%
        freq.freq = freq_.freq; %%%%%
        freq.time = freq_.time; %%%%%
        freq.cfg = freq_.cfg; %%%%%
    else
      if trials_flag
        freq.powspctrm(:,:,k,:) = freq_.powspctrm(:,:,1,:);
        freq.spctrcmpx(:,:,k,:) = freq_.spctrcmpx(:,:,1,:);
        if csdflag, freq.crsspctrm(:,:,k,:) = freq_.crsspctrm(:,:,1,:); end
      else
%%%%% Dimensions of freq are channels, foi, & time.
        freq.powspctrm(:,k,:) = freq_.powspctrm(:,1,:);
        freq.spctrcmpx(:,k,:) = freq_.spctrcmpx(:,1,:);
        if csdflag, freq.crsspctrm(:,k,:) = freq_.crsspctrm(:,1,:); end
      end
    end
% memuse = whos('freq'); message = sprintf('Memory used for freq reached %3.3f Gig.',memuse.bytes/1e9); display(message);   %%%%%
% whos %%%%%

end

% data %%%%%


freq.cfg  = cfg;
freq.freq = orig_foi;

v1 = find(freq.time == orig.time{1}(1));
v2 = find(freq.time == orig.time{1}(length(orig.time{1})));

freq.time = orig.time{1};
if trials_flag
  freq.powspctrm = freq.powspctrm(:,:,:,v1:v2);
  freq.spctrcmpx = freq.spctrcmpx(:,:,:,v1:v2);
  if csdflag, freq.crsspctrm = freq.crsspctrm(:,:,:,v1:v2); end
else
  freq.powspctrm = freq.powspctrm(:,:,v1:v2);
  freq.spctrcmpx = freq.spctrcmpx(:,:,v1:v2);  
  if csdflag, freq.crsspctrm = freq.crsspctrm(:,:,v1:v2); end
end
% whos ('freq','data','dta') %%%%%

return