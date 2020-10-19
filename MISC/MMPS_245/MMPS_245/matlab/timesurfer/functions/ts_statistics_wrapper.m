function [stat_data] = ts_statistics_wrapper(data,varargin)
% Purpose: perform statistical analysis using FieldTrip functions, timelockstatistics and freqstatistics
%
% Usage: stat_data = ts_statistics( data, 'key1', val1, 'key2', val2, ...);
%   
% Example: stat_data = ts_statistics_wrapper(epoch_data,...
%                                   'events'      ,comparisons,...
%                                   'neighbours'  ,'yes',...
%                                   'neighbthresh',30,...
%                                   'verbose'     ,show_stdout_flag);
%
% Inputs:
%    data: epoch_data or avg_data, structure must contain the conditions tocompare
%    events: [n1 n2] will compare event codes n1 and n2.  Also, [n1 n2] [n1n3] ...}
%    cfg:  a structure specifying the configuration for MC clusterin
%    subj: a string specifying the subject name
%    hemi: a string specifying the hemisphere, i.e. 'lh' or 'rh'
%
% Outputs:
%   stat_data: a statistics data structure containing the following fields:
%         num_sensors: number of channels
%         sensor_info: information on the channels
%         sfreq: sampling frequency
%         stats: a structure array with following:
%            condition: the condition number in relation to the original data (-1 for between conditions)
%            svent_code: the event code (-1 for between conditions)
%            prob: the result of the test p-value - a matrix [channels x time points]
%            mask: mask (0 for bad channels)
%            time: vector of time points 
%            parms:
%                 method: the trip method used
%                 statistic: the trip statistic used
%                 alpha: alpha value used
%            num_trials: number of trials
%            norm_bl: result of Kolmogorov-Smirnov test for baseline data (0 for between conditions)
%            posclusters: a structure array with statistical results for positive clusters (sorted by p-value)
%            negclusters: a structure array with statistical results for negative clusters (sorted by p-value)
%            posclusterslabelmat: matrix [channels x time points] indicating pos clusters to which (channel,time)-pairs belong
%            negclusterslabelmat: matrix [channels x time points] indicating neg clusters to which (channel,time)-pairs belong
% 
% Parameters: default : options :  description
%     conditions: [] :: condition number (not event code) used to index data inside structure
%     type: 'both' : 'between','within','both' :
%     method: 'montecarlo' : 'montecarlo','analytic','stats','glm' : statistical method for comparing conditions
%     statistic: 'indepsamplesT' : 'indepsamplesT','indepsamplesF','indepsamplesregrT',...
%                'indepsamplesZcoh','depsamplesT','depsamplesF' ,'depsamplesregrT','actvsblT','ttest','ttest2','paired-ttest',...
%                'anova1','kruskalwallis' : sample level statistic for thresholding prior to clustering
%     actvsbl_method: 'stats' : 'montecarlo','analytic','stats','glm':  
%     actvsbl_statistic: 'ttest' : 'ttest','ttest2','actvsblT':
%     correctm: 'cluster' : 'no','max','cluster','bonferoni','fdr' :
%     parameter: 'individual' :: 
%     alpha: 0.05 :: alpha level of statistical test
%     tail: 0 : -1, 1, 0 : tail of statistical test
%     ivar: 1 ::
%     uvar: [] ::
%     wvar: [] ::
%     feedback: 'textbar' : 'gui', 'text', 'textbar', 'no' :
%     clusterstatistic: 'maxsum' : 'maxsum','maxsize','wcm' : cluster level statistic (corrects for multiple comparisons)
%     clusterthreshold: 'parametric' : 'parametric','nonparametric' :
%     clusteralpha:  0.05 : ::
%     clustercrtival: [] ::
%     latency: 'all' :: times to be analyzed [begining end]
%     blcorrection: 'no' : 'no','absolute','relative','relchange','zscore','zscores','yes', : whether to do baseline correction
%     bl_window: [] ::  [begin end] in seconds 
%     bl_tile: 0 : 0,1 :: 
%     avgoverchan: 'no' : 'no','yes' : whether to average over channels
%     avgovertime: 'no' : 'no','yes' : whether to average over time
%     numrandomization: 500 :: the number of monte carlo iterations
%     minnbchan: [] :  minimum number of channels to define a spatial cluster
%     neighbours: [] : neighbourhoods
%     neighbthresh: 25 :: 
%     subjdir: '' :: getenv(''SUBJECTS_DIR'')
%     design: [] :: design matrix (see FieldTrip documentation for timelockstatistics)
%     lpfilter: 'no' : 'yes','no' lowpass filter
%     hpfilter: 'no' : 'yes','no' highpass filter
%     bpfilter: 'no' : 'yes','no' bandpass filter
%     bsfilter: 'no' : 'yes','no' bandstop filter
%     lnfilter: 'no' : 'yes','no' line noise removal using notch filter
%     dftfilter: 'no' : 'yes','no' line noise removal using discrete fourier transform
%     lpfreq: [] :  lowpass  frequency in Hz
%     hpfreq: [] :  highpass  frequency in Hz
%     bpfreq: [] :  bandpass frequency range, specified as low high in Hz
%     lnfreq: [] :  line noise frequency in Hz 
%     dftfreq: [] :  line noise frequencies for DFT filter 
%     medianfiltord: [] :  length of median filter
%     blc: 'no' : 'yes', 'no' : whether to do baseline correction
%     blcwindow: [] ::  in seconds [begin end]
%     removemean: 'no': 'yes','no' :
%     detrend: 'no' : 'yes','no' : linear detrending
%     fieldtrip_flag: 1 : 0,1 :
%     precision: 'double' : 'single','double' : whether to save data with single precision instead of double
%     datatype: 'epoch_data' ::
%     chantype: 'all' : 'all' 'mag' 'grad1' 'grad2' 'eeg' 'other' 'grad' 'meg':
%     removebadchans: 1 ::
%     verbose: 1 : 0,1 :
%     logfile: [] : filename for log file (stdout will be appended)
%     logfid: 1 : file FID for log file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Mofified:  SRD 3/14/12 - Test for attempted use of wrong statistic, and return.
%  Correction: SRD 5/3/12 - Parameter was not being referenced properly in last fix.

parms = mmil_args2parms(varargin,...
                      {'conditions' , [],[],...
                       'events'     , [],[],...
                       'type'       , 'both',{'between','within','both'},...
                       'method'     , 'montecarlo',{'montecarlo','analytic','stats','glm'},...
                       'statistic'  , 'indepsamplesT',{'indepsamplesT','indepsamplesF','indepsamplesregrT',...
                                          'indepsamplesZcoh','depsamplesT','depsamplesF' ,...
                                          'depsamplesregrT','actvsblT','ttest','ttest2','paired-ttest',...
                                          'anova1','kruskalwallis'},...
                       'actvsbl_method','stats',{'montecarlo','analytic','stats','glm'},...
                       'actvsbl_statistic','ttest',{'ttest','ttest2','actvsblT'},...                   
                       'correctm'   , 'cluster',{'no','max','cluster','bonferoni','fdr'},...
                       'parameter', 'individual',[],...
                       'alpha'      , 0.05, [],...
                       'tail'       ,    0, {-1, 1, 0},...
                       'ivar'       ,   1,[],...
                       'uvar'       ,   [],[],...
                       'wvar'       ,   [],[],...
                       'feedback'   , 'textbar',{'gui', 'text', 'textbar', 'no'},...
                       'clusterstatistic', 'maxsum', {'maxsum','maxsize','wcm'},...
                       'clusterthreshold', 'parametric', {'parametric','nonparametric'},...
                       'clusteralpha'    ,        0.05 , [],...
                       'clustercrtival'  ,        [],[],...
                       'clustertail',0,{-1, 1, 0},...
                       'latency'    , 'all',[],...
 											 'blcorrection','no',{'no','absolute','relative','relchange','zscore','zscores','yes',''},...
											 'bl_window',[],[],...		
                       'bl_tile',0,{0,1},...
                       'avgoverchan', 'no',[],...
                       'avgovertime', 'no',[],...
                       'numrandomization',500,[],...
                       'minnbchan',[],[],...
                       'neighbours',[],[],...
                       'neighbthresh',25,[],...
                       'subjdir','',[],...
                       'subj','',[],...
                       'hemi','',[],...
                       'design',[],[],...
                       'cfg',[],[],...
                       'lpfilter','no',{'yes','no'},...
                       'hpfilter','no',{'yes','no'},...
                       'bpfilter','no',{'yes','no'},...
                       'bsfilter','no',{'yes','no'},...
                       'lnfilter','no',{'yes','no'},...
                       'dftfilter','no',{'yes','no'},...
                       'medianfilter','no',{'yes','no'},...
                       'lpfreq',[],[],...
                       'hpfreq',[],[],...
                       'bpfreq',[],[],...
                       'lnfreq',[],[],...
                       'dftfreq',[],[],... 
                       'medianfiltord',[],[],...
                       'blc','no',{'yes', 'no'},...
                       'blcwindow',[],[],...
                       'removemean','no',{'yes','no'},...
                       'detrend','no',{'yes','no'},...      
                       'fieldtrip_flag',1,{0,1},...
                       'precision','double',{'double','single'},...
                       'datatype','epoch_data',[],...
                       'chantype','all',[],...
                       'removebadchans',1,[],...
                       'verbose',1,{0,1},...
                       'logfile',      [],[],...
                       'logfid',       [1],[], ...                          
                      },...
                      false);

if iscell(parms.events)
  data = ts_checkdata_header(data,'precision',parms.precision,'events',[parms.events{:}],'datatype',parms.datatype);
else
  data = ts_checkdata_header(data,'precision',parms.precision,'events',parms.events,'datatype',parms.datatype);
end
if any(find(ismember(varargin(1:2:end),'removebadchans')))
  rmix = find(ismember(varargin(1:2:end),'removebadchans'));
  rmix = rmix*2-1;
  rmix = [rmix rmix+1];
  varargin(rmix) = [];
end
data = ts_data_selection(data,varargin{:},'removebadchans',parms.removebadchans);
[datatype datafield] = ts_object_info(data);

% 3/14/12 SRD - Check for user attempting wrong test for the number of conditions.  5/3/12 SRD - Epicycle to change statistic below to parms.statistic
if ( (length(parms.events) > 2) && ~strcmp(parms.statistic,'indepsamplesF') && ~strcmp(parms.statistic,'depsamplesF') && ~strcmp(parms.statistic,'kruskalwallis') ) 
  fprintf('%s: Cannot do this type of test on more than 2 conditions, try an (in)depsamplesF or kruskalwallis statistic.\n',mfilename);
  return;
end  % SRD 3/14/12  SRD 5/3/12 

if ~parms.fieldtrip_flag
  data = ts_preproc(data,varargin{:});
  parms.detrend = 'no';
  parms.blc     = 'no';
  parms.lpfilter = 'no';
  parms.hpfilter = 'no';
  parms.bpfilter = 'no';
else
  if isequal(parms.detrend,1) , parms.detrend  = 'yes'; end
  if isequal(parms.blc,1)     , parms.blc      = 'yes'; end
  if isequal(parms.lpfilter,1), parms.lpfilter = 'yes'; end
  if isequal(parms.hpfilter,1), parms.hpfilter = 'yes'; end
  if isequal(parms.bpfilter,1), parms.bpfilter = 'yes'; end
end

args = mmil_parms2args(parms);
chans = 1:data.num_sensors;
stat_data = [];
if (ischar(parms.neighbours) && strcmp(parms.neighbours,'yes')) || ...
  isnumeric(parms.neighbours) && ~isempty(parms.neighbours) && parms.neighbours == 1
mmil_logstr(parms,'Performing Monte Carlo statistical comparison with spatiotemporal clustering\n');  
% fprintf('Performing Monte Carlo statistical comparison with spatiotemporal clustering\n');
% % get angles
  zshift = 0;
  thresh = parms.neighbthresh;
  sens   = data.sensor_info;
  if ~isempty(parms.subj)
    % source space
    parms = ts_vertex_neighbours (parms,parms.subjdir,parms.subj,parms.hemi);
  else
    % sensor space
    parms = ts_MEG_sensor_neighbours (parms,sens,thresh,zshift);
  end
  if isempty(parms.minnbchan), parms.minnbchan = 0; end  
  args = mmil_parms2args(parms);
  stat_data = ts_statistics(data,args{:});
else
  for k = 1:length(chans)
    ch = chans(k);
    if data.sensor_info(ch).badchan
      mmil_logstr(parms,'skipping bad channel: %s\n',data.sensor_info(ch).label);
%       fprintf('skipping bad channel: %s\n',data.sensor_info(ch).label);
      continue;
    end
    dat = ts_data_selection(data,'channels',ch,'removebadchans',1,'verbose',0);
    if ~isstruct(data)
      mmil_logstr(parms,'skipping empty channel: %s\n',data.sensor_info(ch).label);
%       fprintf('skipping empty channel: %s\n',data.sensor_info(ch).label);
      continue;
    end
    mmil_logstr(parms,'%s: running stats on channel %g of %g\n',mfilename,k,length(chans)); 
%     fprintf('%s: running stats on channel %g of %g\n',mfilename,k,length(chans));
    stats = ts_statistics(dat,args{:});
    if isempty(stat_data)
      n  = data.num_sensors;
      nt = length(stats.stats(1).time);
      nc = length(stats.stats);
      stat_data                           = rmfield(stats,{'sensor_info','stats'});
      stat_data.sensor_info               = init_sensor_info(n);
      stat_data.num_sensors               = length(chans);
      [stat_data.stats(1:nc).mask]        = deal(false(n,nt));
      [stat_data.stats(1:nc).prob]        = deal(zeros(n,nt));
      [stat_data.stats(1:nc).stat]        = deal(zeros(n,nt));
      [stat_data.stats(1:nc).dimord]      = deal(stats.stats(1).dimord);
      [stat_data.stats(1:nc).posclusters] = deal([]);
      [stat_data.stats(1:nc).negclusters] = deal([]);
      [stat_data.stats(1:nc).posclusterslabelmat] = deal(zeros(n,nt));
      [stat_data.stats(1:nc).negclusterslabelmat] = deal(zeros(n,nt));          
      for j = 1:nc
        stat_data.stats(j).time         = stats.stats(j).time;
        try [stat_data.stats(j).label{1:n}] = deal(''); end
        try 
          stat_data.stats(j).parms      = stats.stats(j).parms;
        catch
          stat_data.stats(j).parms      = parms;
        end
        stat_data.stats(j).event_code   = stats.stats(j).event_code;
      end
    end
    try
      stat_data.sensor_info(ch)  = stats.sensor_info;
    catch
      stats.sensor_info = orderfields(stats.sensor_info,stat_data.sensor_info);
      stat_data.sensor_info(ch)  = stats.sensor_info;
    end
    for j = 1:length(stats.stats)
      try stat_data.stats(j).mask(ch,:) = stats.stats(j).mask;  end
      try stat_data.stats(j).prob(ch,:) = stats.stats(j).prob;  end
      try stat_data.stats(j).stat(ch,:) = stats.stats(j).stat;  end
      try stat_data.stats(j).label(ch)  = stats.stats(j).label; end
      try stat_data.stats(j).posclusters = cat(2,stat_data.stats(j).posclusters,stats.stats(j).posclusters); end
      try stat_data.stats(j).negclusters = cat(2,stat_data.stats(j).negclusters,stats.stats(j).negclusters); end    
      try stat_data.stats(j).posclusterslabelmat(ch,:) = stats.stats(j).posclusterslabelmat; end
      try stat_data.stats(j).negclusterslabelmat(ch,:) = stats.stats(j).negclusterslabelmat; end    
    end
    clear dat stats;
  end  
end

if ~isfield(stat_data,'parms')
  stat_data.parms = parms;
end
% if ~strcmp(parms.type,'within')
%   stat_data.stats(1).event_code = -1;
% end

function sens = init_sensor_info(nchan)
[sens(1:nchan).label]      = deal('');
[sens(1:nchan).typestring] = deal('');
[sens(1:nchan).type]       = deal(1);
[sens(1:nchan).kind]       = deal(2);
[sens(1:nchan).badchan]    = deal(1);
[sens(1:nchan).lognum]     = deal([]);
[sens(1:nchan).loc]        = deal([]);  
