function [data reject_data] = ts_reject(data,varargin)
% Purpose: run automatic threshold-based rejection, manual visual
%          rejection, Manual ICA, Automatic ICA, and/or manual channel rejection
%
% Usage: [data reject_data] = ts_reject(data,'key1',val1, ... )
%
% Example:  [epoch_data,reject_data] = ts_reject(epoch_data, ...
%                                       'reject_auto_flag'        ,1,...
%                                       'reject_manual_flag'      ,0,...
%                                       'ICA_auto_flag'           ,0,...
%                                       'ICA_manual_flag'         ,0,...
%                                       'keepbadtrials_flag'      ,0,...
%                                       'write_reject_log'        ,1,...
%                                       'write_summary_reject_log',1,...
%                                       'write_events_file'       ,1,...
%                                       'reject_grad'             ,1000,...
%                                       'reject_eeg'              ,150,...
%                                       'reject_eog'              ,150,...
%                                       'rootoutdir'              ,'/home/user/rootoutdir/',...
%                                       'prefix'                  ,'prefix',...
%                                       'verbose'                 ,1);
%
% Imputs:
%     data: timesurfer data struct of epoch, everage, or time-frequency data
%
% Outputs: 
%     data: timesurfer data struct of epoch, everage, or time-frequency data 
%           with rejected trials removed (if keepbadtrials_flag is set to 0)
%       -OR-
%     data: timesurfer data struct of epoch, everage, or time-frequency data 
%           with rejected trials marked in the trial_info sectuion of the data struct
%           (if keepbadtrials_flag is set to 1)
%     reject_data: timesurfer data struct of epoch, everage, or
%                  time-frequency data containing only the rejected trials
%
% Parameters: default : options :  description 
%     keepbadtrials_flag: 1 : 0,1 : whether to mark trials as bad in the
%                                    trial info section of the data struct (1) 
%                                    or remove the trials alltogether (0).
%     reject_auto_flag: 0 : 0,1 : whether to perform automatic threshold-based rejection
%     reject_manual_flag: 0 : 0,1 : whether to perform manual visual rejection
%     ICA_auto_flag: 0 : 0,1 : whether to perform automatic ICA
%     ICA_manual_flag: 0 : 0,1 : whether to perform manual ICA
%     calc_ncov_flag: 0 : 0,1 : whether to calculate/recalculate noise
%                               covariance matrix after epoching/re-epoching
%     write_events_file: 0 : 0,1 : whether to generate a .txt file
%                                  containing events information
%     ICA_ref_chan: 'EOG061' ::
%     ICA_chantype: 'all' : 'all','mag','grad1','grad2' 'eeg','other','grad','meg' :
%     ICA_maxsteps: 20 ::
%     ICA_ntrial: 5 ::
%     ICA_start_trial: 1 ::
%     ICA_ncomponents: 80 ::
%     ICA_plotype: 'activations', 'activations','alltrials'
%     ICA_rescale_flag: 1 : 0,1 :
%     ICA_sorttrials: 0 : 0,1 :
%     ICA_saveout_flag: 0 : 0,1 :
%     reject_method: 'summary' ::
%     reject_metric: 'var' ::
%     reject_mag: 6000 :: automatic rejection threshold for magnetometer channels (fT)
%                         if 0, rejection based on magnetometers is disabled
%     reject_grad: 3000 :: automatic rejection threshold for gradiometer channels (fT/cm)
%                          if 0, rejection based on gradiometers is disabled
%     reject_eeg: 0 :: automatic rejection threshold for eeg channels (uV)
%                      if 0, rejection based on eeg is disabled
%     reject_eog: 200 :: automatic rejection threshold for eog channel (uV)
%                        if 0, rejection based on eog is disabled
%     reject_ieeg: [] :: automatic rejection threshold for ieeg channels (uV)
%                        if 0, rejection based on eeg is disabled
%     ignorebadchans_flag: 0 : 0,1 : ignore channels marked in the sensorinfo as bad when applying automatic rejection thresholds
%     prescale_mag: 10^15 :: scale to fT
%     prescale_grad: 10^13 :: scale to fT/cm
%     prescale_eeg: 10^6 :: scale to uV
%     prescale_eog: 10^6 :: scale to uV
%     prefix: [] :: specifies a prefix to use for all output file names
%     rootoutdir: pwd :: All output files will be saved in the rootoutdir directory (default is working directory)
%     verbose: 1 : 1,0 : whether to produce fulll output
%     logfile: [] :: filename for log file (stdout will be appended)
%     logfid: 1 :: file FID for log file
%     rejectfile: [] ::
%     filename: [] ::
%     events: [] :: array of event codes to process
%     badchanfile: [] :: full file path to .txt file containing a line delimited 
%                        list of bad channel labels
%     badchans: [] :: vector of channel indices for rejected channels
%     toilim: [] :: limits on the times of interest, [begin end] in seconds
%     write_reject_only_log: 0 : 0,1 : whether to create a .txt logfile containing
%                                      only rejected trial and channel information
%                                      in Rajan's format
%     write_reject_log: 0 : 0,1 : whether to create a .log logfile containing 
%                                 events (rejected and non-rejected) information 
%                                 in Sanja's format
%     write_summary_reject_log: 0 : 0,1 : whether to create .log logfile
%                                         of rejections summarized by event code 
%                                         in Burke's format

% added ICA_start_trial 09/01/11 -BQR




% data  = ts_checkdata_header(data);
parms = mmil_args2parms(varargin,...
             {'keepbadtrials_flag'     ,1,{0,1},...
             'calc_ncov_flag',false,[false true],...
             'ICA_auto_flag'           ,0,{0,1},...
             'ICA_manual_flag'         ,0,{0,1},...
             'reject_auto_flag'        ,0,{0,1},...
             'reject_manual_flag'      ,0,{0,1},...
             'ICA_ref_chan'            ,'EOG061',[],...
             'ICA_chantype'            ,'all',[],...
             'ICA_maxsteps'            ,20,[],...
             'ICA_ntrial'              ,5,[],...
             'ICA_start_trial'         ,1,[],...
             'ICA_ncomponents'         ,80,[],...
             'ICA_plottype'            ,'activations',{'activations','alltrials'},...
             'ICA_rescale_flag'        ,1,{0,1},...
             'ICA_sorttrials'          ,0,{0,1},... 
             'ICA_saveout_flag'        ,0,[]...
             'reject_method'           ,'summary',[],...
             'reject_metric'           ,'var',[],...     
             'rejectfile'              ,[],[],...
             'reject_data'             ,[],[],...
             'reject_mag'              ,6000,[],...
             'reject_grad'             ,3000,[],...
             'reject_eeg'              ,0,[],...
             'reject_eog'              ,200,[],...      
             'reject_ieeg'             ,[],[],...
             'ignorebadchans_flag'     ,0,{0,1},...
             'prescale_mag'            ,10^15,[],...  % to fT
             'prescale_grad'           ,10^13,[],... % to fT/cm
             'prescale_eeg'            ,10^6,[],...   % to uV
             'prescale_eog'            ,10^6,[],...   % to uV
             'visualreject'            ,[],[],...
             'prefix','proc'           ,[],...
             'events'                  ,[],[],...
			 'conditions'              ,[],[],...
             'badchanfile'             ,[],[],...
             'badchans'                ,[],[],...
             'toilim'                  ,[],[],...
             'verbose'                 ,1,{0,1},...
             'logfile'                 ,[],[],...
             'logfid'                  ,[1],[], ...  
             'filename'                ,[],[],...
             'write_reject_log'        ,0,{0,1},...
             'write_summary_reject_log',0,{0,1},...
             'write_events_file'       ,0,{0,1},...
          	 'write_reject_only_log'   ,0,{0,1},...
             'rootoutdir'              ,pwd,[],...
						},false);

parms = backcompatible(parms,varargin{:});      
data  = ts_checkdata_header(data,'events',parms.events);

if exist(parms.badchanfile,'file') || ~isempty(parms.badchans) || ~isempty(parms.toilim) || ~isempty(parms.events) || ~isempty(parms.rejectfile) || ~isempty(parms.reject_data)
  data = ts_data_selection(data,'rejectfile',parms.rejectfile,'reject_data',parms.reject_data,'badchanfile',parms.badchanfile,'badchans',parms.badchans,'toilim',parms.toilim,'events',parms.events);
end
[datatype,datafield,dataparam] = ts_object_info(data);

nchan  = data.num_sensors;
ncond  = length(data.(datafield));
[allbadtrials{1:ncond}] = deal([]);
[allbadsamps{1:ncond}]  = deal([]);
allbadchans             = [];
%% Automatic Threshold Rejection
if parms.reject_auto_flag
  % determine which channels are which type
  typestring = {data.sensor_info.typestring};
  mag_i   = strcmp ('mag' ,typestring);
  grad_i  = strncmp('grad',typestring,4);
  eeg_i   = strcmp ('eeg' ,typestring);
  eog_i   = strcmp ('eog' ,typestring);
  % generate rejection thresholds for each channel
  reject_thresh         = zeros(nchan,1);
  reject_thresh(mag_i)  = parms.reject_mag/parms.prescale_mag;
  reject_thresh(grad_i) = parms.reject_grad/parms.prescale_grad;
  reject_thresh(eeg_i)  = parms.reject_eeg/parms.prescale_eeg;
  reject_thresh(eog_i)  = parms.reject_eog/parms.prescale_eog;
  reject_thresh(reject_thresh<=0) = inf;
  if parms.ignorebadchans_flag % Do not apply threshold to bad channels added 02/02/2012 -BQR
    reject_thresh(logical([data.sensor_info.badchan])) = inf; 
  end
  badchanlabels = cell(ncond,1);
  excess_max = cell(ncond,1);
  badchanthreshs = cell(ncond,1);
  for c = 1:ncond
    ntrl      = data.(datafield)(c).num_trials;
    nsamp     = length(data.(datafield)(c).time);
    tmpthresh = reject_thresh * ones(1,nsamp);
    rej_idx   = cellfun(@(x)find(abs(data.(datafield)(c).data(:,:,x))>tmpthresh),num2cell(1:ntrl),'uniformoutput',false);
    badtrl_ix = find(~cellfun(@isempty,rej_idx));
    badchan_ix = nan(1,length(badtrl_ix));
    badchanlabels{c,:} = nan(1,length(badtrl_ix));
    excess_max{c,:} = nan(1,length(badtrl_ix));
    badchanthreshs{c,:} = nan(1,length(badtrl_ix));
    for b = 1:length(badtrl_ix)
      trl = badtrl_ix(b);
      [chan_ix,samp_ix] = ind2sub([nchan nsamp],rej_idx{trl});
      chan_ix = unique(chan_ix);
      amps   = max(abs(data.(datafield)(c).data(chan_ix,:,trl)),[],2);  % maximum values for channels exceeding threshold
      ord    = floor(log10(abs(reject_thresh(chan_ix))));               % order of magnitude of rejection thresholds for channels > threshold
      excess = amps- reject_thresh(chan_ix);                            % how much each chan exceeds threshold
      amps   = (excess)./(10.^ord);                                     % how much each chan exceeds threshold scaled by order of magnitude
      [amps,I] = sort(amps,'descend');                                  % normalized excesses sorted largest to smallest
      chan_ix  = chan_ix(I);                                            % channels sorted with the one exceeding the threshold by the most listed first
      excess = excess(I);
      badchan_ix(b) = chan_ix(1);
      excess_max{c,:}(b) = excess(1);
      badchanthreshs{c,:}(b) = reject_thresh(badchan_ix(b));
    end
    badchanlabels{c,:} = {data.sensor_info(badchan_ix).label};  
    allbadtrials{c} = [allbadtrials{c} badtrl_ix];
    if issubfield(data.(datafield),'trial_info.latency')
      allbadsamps{c} = [allbadsamps{c}  data.(datafield)(c).trial_info.latency(badtrl_ix)];
    end
    badchantypes    = typestring(badchan_ix);
    tmp       = data.(datafield)(c).num_rejects;
    tmp.mag   = tmp.mag  + length(find(strmatch('mag' ,badchantypes)));
    tmp.grad  = tmp.grad + length(find(strmatch('grad',badchantypes)));
    tmp.eeg   = tmp.eeg  + length(find(strmatch('eeg' ,badchantypes)));
    tmp.eog   = tmp.eog  + length(find(strmatch('eog' ,badchantypes)));
    data.(datafield)(c).num_rejects = tmp;
    clear tmp
    data.(datafield)(c).trial_info.badtrial(badtrl_ix) = 1;
  end
end

%% Manual Visual Rejection

% ts_browseraw()
% ts_browseraw_graph()

% ts_visual_reject
if parms.reject_manual_flag
  [data, reject_data] = ts_manual_reject(data,varargin{:}); % Andrew's program
%   remove_bad_trials; % remove bad trials before visual rejection
%   ix = find(cellfun(@(x)isequal(x,'method'),varargin)); if any(ix), varargin([ix ix+1])=[]; end
%   ix = find(cellfun(@(x)isequal(x,'metric'),varargin)); if any(ix), varargin([ix ix+1])=[]; end
%   [data,badchans,badtrials] = ts_visual_reject(data,varargin{:},'method',parms.reject_method,'metric',parms.reject_metric);
%   allbadchans  = [allbadchans badchans];
%   for k = 1:ncond
%     if isempty(badtrials{k}), continue; end
%     allbadtrials{k} = [allbadtrials{k} badtrials{k}];
%     allbadsamps{k}  = [allbadsamps{k}  data.(datafield)(k).trial_info.latency(badtrials{k})];
%     data.(datafield)(k).trial_info.number(badtrials{k})        = [];
%     data.(datafield)(k).trial_info.latency(badtrials{k})       = [];
%     data.(datafield)(k).trial_info.badtrial(badtrials{k})      = [];
%     data.(datafield)(k).trial_info.event_code(badtrials{k})    = [];
%     data.(datafield)(k).trial_info.duration(badtrials{k})      = [];
%     data.(datafield)(k).trial_info.datafile(badtrials{k})      = [];
%     data.(datafield)(k).trial_info.events_fnames(badtrials{k}) = [];    
%   end
end

%% Automtic ICA
if parms.ICA_auto_flag
  % perform automatic ICA rejection for blinks
  data = ts_autoICA(data,'ICA_ref_chan',parms.ICA_ref_chan,'chantype',...
                    parms.ICA_chantype,'rescale',parms.ICA_rescale_flag);
end

%% Manual ICA
if parms.ICA_manual_flag
  % perform manual ICA rejection for blinks and/or EKG
  data = ts_manualICA(data,'maxsteps',parms.ICA_maxsteps,'ntrial',parms.ICA_ntrial,'start_trial',parms.ICA_start_trial,...
                'ncomponents',parms.ICA_ncomponents,'chantype',parms.ICA_chantype,'plottype',parms.ICA_plottype,...
                'rescale',parms.ICA_rescale_flag,'sorttrials',parms.ICA_sorttrials,'allconditions',1,'ICA_saveout_flag',0);
end

%% Create and save Reject Data (*.mat, *.log)
if ~exist('reject_data') % if 1
    if ~isempty(parms.rejectfile)
        load(parms.rejectfile)
    else
  reject_data = [];
  reject_data.badchans      = allbadchans;
  reject_data.badchanlabels = {data.sensor_info(allbadchans).label};
  reject_data.badtrials= [];
  for i = 1:length(data.(datafield))
    reject_data.badtrials{i}  = allbadtrials{i};
    reject_data.badsamples{i} = allbadsamps{i}; %data.(datafield)(i).trial_info.latency(allbadsamps{i});
    reject_data.event_code(i) = data.(datafield)(i).event_code;
    flds = fieldnames(data.(datafield)(i).trial_info);
    for f = 1:length(flds)
      reject_data.trial_info.(flds{f}) = data.(datafield)(i).trial_info.(flds{f})(allbadtrials{i});
    end
  end
% else
%   for i = 1:length(data.(datafield))
%     reject_data.badtrials{i} = unique([reject_data.badtrials{i} allbadtrials{i}]);
%     if isfield(reject_data,'badsamples')
%       reject_data.badsamples{i} = unique([reject_data.badsamples{i} badsamples{i}]);
%     else
%       reject_data.badsamples{i} = allbadsamps{i};
%     end
%   end
    end
end

if parms.write_reject_only_log
  % Save badchannel and trials as matlab file
  tmp        = rmfield(data,'epochs');
  tmp.(datafield) = rmfield(data.(datafield),'data');
  save_reject(reject_data,parms,tmp);  
  clear tmp
end

if parms.write_reject_log
  write_reject_log;
end

if parms.write_summary_reject_log
  write_summary_reject_log;
end

if parms.write_events_file
    write_events_file;
end

if ~parms.keepbadtrials_flag
  remove_bad_trials;
end

% calculate noise covariance matrix
if parms.calc_ncov_flag
  data = ts_calc_ncov(data,varargin{:}); 
end

% remove trial_info if it wasn't originally present
if ~issubfield(data.(datafield),'trial_info.latency')
  data.(datafield) = rmfield(data.(datafield),'trial_info');
end
    
%%
  function remove_bad_trials
    for c = 1:ncond
      % remove rejects from data
      idx = find(data.(datafield)(c).trial_info.badtrial == 1);
      if isempty(idx) || all(idx==0), continue; end
      if issubfield(data.(datafield),'trial_info.latency')
        data.(datafield)(c).trial_info.number(idx)        = [];
        data.(datafield)(c).trial_info.latency(idx)       = [];
        data.(datafield)(c).trial_info.badtrial(idx)      = [];
        data.(datafield)(c).trial_info.event_code(idx)    = [];
        data.(datafield)(c).trial_info.duration(idx)      = [];
        data.(datafield)(c).trial_info.datafile(idx)      = [];
        data.(datafield)(c).trial_info.events_fnames(idx) = [];
      end
      data.(datafield)(c).data(:,:,idx) = [];
      data.(datafield)(c).num_trials = data.(datafield)(c).num_trials - length(idx);
    end
    rmix                    = find([data.(datafield).num_trials]==0);
    ncond                   = ncond - length(rmix);
    allbadtrials(rmix)      = [];
    allbadsamps(rmix)       = [];
    data.(datafield)(rmix)  = [];
  end
%% 
   function write_reject_log 
    if ~isfield(data.(datafield),'trial_info')
      fprintf('Failed to write reject log; data structure must contain trial_info.\n');
      return;
    end
    logfile = sprintf('%s/%s_reject.log',parms.rootoutdir,parms.prefix);
 
    % write the log 
    fid     = fopen(logfile,'wt'); 
    fprintf(fid, '%-4s\t%-4s\t%-4s\t%-4s\t%-4s\n','trial','event','latency','type','comment');  
    logcell = cell(5,sum([data.(datafield)(:).num_trials]));
    l_cntr = 0;
     for iepo = 1:ncond
         nlat = length(data.(datafield)(iepo).trial_info.latency)+l_cntr;
         t_cntr = 0;
         c_cntr = 1;
         for ilat = 1+l_cntr:nlat
             logcell{2,ilat} = data.(datafield)(iepo).trial_info.event_code(ilat-l_cntr);
             logcell{3,ilat} = data.(datafield)(iepo).trial_info.latency(ilat-l_cntr);
             if data.(datafield)(iepo).trial_info.badtrial(ilat-l_cntr) &&  ~parms.reject_manual_flag
                logcell{4,ilat} = 'reject';
                 logcell{1,ilat} = [];
                 if parms.reject_auto_flag && ~parms.reject_manual_flag
                     logcell{5,ilat} = ...
                         sprintf('Channel %s exceeds %i by %i',...
                         badchanlabels{iepo}{c_cntr},...
                         badchanthreshs{iepo,:}(c_cntr),...
                         excess_max{iepo,:}(c_cntr));
                 end
             elseif find(reject_data.badtrials{iepo}==(ilat-l_cntr)) & (parms.reject_manual_flag || ~isempty(parms.rejectfile))% added 09.19.11 -BQR
                     logcell{5,ilat} = 'Manually approved rejection';
%                  elseif parms.reject_auto_flag && parms.reject_manual_flag %%% TODO only works if same trials are not affected
%                      if sum(ismember(reject_data(iepo).badtrials{1},data.(datafield)(iepo).trial_info.number(ilat-l_cntr)))
%                         logcell{5,ilat} = 'Manually rejected';
%                         if autodata.(datafield)(iepo).trial_info.badtrial(ilat-l_cntr) ==...
%                                 data.(datafield)(iepo).trial_info.badtrial(ilat-l_cntr) % if not also automaticaly rejected
%                          c_cntr = c_cntr-1;
%                         end
%                      else
%                         logcell{5,ilat} = ...
%                              sprintf('Channel %s exceeds %i by %i',...
%                              badchanlabels{iepo}{c_cntr},...
%                              badchanthreshs{iepo,:}(c_cntr),...
%                              excess_max{iepo,:}(c_cntr)); 
%                      end
                     
             
                 c_cntr = c_cntr+1;
             else
                 logcell{4,ilat} = 'ok';
                 logcell{1,ilat} = 1 + t_cntr;
                 t_cntr = t_cntr + 1;
                 logcell{5,ilat} = '';
             end
         end     
     
         l_cntr = l_cntr + length(data.(datafield)(iepo).trial_info.latency);
     end
   [jnk,sortidx] = sort(cell2mat(logcell(3,:)),'ascend'); % sort by latency
   for ilog =  1:size(logcell,2)
       fprintf(fid,'%-4i\t%-4i\t%-4i\t%-4s\t%-4s\n',...
           logcell{1,sortidx(ilog)},... trial
           logcell{2,sortidx(ilog)},... event
           logcell{3,sortidx(ilog)},... latency
           logcell{4,sortidx(ilog)},... type
           logcell{5,sortidx(ilog)});%  comment
   end
    fclose(fid); 
   end

    function write_summary_reject_log
       if ~isfield(data.(datafield),'trial_info')
          fprintf('Failed to write reject log; data structure must contain trial_info.\n');
          return;
       end
      logfile = sprintf('%s/%s_summary_reject.log',parms.rootoutdir,parms.prefix);

      % write the log
      fid     = fopen(logfile,'wt');
      fprintf(fid,'event_code\tnum_trials\tnum_rejects\tpercent_rej');

      for ievnt = 1:length(data.(datafield))
         event_code = data.(datafield)(ievnt).event_code;
         num_trials = data.(datafield)(ievnt).num_trials;
         num_rejects = cell2mat(struct2cell(data.(datafield)(ievnt).num_rejects));
         num_rejects = sum(num_rejects);
         %tot_trials = num_trials + num_rejects;
         percent_rej = num_rejects / num_trials * 100; 
         fprintf(fid,'\n%g\t%g\t%g\t%g',event_code,num_trials,num_rejects,percent_rej);
      end
      fclose(fid);
    end

    function  write_events_file
      if ~isfield(data.(datafield),'trial_info')
        fprintf('Failed to write reject log; data structure must contain trial_info.\n');
        return;
      end

      logfile = sprintf('%s/%s_events.txt',parms.rootoutdir,parms.prefix);
      fid     = fopen(logfile,'wt');
       fprintf(fid, '%-4s\t %-4s\t %-4s\t %-4s\t\n','type','latency','condition','duration');
      logcell = cell(3,sum([data.(datafield)(:).num_trials]));

      % collect all trial info for all conditions
      nsamp = sum([data.(datafield).num_trials]);
      type  = repmat({'Trigger'},1,nsamp);
      lat   = arrayfun(@(x)x.trial_info.latency,data.(datafield),'uniformoutput',false);
      cond  = arrayfun(@(x)x.trial_info.event_code,data.(datafield),'uniformoutput',false);
      badtrials = arrayfun(@(x)x.trial_info.badtrial,data.(datafield),'uniformoutput',false);
      
      code = [];
      for k = 1:ncond
        code = [code data.(datafield)(k).event_code*ones(1,data.(datafield)(k).num_trials)];
      end
      
      % convert cell arrays into numeric arrays
      lat = [lat{:}];
      cond = [cond{:}];
      badtrials = [badtrials{:}];

      % sort numeric arrays by latency and label rejects
      [lat,I] = sort(lat);
      cond = cond(I);
      code = code(I);
      badtrials = find(badtrials(I));
      [type{badtrials}] = deal('Reject');

%       % remove duplicate samples
%       [lat,I,J] = unique(lat);
%       cond = cond(I);
%       type = type(I);
      
      for k = 1:length(lat)
        fprintf(fid,'%-4s\t%-4i\t%-4i\t%-4i\n',...
             type{k},... type
             lat(k),... latency
             code(k),... condition (cond or code variable?)
             1);                        % duration %% TODO: Is this meaningful?
      end
      fclose(fid);
    
%     l_cntr = 0;
%      for iepo = 1:ncond
%          nlat = length(data.(datafield)(iepo).trial_info.latency)+l_cntr;
%          t_cntr = 0;
%          c_cntr = 1;
%          for ilat = 1+l_cntr:nlat
%              logcell{1,ilat} = data.(datafield)(iepo).trial_info.event_code(ilat-l_cntr);
%              logcell{2,ilat} = data.(datafield)(iepo).trial_info.latency(ilat-l_cntr);
%              if data.(datafield)(iepo).trial_info.badtrial(ilat-l_cntr)
%                 logcell{3,ilat} = 'reject';
%                  c_cntr = c_cntr+1;
%              else
%                  logcell{3,ilat} = 'ok';
%                  t_cntr = t_cntr + 1;
%              end     
%          end
%          l_cntr = l_cntr + length(data.(datafield)(iepo).trial_info.latency);
%      end
%    [jnk,sortidx] = sort(cell2mat(logcell(2,:)),'ascend'); % sort by latency
%    
%    for ilog =  1:size(logcell,2)
%     if strmatch(logcell{3,sortidx(ilog)},'reject')
%         logcell{3,sortidx(ilog)} = 'Reject';
%     else
%            logcell{3,sortidx(ilog)} = 'Trigger';
%     end
%        fprintf(fid,'%-4s\t%-4i\t%-4i\t%-4i\n',...
%            logcell{3,sortidx(ilog)},... type
%            logcell{2,sortidx(ilog)},... latency
%            logcell{1,sortidx(ilog)},... condition
%            1);                        % duration %% TODO: Is this meaningful?
%    end
%     fclose(fid);
%     end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_reject(reject_data,parms,hdr)
  % save reject_data
  % Save reject data
  if ~isempty(reject_data)
    % save info to mat file
    if isempty(parms.filename)
      args           = mmil_parms2args(parms); 
      parms.filename = ts_create_filename('ts_reject',args{:});
    end
    if iscell(parms.filename)
      filename = parms.filename{1};
    else
      filename = parms.filename;
    end
    filename = strrep(filename,'epoch_data_rej','reject_data');
    mmil_logstr(parms,'%s: saving reject data: %s\n',mfilename,filename);
    save(filename,'reject_data');

    % save info to text file (copied from Rajan's code in ts_iEEG_ProcessEEG)
    filename = strrep(filename,'.mat','.txt');
    fid = fopen(filename,'w+');
    badchans = find([hdr.sensor_info.badchan]);  
    fprintf(fid,'Rejected Channels\n\n');
    if ~isempty(badchans)
      for j = 1:length(badchans)
        fprintf(fid,'%s\n',hdr.sensor_info(badchans(j)).label);
      end
    else
      fprintf(fid,'No bad channels.\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'Trial Info\n\n');
    fprintf(fid,'         \t Good\t    Rejected     \t Bad\n');
    fprintf(fid,'Condition\tTrials\tEEG File\tManual\tTrials\n\n');
    for j = 1:length(hdr.(datafield))
      if j > length(reject_data.badtrials), break; end
     fprintf(fid,'%-9s\t%-6s\t%-8s\t%-6s\t%s\n',num2str(j),num2str(hdr.(datafield)(j).num_trials),...
                                                       num2str(hdr.(datafield)(j).num_rejects.eeg),...
                                                       num2str(hdr.(datafield)(j).num_rejects.manual),...
                                                       num2str(reject_data.badtrials{j}));
    end
    fclose(fid);  
  end
end
%%
function parms = backcompatible(parms,varargin)

opt = mmil_args2parms(varargin,{...
     'ICA_ntrials','',[],...   
     'method',[],[],...
     'metric',[],[],...
     'visualreject',[],[],...
     'visualreject_flag',[],[],...
     },false);

if ~isempty(opt.ICA_ntrials)
  parms.ICA_ntrial = opt.ICA_ntrials;
end
if ~isempty(opt.method)
  parms.reject_method = opt.method;
end
if ~isempty(opt.metric)
  parms.reject_metric = opt.metric;
end
if ~isempty(opt.visualreject_flag)
  parms.reject_manual_flag = opt.visualreject_flag;
end
if ~isempty(opt.visualreject)
  parms.reject_manual_flag = opt.visualreject;
end
end