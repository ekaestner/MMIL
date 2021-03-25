function [data_out, badtrials] = ts_freqanalysis_bl_reject (data_in,varargin)
% function [data_out, badtrials] = ts_freqanalysis_bl_reject (data_in,varargin)
%
% Usage: [data_out, badtrials] = ts_freqanalysis_bl_reject (data_in,'option',value...
%
% Rejects trials on a single channel of timefreq data.
%
% Rejection threshold: the threshold times the distance from the median to
% the first decile above the median determines the cut off
%
% Required Inputs:
%
% baseline_data - a valid timefreq data structure
%
% Optional Parameters
%
% threshold - default = 3
% reject_exclude - Not Yet Implemented here (only done by channel)
% 
% Outputs:
%
% baseline_data - data averaged across trials after removing bad trials
% badbaselines  - the list of bad baselines
% 
% Created on 04/17/08 by Rajan Patel
% Last Modified on 04/18/08 by Rajan Patel
%% Check Inputs

if nargin < 1, help(mfilename); end

error(nargoutchk(1,2,nargout,'string')); 

parms = mmil_args2parms(varargin,...
                        {...
                         'bl_threshold',3,[],...
                         'reject_exclude','bychannel',{'inall','groupchannels','bychannel'},...
                         'debug_mode',false ,sort([false true]),...
                         'logfile',[],[],...
                         'logfid',1,[]...
                         },...
                         false); 

errors = ts_checkdata(data_in);
if ~isempty(errors)
    mmil_error(parms,'Errors in supplied data: %s.', sprintf('\t%s\n', errors{:}));
end; 
                     
if ~isfield(data_in,'timefreq')
    mmil_error(parms,'Input data must be timefreq data.');
end

if length(data_in.timefreq) ~= 1
    mmil_error(parms,'Baseline data should be just have one condition.');
end

if ndims(data_in.timefreq.power) ~= 4
    mmil_error(parms,'Baseline data should have trial information for rejection.');
end

if parms.debug_mode
    outdir = fullfile(pwd,'histograms');
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
end


%% Initialize data_out

data_out                = data_in;
data_out.timefreq.power = [];
data_out.timefreq.data  = [];

%% Perform rejection on a channel by channel basis

for ch = 1:length(data_in.sensor_info) 
    badtrials{ch} = [];
    if parms.debug_mode
      % fprintf('%s: Channel %s - number of bad trials:    ',mfilename,data_in.sensor_info(ch).label);
    end
    for f = 1:size(data_in.timefreq.power,3)
       %% Mean across the baseline period for each trial
       rej_power           = squeeze(mean(data_in.timefreq.power(ch,:,f,:),2));
       %% Sort the trials by power
       rej_power           = permute(rej_power,[2 1]);
       [rej_power,trials]  = sort(rej_power);
       %% Get the median
       med_power           = nanmedian(rej_power);
       %% Find the first 10% of trials
       decile        = round(.10 * length(rej_power));
       tmp_badtrials = [];...1:decile; %% DO we include the first decile?? causes a LOT of bad trials
       %% Find the top cut off as 
       cut_off       = med_power + (parms.bl_threshold*(med_power - rej_power(decile)));
       tmp_badtrials = [tmp_badtrials find(rej_power >= cut_off)];
       badtrials{ch} = unique([badtrials{ch} trials(tmp_badtrials)]);
       if parms.debug_mode
           %% DEBUG PLOTTING of distribution and cut offs.           
           %fprintf('\b\b\b%3.0f',length(badtrials{ch}));        
          if ismember(data_in.timefreq.frequencies(f),[10 30 80]) 
           figure    
           hist(rej_power,[0:(cut_off*5)/500:cut_off*5])        
           xlim([0 cut_off*5]);
           %[y,p] = ksdensity(rej_power);
           %plot(p,y,'-b');
           title(sprintf('Channel %s, Frequency %d Hz, # bad trials: %d\n',...
                   data_in.sensor_info(ch).label,...
                   data_in.timefreq.frequencies(f),...
                   length(tmp_badtrials))...
                   ...
                 );
           limits = ylim;
           hold on
           plot(repmat(med_power,1,2),[limits(1) limits(2)],'-.k');
           plot(repmat(rej_power(decile),1,2),[limits(1) limits(2)],'-.c');
           plot(repmat(cut_off,1,2),[limits(1) limits(2)],'-.r');
           hold off          
           outfile = fullfile(outdir,sprintf('%s-%dHz.jpg',data_in.sensor_info(ch).label,data_in.timefreq.frequencies(f)));
           tag = 1;
           while exist(outfile,'file')
               outfile = fullfile(outdir,sprintf('%s-%dHz-%d.jpg',data_in.sensor_info(ch).label,data_in.timefreq.frequencies(f),tag));
               tag = tag+1;
           end
           print('-djpeg','-r300',outfile);
           close all
          end
       end
    end
    %fprintf('\n');
    if size(data_in.timefreq.power,4) == badtrials{ch}
        mmil_logstr(parms,'WARNING: There are no more trials left for the baseline of channel: %s.',data_in.sensor_info(ch).label);
        data_out.timefreq.power(ch,:,:) = nan(1,...
                                              size(data_out.timefreq.power,2),...
                                              size(data_out.timefreq.power,3));
        data_out.timefreq.data (ch,:,:) = nan(1,...
                                              size(data_out.timefreq.data,2),...
                                              size(data_out.timefreq.data,3));
    else
      data_out.timefreq.power(ch,:,:) = squeeze(mean(data_in.timefreq.power(ch,:,:,...
                                                   setdiff(1:size(data_in.timefreq.power,4),...
                                                   badtrials{ch})),4));
      data_out.timefreq.data (ch,:,:) = squeeze(mean(data_in.timefreq.data(ch,:,:,...
                                                   setdiff(1:size(data_in.timefreq.power,4),...
                                                   badtrials{ch})),4));
    end
    %mmil_logstr(parms,'Removed %.0f trials from channel %s.',length(badtrials{ch}),data_in.sensor_info(ch).label);                                            
end

% Do not return a cell of length == 1
if length(badtrials)==1, badtrials = badtrials{1}; end
