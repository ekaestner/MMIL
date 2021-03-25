function [data_out,varargout] = ts_timefreq_reject (data_in,varargin)

% Usage: [timefreq_data] = ts_timefreq_reject (timefreq_data,'option',value ...
%
% Performs artifact rejection in the timefrequency domain using the options
% specified.
%

%% Verify Inputs

if (~mmil_check_nargs(nargin, 1))
    return;
end

opt = mmil_args2parms(varargin,...
                      {...
                       'logfile',[],[],...
                       'logfid',1,[]...
                      },...
                    false);
                  
errors = ts_checkdata(data_in);

if (length(errors)~=0)
 mmil_error(opt, 'Errors supplied data: %s.', sprintf('\t%s\n', errors{:}));
end;
                    
if ~isfield(data_in,'timefreq')
 mmil_error(opt,'Input data must be timefreq_data.\n');
end  

error(nargoutchk(1,3,nargout,'string'));

%% Initialize Data Out

data_out = data_in;
for i = 1:length(data_out.epochs)
  data_out.timefreq(i).power = [];   % clear memory
end

badchannels = [];
badtrials   = [];
goodtrials  = [];

%% Rejections

for j = 1:length(data_in.timefreq)
  mmil_logstr(opt,'Rejecting artifacts for event %s.',num2str(timefreq_data.timefreq(j).event_code));
  
end

%% Assign Extra Output

if nargout > 1, varargout{1} = badchannels; end
if nargout > 2, varargout{2} = badtrials;   end
