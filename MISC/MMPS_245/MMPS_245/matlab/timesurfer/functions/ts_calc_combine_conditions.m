function data_out = ts_calc_combine_conditions (data_in, varargin)

% Use: [data_out] = ts_calc_comb_conditions (data_in, 'conditions', [1 4 3], 'calc', 'weighted')
%
% Purpose:
%
% Function to create weighted averages, regular averages, or subtraction
% conditions for existing avg_data.  For epochs it will combine multiple
% events into a single new event no matter what is specified in 'calc'.
% The new condition will be added as the last condition with an event_code
% = condition number or you can specify an eventcode.
%
% Required input:
% 
%   data_in    - avg_data,epoch_data,or timefreq_data TimeSurfer data structure
%   conditions - the list of conditions (not event codes) to be used for averaging; 
%                in the case of subtraction conditions there should be only two 
%                conditions listed
%                (eg: to get condition 1 - condition 2, conditions = [1 2])
%   events     - list of event codes to be used (this option will override
%                'conditions' option).
%
% Optional input:
% 
%   calc       - how to calculate the new combined condition {default = 'weighted'}
%               Only valid for use with avg_data:
%               'weighted'     - computes a weighted average
%               'avg'          - computes a straight average
%               'sum'          - computes a sum
%               'sub'          - computes a subtraction condition
%   eventcode  - event_code for the new condition {default: eventcode =
%                condition number)
%
% Output:
%
%   data_out   - avg_data or epoch_data depending on what was supplied;
%                the new structure will have the new condition as the last
%                condition in the data structure, retaining all old
%                conditions.
%
% Created:       10/08/07 - Rajan Patel
% Last Modified: 01/06/11 - Jason Sherfey (added trial_info)
%

%% Check inputs

if ~mmil_check_nargs(nargin, 1), return; end;

opt = mmil_args2parms(varargin,{...
  'conditions', [],[],...
  'events'    , [],[],...
  'calc'      ,'weighted',{'weighted','avg','sum','sub'},...
  'eventcode' , [],[],...
});
                    
type = opt.calc;

if ismember('averages',fieldnames(data_in))
  datafield = 'averages';
elseif ismember('epochs',fieldnames(data_in))
  datafield = 'epochs';
elseif ismember('timefreq',fieldnames(data_in))
  datafield = 'timefreq';
else
  error('invalid data supplied');
end

data_out = data_in; % all common variables set and old conditions left intact

if ~isempty(opt.events)  % if event codes given convert to condition numbers
  if ~isempty(opt.conditions)
    fprintf('%s: WARNING: Event values given will be used instead of condition numbers.\n',mfilename);
  end
  conditions = [];
  evcodes = cell2mat({data_in.(datafield).event_code});
  for k=1:length(opt.events) % convert events to conditions without changing order
    c = find(opt.events(k)==evcodes);
    if ~isempty(c)
      try
        conditions = cat(1,conditions,c); 
      catch
        conditions = cat(2,conditions,c); 
      end;
    end
  end;
  if length(conditions) ~= length(opt.events)                      % make sure found all events
    fprintf('%s: WARNING: some event code(s) [ %s] not found\n',mfilename,sprintf('%d ',opt.events));
    return;
  end
elseif ~isempty(opt.conditions)
  conditions = opt.conditions;   
end

if ~exist('conditions','var')
  error('no condition numbers or event codes supplied');
end

for i = 1:length(conditions)                                          % Quick check if conditions given are valid
  if conditions(i) > length(data_in.(datafield))
    error('%s is an invalid condition for the given data',num2str(conditions(i)));
  end
end

if isfield(opt,'eventcode'), eventcode = opt.eventcode;
else eventcode = []; end

%% Calculate Combined Averages

if strcmp(datafield,'averages')
    switch type
        case {'weighted' 'avg' 'sum'}
            numerator = 0;
            stdev_num = zeros(length(data_in.sensor_info),length(data_in.(datafield)(1).time));
            mag    = 0;
            grad   = 0;
            eog    = 0;
            eeg    = 0;
            manual = 0;
            skip   = 0;
            ntrial = 0;
            for j = 1:length(conditions) 
                switch type
                  case 'weighted'
                    numerator = numerator + (data_in.(datafield)(conditions(j)).data.*data_in.(datafield)(conditions(j)).num_trials);
                  case {'avg' 'sum'}
                    numerator = numerator + (data_in.(datafield)(conditions(j)).data);
                end
                stdev_num   = stdev_num + ...
                              ((data_in.(datafield)(conditions(j)).num_trials-1).* ...
                               (data_in.(datafield)(conditions(j)).stdev.^2));
                ntrial = ntrial + data_in.(datafield)(conditions(j)).num_trials;
                mag    = mag + data_in.(datafield)(conditions(j)).num_rejects.mag;
                grad   = grad + data_in.(datafield)(conditions(j)).num_rejects.grad;
                eog    = eog + data_in.(datafield)(conditions(j)).num_rejects.eog;
                eeg    = eeg + data_in.(datafield)(conditions(j)).num_rejects.eeg;
                manual = manual + data_in.(datafield)(conditions(j)).num_rejects.manual;
                skip   = skip + data_in.(datafield)(conditions(j)).num_rejects.skip;
            end
            data_out.(datafield)(end+1).time             = data_in.(datafield)(1).time;
            data_out.(datafield)(end).num_trials         = ntrial;
            if ~isempty(eventcode), data_out.(datafield)(end).event_code = eventcode;
            else                    data_out.(datafield)(end).event_code = length(data_out.(datafield));
            end
            switch type
              case 'weighted'
                denominator = ntrial;
              case 'avg'
                denominator = length(conditions);
              otherwise
                denominator = 1;
            end
            data_out.(datafield)(end).data               = numerator./denominator;
            data_out.(datafield)(end).num_rejects.mag    = mag;
            data_out.(datafield)(end).num_rejects.grad   = grad;
            data_out.(datafield)(end).num_rejects.eog    = eog;
            data_out.(datafield)(end).num_rejects.eeg    = eeg;
            data_out.(datafield)(end).num_rejects.manual = manual;
            data_out.(datafield)(end).num_rejects.skip   = skip;
            %% todo: stdev is wrong for weighted average
            data_out.(datafield)(end).stdev              = sqrt(stdev_num./(denominator-length(conditions)));
           
        case 'sub'
          if length(conditions)~=2,
            error('for subtractions, must have only two condition numbers');
          end
            data_out.(datafield)(end+1).time        = data_in.(datafield)(1).time;
            % when subtracting two conditions, the stdev increases by a factor of sqrt(2)
            %   which would be equivalent to reducing the number of trials by half
            %
            % with different number of conditions for the two conditions, maybe it will work to
            %   use the average number of trials and then divide by 2 (or just sum and divide by 4)
            % This is only appropriate when there are the same number of conditions
            data_out.(datafield)(end).num_trials    =  round((data_in.(datafield)(conditions(1)).num_trials + ...
                                                              data_in.(datafield)(conditions(2)).num_trials)/4);
            fprintf('%s: number of trials for first condition: %d\n',mfilename,data_in.(datafield)(conditions(1)).num_trials);
            fprintf('%s: number of trials for second condition: %d\n',mfilename,data_in.(datafield)(conditions(2)).num_trials);
            fprintf('%s: effective number of trials for new subtraction condition: %d\n',mfilename,data_out.(datafield)(end).num_trials);
            if ~isempty(eventcode), data_out.(datafield)(end).event_code = eventcode;
            else                    data_out.(datafield)(end).event_code = length(data_out.(datafield));
            end
            data_out.(datafield)(end).data               = data_in.(datafield)(conditions(1)).data -...
                                                           data_in.(datafield)(conditions(2)).data;
            data_out.(datafield)(end).num_rejects.mag    = data_in.(datafield)(conditions(1)).num_rejects.mag + ...
                                                           data_in.(datafield)(conditions(2)).num_rejects.mag;
            data_out.(datafield)(end).num_rejects.grad   = data_in.(datafield)(conditions(1)).num_rejects.grad + ...
                                                           data_in.(datafield)(conditions(2)).num_rejects.grad;
            data_out.(datafield)(end).num_rejects.eog    = data_in.(datafield)(conditions(1)).num_rejects.eog + ...
                                                           data_in.(datafield)(conditions(2)).num_rejects.eog;
            data_out.(datafield)(end).num_rejects.eeg    = data_in.(datafield)(conditions(1)).num_rejects.eeg + ...
                                                           data_in.(datafield)(conditions(2)).num_rejects.eeg;
            data_out.(datafield)(end).num_rejects.manual = data_in.(datafield)(conditions(1)).num_rejects.manual + ...
                                                           data_in.(datafield)(conditions(2)).num_rejects.manual;
            data_out.(datafield)(end).num_rejects.skip   = data_in.(datafield)(conditions(1)).num_rejects.skip + ...
                                                           data_in.(datafield)(conditions(2)).num_rejects.skip;  
            data_out.(datafield)(end).stdev              = sqrt(...
              ((data_out.(datafield)(conditions(1)).stdev.^2)./data_out.(datafield)(conditions(1)).num_trials) + ...
              ((data_out.(datafield)(conditions(2)).stdev.^2)./data_out.(datafield)(conditions(2)).num_trials));
        otherwise
          error('unknown type of calculation');
    end
    
%% Combine multiple epochs

elseif strcmp(datafield,'epochs')
    data_out.(datafield)(end+1).time             = data_in.(datafield)(1).time;
    data_out.(datafield)(end).num_rejects.mag    = 0;
    data_out.(datafield)(end).num_rejects.grad   = 0;
    data_out.(datafield)(end).num_rejects.eog    = 0;
    data_out.(datafield)(end).num_rejects.eeg    = 0;
    data_out.(datafield)(end).num_rejects.manual = 0;
    data_out.(datafield)(end).num_rejects.skip   = 0;
    data_out.(datafield)(end).num_trials         = 0;
    data_out.(datafield)(end).data               = [];
    if ~isempty(eventcode), data_out.(datafield)(end).event_code = eventcode;
    else                    data_out.(datafield)(end).event_code = length(data_out.(datafield));
    end
    for j = 1:length(conditions)
        data_out.(datafield)(end).data               = cat (3, data_out.(datafield)(end).data, ...
                                                               data_in.(datafield)(conditions(j)).data);
        data_out.(datafield)(end).num_rejects.mag    = data_out.(datafield)(end).num_rejects.mag + ...
                                                       data_in.(datafield)(conditions(j)).num_rejects.mag;
        data_out.(datafield)(end).num_rejects.grad   = data_out.(datafield)(end).num_rejects.grad + ...
                                                       data_in.(datafield)(conditions(j)).num_rejects.grad;
        data_out.(datafield)(end).num_rejects.eog    = data_out.(datafield)(end).num_rejects.eog + ...
                                                       data_in.(datafield)(conditions(j)).num_rejects.eog;
        data_out.(datafield)(end).num_rejects.eeg    = data_out.(datafield)(end).num_rejects.eeg + ...
                                                       data_in.(datafield)(conditions(j)).num_rejects.eeg;
        data_out.(datafield)(end).num_rejects.manual = data_out.(datafield)(end).num_rejects.manual + ...
                                                       data_in.(datafield)(conditions(j)).num_rejects.manual;
        data_out.(datafield)(end).num_rejects.skip   = data_out.(datafield)(end).num_rejects.skip + ...
                                                       data_in.(datafield)(conditions(j)).num_rejects.skip;
        data_out.(datafield)(end).num_trials         = data_out.(datafield)(end).num_trials + ...
                                                       data_in.(datafield)(conditions(j)).num_trials;
        data_in.(datafield)(conditions(j)).data = [];                                             
    end 

%% Combine multiple tf
elseif strcmp(datafield,'timefreq')
   datflag = 0; powflag = 0; cpxflag = 0;
   if issubfield(data_in,'timefreq.data'),  datflag = 1; datnum = 0; zparam = 'data';  end
   if issubfield(data_in,'timefreq.power'), powflag = 1; pownum = 0; zparam = 'power'; end
   if issubfield(data_in,'timefreq.cmplx'), cpxflag = 1; cpxnum = 0; zparam = 'cmplx'; end
   if cpxflag && ~powflag
     for cc = 1:length(data_in.timefreq)
       data_in.timefreq(cc).power = 2*abs(data_in.timefreq(cc).cmplx).^2;
     end
     data_out = data_in;
     powflag = 1; pownum = 0; zparam  = 'power';
   end
   switch type
       case {'weighted','avg','sum'}
%             pow_num = 0;
%             data_num = 0;
            mag    = 0;
            grad   = 0;
            eog    = 0;
            eeg    = 0;
            manual = 0;
            skip   = 0;
            ntrial = 0;
            for j = 1:length(conditions)               
                if length(data_in.(datafield)(conditions(j)).num_trials)==1
                    weights = data_in.(datafield)(conditions(j)).num_trials ...
                              * ones(size(data_in.(datafield)(conditions(j)).(zparam)));
                else
                    weights = data_in.(datafield)(conditions(j)).num_trials;
                end
%                     data_in.(datafield)(conditions(j)).num_trials = repmat(data_in.(datafield)(conditions(j)).num_trials,...
%                                                                            1,data_in.num_sensors);
%                 end
                switch type
                  case 'weighted'
%                     wt_fact  = repmat(permute(data_in.(datafield)(conditions(j)).num_trials,[2 1]),...
%                                               [ 1 ... 
%                                                 size(data_in.(datafield)(conditions(j)).power,2)...
%                                                 size(data_in.(datafield)(conditions(j)).power,3) ] );
                    if datflag, datnum = datnum + (data_in.(datafield)(conditions(j)).data.*weights);  end
                    if powflag, pownum = pownum + (data_in.(datafield)(conditions(j)).power.*weights); end
                    if cpxflag, cpxnum = cpxnum + (data_in.(datafield)(conditions(j)).cmplx.*weights); end
%                     data_num = data_num + (data_in.(datafield)(conditions(j)).data.*wt_fact);
%                     pow_num  = pow_num  + (data_in.(datafield)(conditions(j)).power.*wt_fact);
                  case {'avg' 'sum'}
                    if datflag, datnum = datnum + (data_in.(datafield)(conditions(j)).data);  end
                    if powflag, pownum = pownum + (data_in.(datafield)(conditions(j)).power); end
                    if cpxflag, cpxnum = cpxnum + (data_in.(datafield)(conditions(j)).cmplx); end                    
%                     data_num = data_num + (data_in.(datafield)(conditions(j)).data);
%                     pow_num  = pow_num + (data_in.(datafield)(conditions(j)).power);
                end
                ntrial = ntrial + data_in.(datafield)(conditions(j)).num_trials;
                mag    = mag + data_in.(datafield)(conditions(j)).num_rejects.mag;
                grad   = grad + data_in.(datafield)(conditions(j)).num_rejects.grad;
                eog    = eog + data_in.(datafield)(conditions(j)).num_rejects.eog;
                eeg    = eeg + data_in.(datafield)(conditions(j)).num_rejects.eeg;
                manual = manual + data_in.(datafield)(conditions(j)).num_rejects.manual;
                skip   = skip + data_in.(datafield)(conditions(j)).num_rejects.skip;
            end
            data_out.(datafield)(end+1).time               = data_in.(datafield)(1).time;
            data_out.(datafield)(end).frequencies          = data_in.(datafield)(1).frequencies;
            if ~isempty(eventcode), data_out.(datafield)(end).event_code = eventcode;
            else                    data_out.(datafield)(end).event_code = length(data_out.(datafield));
            end
            if length(data_in.(datafield)(conditions(j)).num_trials)==1
              data_out.(datafield)(end).num_trials = ntrial(1);
            else
              data_out.(datafield)(end).num_trials = ntrial;
            end            
            switch type
              case 'weighted'
                denominator = ntrial;
              case 'avg'
%                 %% todo: code was changed without testing (DH)
%                 denominator = length(conditions)*ones(1,data_in.num_sensors);
                denominator = length(conditions); 
                  % changed to scalar (JS) 
                  % why vectorize? - isn't this the same for all channels?
              otherwise
%                 denominator = ones(1,data_in.num_sensors);
                denominator = 1; % changed to scalar (JS)
            end
%             denominator = repmat(permute(denominator,[2 1]),[1 size(pow_num,2) size(pow_num,3)]);
            if datflag, data_out.(datafield)(end).data  = datnum ./ denominator; end
            if powflag, data_out.(datafield)(end).power = pownum ./ denominator; end
            if cpxflag, data_out.(datafield)(end).cmplx = cpxnum ./ denominator; end
%             data_out.(datafield)(end).data               = data_num./denominator;
%             data_out.(datafield)(end).power              = pow_num./denominator;
            data_out.(datafield)(end).num_rejects.mag    = mag;
            data_out.(datafield)(end).num_rejects.grad   = grad;
            data_out.(datafield)(end).num_rejects.eog    = eog;
            data_out.(datafield)(end).num_rejects.eeg    = eeg;
            data_out.(datafield)(end).num_rejects.manual = manual;
            data_out.(datafield)(end).num_rejects.skip   = skip;
       case 'sub'
           if length(conditions)~=2,
                error('for subtraction, must have only two condition numbers');
           end          
%             datflag = 0; powflag = 0; cpxflag = 0;
%             if issubfield(data_in,'timefreq.data'),  datflag = 1; end
%             if issubfield(data_in,'timefreq.power'), powflag = 1; end
%             if issubfield(data_in,'timefreq.cmplx'), cpxflag = 1; end            
                data_out.(datafield)(end+1).time        = data_in.(datafield)(1).time;
                data_out.(datafield)(end).frequencies   = data_in.(datafield)(1).frequencies;
                data_out.(datafield)(end).num_trials    = round(mean(cat(1,data_in.(datafield)(conditions(1)).num_trials,...
                                                                  data_in.(datafield)(conditions(2)).num_trials),1)/2);
                fprintf('%s: Effective number of trials for new subtraction condition: %d\n',mfilename,data_out.(datafield)(end).num_trials);
                                                              
                if ~isempty(eventcode), data_out.(datafield)(end).event_code = eventcode;
                else                    data_out.(datafield)(end).event_code = length(data_out.(datafield));
                end
                if datflag, data_out.(datafield)(end).data   = data_in.(datafield)(conditions(1)).data -...
                                                               data_in.(datafield)(conditions(2)).data; 
                end
                if powflag, data_out.(datafield)(end).power  = data_in.(datafield)(conditions(1)).power -...
                                                               data_in.(datafield)(conditions(2)).power;    
                end
                if cpxflag, data_out.(datafield)(end).cmplx  = data_in.(datafield)(conditions(1)).cmplx -...
                                                               data_in.(datafield)(conditions(2)).cmplx;    
                end                
                data_out.(datafield)(end).num_rejects.mag    = data_in.(datafield)(conditions(1)).num_rejects.mag + ...
                                                               data_in.(datafield)(conditions(2)).num_rejects.mag;
                data_out.(datafield)(end).num_rejects.grad   = data_in.(datafield)(conditions(1)).num_rejects.grad + ...
                                                               data_in.(datafield)(conditions(2)).num_rejects.grad;
                data_out.(datafield)(end).num_rejects.eog    = data_in.(datafield)(conditions(1)).num_rejects.eog + ...
                                                               data_in.(datafield)(conditions(2)).num_rejects.eog;
                data_out.(datafield)(end).num_rejects.eeg    = data_in.(datafield)(conditions(1)).num_rejects.eeg + ...
                                                               data_in.(datafield)(conditions(2)).num_rejects.eeg;
                data_out.(datafield)(end).num_rejects.manual = data_in.(datafield)(conditions(1)).num_rejects.manual + ...
                                                               data_in.(datafield)(conditions(2)).num_rejects.manual;
                data_out.(datafield)(end).num_rejects.skip   = data_in.(datafield)(conditions(1)).num_rejects.skip + ...
                                                               data_in.(datafield)(conditions(2)).num_rejects.skip;  
       otherwise
         error('unknown type of calculation');    
   end
   
end
if ~isfield(data_in.(datafield),'trial_info')
  return;
end
% combine trial_info
trial_arr = [data_in.(datafield)(conditions).trial_info];
trl   = trial_arr(1);
for i = 2:length(trial_arr)
  trl.datafile      = {trl.datafile{:} trial_arr(i).datafile{:}};
  trl.events_fnames = {trl.events_fnames{:} trial_arr(i).events_fnames{:}};
  trl.latency       = [trl.latency trial_arr(i).latency];
  trl.event_code    = [trl.event_code trial_arr(i).event_code];
  trl.duration      = [trl.duration trial_arr(i).duration];
  trl.badtrial      = [trl.badtrial trial_arr(i).badtrial];
  trl.number        = [trl.number trial_arr(i).number];
  if isfield(trl,'epoch_num')
    trl.epoch_num   = [trl.epoch_num trial_arr(i).epoch_num];
  end
end
data_out.(datafield)(end).trial_info = trl;