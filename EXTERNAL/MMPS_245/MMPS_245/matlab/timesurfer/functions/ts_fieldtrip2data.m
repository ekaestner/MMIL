function data = ts_fieldtrip2data(ft_data, outtype, template_data)
%  function data = ts_fieldtrip2data(ft_data, outtype, template_data)
%
%  Creates a TimeSurfer data structure from FieldTrip data.  The template
%  data is a TimeSurfer data structure that was used to initially create
%  the FieldTrip data. 
%
%  ft_data should be a cell array of ft_data structures if you want to
%  convert more than one ft_data structure into one TimeSurfer data
%  structure.
%
%  Using no template_data TimeSurfer data set will set certain parameters
%  to defaults
% 
%  Created by Ben C. ???
%  Last Modified by Jason Sherfey on 25-Jan-2011

%  04/10/08 - fixed timefreq conv broken by last mod
%  05/15/08 - Added support for avg_data and epoch_data
%  03/21/09 - added support for cross-spectra
%  01/25/11 - added support for trial_info
%  TO DO:
%
%    - Ability to specify which conditions or events the ft_data sets
%      belong to and then search for them in the template_data
%    - Input list of bad chans or text file 

%  05/15/08 - Correctly exanded out stdev to account for missing chans

%% Check Inputs

mmil_check_nargs(nargin, 2);

if ~iscell(ft_data), ft_data = {ft_data}; end

if exist('template_data','var')
    % Copy fields, bit-by-bit.
    data.num_sensors = template_data.num_sensors;
    data.sfreq       = template_data.sfreq;
    data.sensor_info = template_data.sensor_info;
    data.coor_trans  = template_data.coor_trans;
    data.noise       = template_data.noise;
    if (isfield(template_data,'epochs'))
        DATA_FIELD = 'epochs';
    elseif (isfield(template_data,'averages'))
        DATA_FIELD = 'averages';
    elseif (isfield(template_data,'timefreq'))
        DATA_FIELD = 'timefreq';
    else
        mmil_error(parms, 'Unknown input datatype');
    end;
else % Use the information from the first ft_data set
    data.num_sensors = length(ft_data{1}.label);
    data.sfreq       = ft_data{1}.fsample;
    for s = 1:length(ft_data{1}.label)
      if iscell(ft_data{1}.label(s))
          data.sensor_info(s).label = ft_data{1}.label{s};
      else
          data.sensor_info(s).label = ft_data{1}.label(s);
      end
      if isfield(ft_data,'elec')
          if strmatch(avg_data.sensor_info(s).label,ft_data{1}.elec.label,'exact')
              data.sensor_info(s).typestring = 'eeg';
              data.sensor_info(s).type       = 1;
              data.sensor_info(s).kind       = 2;
          elseif strmatch(avg_data.sensor_info(s).label,ft_data{1}.grad.label,'exact')
              data.sensor_info(s).typestring = 'grad1';
              data.sensor_info(s).type       = 0;  %% FIX!!
              data.sensor_info(s).kind       = 0;  %% FIX!!
          end
      else %% make it a default eeg channel
          data.sensor_info(s).typestring = 'eeg';
          data.sensor_info(s).type       = 1;
          data.sensor_info(s).kind       = 2;
      end
      data.sensor_info(s).lognum     = s;
      data.sensor_info(s).loc        = eye(4);
      data.sensor_info(s).badchan    = 0;
    end
    data.coor_trans.device2head = eye(4);
    data.coor_trans.mri2head    = [];
    data.noise.num_trials       = [];
    data.noise.num_samples      = [];
    data.noise.covar            = zeros(data.num_sensors,data.num_sensors);
end


% FT data sampling frequency should take precedence
if isfield(ft_data{1},'fsample')
    data.sfreq       = ft_data{1}.fsample;
end


% Backward compatable
if strcmp(outtype,'average')
    outtype = 'averages';
elseif strcmp(outtype,'epoch')
    outtype = 'epochs';
end

%% Convert Data


for f = 1:length(ft_data)
    if exist('template_data','var')
        % Base propreties for each condition
        if (length(template_data.(DATA_FIELD)) < f)
            fprintf('%s: WARNING: There are not enough conditions in the template data.\n',mfilename);
            data.(outtype)(f).event_code         = f;
            data.(outtype)(f).num_trials         = template_data.(DATA_FIELD)(1).num_trials;
            data.(outtype)(f).time               = template_data.(DATA_FIELD)(1).time;
            data.(outtype)(f).num_rejects.mag    = 0;
            data.(outtype)(f).num_rejects.grad   = 0;
            data.(outtype)(f).num_rejects.eog    = 0;
            data.(outtype)(f).num_rejects.eeg    = 0;
            data.(outtype)(f).num_rejects.manual = 0;
            data.(outtype)(f).num_rejects.skip   = 0;
        else
            data.(outtype)(f).event_code  = template_data.(DATA_FIELD)(f).event_code;
            data.(outtype)(f).num_trials  = template_data.(DATA_FIELD)(f).num_trials;
            data.(outtype)(f).time        = template_data.(DATA_FIELD)(f).time;
            data.(outtype)(f).num_rejects = template_data.(DATA_FIELD)(f).num_rejects;
            if isfield(template_data.(DATA_FIELD),'trial_info')
              data.(outtype)(f).trial_info= template_data.(DATA_FIELD)(f).trial_info;
            end
        end;
    else
        data.(outtype)(f).event_code         = f;
        data.(outtype)(f).num_rejects.mag    = 0;
        data.(outtype)(f).num_rejects.grad   = 0;
        data.(outtype)(f).num_rejects.eog    = 0;
        data.(outtype)(f).num_rejects.eeg    = 0;
        data.(outtype)(f).num_rejects.manual = 0;
        data.(outtype)(f).num_rejects.skip   = 0;
    end
    %% FT data time should take precendence
    if iscell(ft_data{f}.time)
      data.(outtype)(f).time          = ft_data{f}.time{1};
    else
      data.(outtype)(f).time          = ft_data{f}.time;
    end
    switch (outtype)
        case 'averages'
           if isfield(ft_data{f},'dimord') 
            switch (ft_data{f}.dimord)
                case 'chan_time'
                    if isfield(ft_data{f},'numsamples')
                        data.(outtype)(f).num_trials = length(ft_data{f}.numsamples);
                    end
                    if isfield(ft_data{f},'avg')
                        data.(outtype)(f).data = ft_data{f}.avg;
                    else
                        if isfield(ft_data{f},'trial')
                            data.(outtype)(f).data = ft_data{f}.trial;
                        else
                            mmil_error('Could not find average data in the field trip data.');
                        end
                    end
                case {'rpt_chan_time','trial_chan_time'}
                   data.(outtype)(f).num_trials = size(ft_data{f}.trial,1);
                   data.(outtype)(f).data       = squeeze(mean(ft_data{f}.trial,1));
                   if ~isfield(ft_data{f},'var')
                       ft_data{f}.var = var(ft_data{f}.trial,0,1);
                   end
                otherwise
                      mmil_error('Unrecognized dimord in field trip data: %s.',ft_data{f}.dimord);
            end
           elseif iscell(ft_data{f}.trial)
               data.(outtype)(f).num_trials = length(ft_data{f}.trial);
              for t = 1:length(ft_data{f}.trial)                 
                   data.(outtype)(f).data(t,:,:)  = ft_data{f}.trial{t};
              end
              data.(outtype)(f).data = squeeze(mean(data.(outtype)(f).data,1));
           else
               mmil_error('There is something wrong with the field trip data.');
           end
           if isfield(ft_data{f},'var')
               data.(outtype)(f).stdev = sqrt(ft_data{f}.var);
           else
               data.(outtype)(f).stdev = zeros(length(data.sensor_info),length(data.(outtype)(f).time));
           end
        case 'epochs'
            if isfield(ft_data{f},'dimord')
                switch (ft_data{f}.dimord)
                    case 'chan_time'
                        mmil_error('There is only average data in the fieldtrip data set.');
                    case {'rpt_chan_time', 'trial_chan_time'}
                        data.(outtype)(f).num_trials = size(ft_data{f}.trial,1);
                        data.(outtype)(f).data       = permute(ft_data{f}.trial,[2 3 1]);
                    otherwise
                        mmil_error('Unrecognized dimord in the field trip data: %s.',ft_data{f}.dimord);
                end
            else
              data.(outtype)(f).num_trials = length(ft_data{f}.trial);
              for t = 1:length(ft_data{f}.trial)                 
                   data.(outtype)(f).data(:,:,t)  = ft_data{f}.trial{t};
              end
            end
        case 'timefreq'
            % Set frequency
            data.(outtype)(f).frequencies   = ft_data{f}.freq;
            if isfield(ft_data{f},'labelcmb')
              data.(outtype)(f).labelcmb = ft_data{f}.labelcmb;
            end
            switch (ft_data{f}.dimord)
                case 'chan_freq_time'
                    % save off coefficients
                    % save off power       
                    if isfield(ft_data{f},'spctrcmpx')
                      data.(outtype)(f).cmplx = ft_data{f}.spctrcmpx; 
                      data.(outtype)(f).cmplx = permute(data.(outtype)(f).cmplx, [1 3 2]);
%                         data.(outtype)(f).power = 2*abs(data.(outtype)(f).cmplx).^2;                      
                    end
                    if (isfield(ft_data{f}, 'powspctrm'))
                        data.(outtype)(f).power = permute(ft_data{f}.powspctrm, [1 3 2]);
                        data.(outtype)(f).power = 2 * data.(outtype)(f).power / data.sfreq;
                    end
                    if (isfield(ft_data{f}, 'coeff'))
                        data.(outtype)(f).data  = permute(ft_data{f}.coeff, [1 3 2]);         
                        data.(outtype)(f).power = permute(2*abs(data.(outtype)(f).data).^2, [1 3 2]);
                    end
%                     else
%                         mmil_error('Could not find any data in the FieldTrip timefrequency data.\n');
%                     end
                    if isfield(ft_data{f},'crsspctrm'), data.(outtype)(f).cross = permute(ft_data{f}.crsspctrm,[1 3 2]); end
                    if isfield(ft_data{f},'cohspctrm'), data.(outtype)(f).coh   = permute(ft_data{f}.cohspctrm,[1 3 2]); end                    
                    if isfield(ft_data{f},'plvspctrm'), data.(outtype)(f).plv   = permute(ft_data{f}.plvspctrm,[1 3 2]); end                                        
                case 'rpttap_chan_freq_time'
                    % we have tapers & trials; want to average.
                    data.(outtype)(f).data = squeeze(sum(ft_data{f}.fourierspctrm, 1));
                    data.(outtype)(f).power = 2*abs(data.(outtype)(f).data).^2;
                    data.(outtype)(f).data  = permute(data.(outtype)(f).data,  [1 3 2]);
                    data.(outtype)(f).power = permute(data.(outtype)(f).power, [1 3 2]);
                % Added by Jason
                case 'rpt_chan_freq_time'
                    % save individual trial information 
                    if isfield(ft_data{f},'spctrcmpx')
                        data.(outtype)(f).cmplx = permute(ft_data{f}.spctrcmpx, [2 4 3 1]); 
                    end
                    if (isfield(ft_data{f}, 'powspctrm'))
                        data.(outtype)(f).num_trials = size(ft_data{f}.powspctrm,1);
                        data.(outtype)(f).power      =  permute(ft_data{f}.powspctrm, [2 4 3 1]);  % chan time freq trial
                        data.(outtype)(f).power      = 2 * data.(outtype)(f).power / data.sfreq;
%                         data.(outtype)(f).data  = nan(size(ft_data{f}.powspctrm));
                    end
                    if (isfield(ft_data{f}, 'coeff'))
                        data.(outtype)(f).num_trials = size(ft_data{f}.coeff,1);
                        data.(outtype)(f).power      = permute(2*abs(ft_data{f}.coeff).^2, [2 4 3 1]);  % chan time freq trial;                        
%                         data.(outtype)(f).data       = squeeze(sum(ft_data{f}.coeff,1));
                    end
%                     else
%                         mmil_error('Could not find any data in the FieldTrip timefrequency data.\n');
%                     end; 
%                     data.(outtype)(f).data = permute(data.(outtype)(f).data, [2 4 3 1]);  % chan time freq trial                    
                    if isfield(ft_data{f},'crsspctrm'), data.(outtype)(f).cross = permute(ft_data{f}.crsspctrm,[2 4 3 1]); end                    
                    if isfield(ft_data{f},'cohspctrm'), data.(outtype)(f).coh   = permute(ft_data{f}.cohspctrm,[2 4 3 1]); end     
                    if isfield(ft_data{f},'plvspctrm'), data.(outtype)(f).plv   = permute(ft_data{f}.plvspctrm,[2 4 3 1]); end     
                otherwise
                    mmil_error('Unknown fieldtrip dimord: %s', ft_data{f}.dimord);
            end;
        otherwise
            mmil_error('Unknown data type: %s', outtype);
    end;
    %%%%%%%%%%
    % dealing with channels
    [a,badchan_idx] =setdiff({data.sensor_info.label}, ft_data{f}.label);
    [a,goodchan_idx]=intersect({data.sensor_info.label}, ft_data{f}.label);
    badchan_idx = sort(badchan_idx);
    goodchan_idx = sort(goodchan_idx);
    % Channels with no data are marked as 'bad'
    [data.sensor_info(badchan_idx).badchan]   = deal(1);
    if     isfield(data.(outtype)(f),'data'), datfield = 'data';
    elseif isfield(data.(outtype)(f),'power'),datfield = 'power';
    elseif isfield(data.(outtype)(f),'cmplx'),datfield = 'cmplx';
    end
    % Data must be expanded out to include all channels
    if ndims(data.(outtype)(f).(datfield)) == 2
        if isfield(data.(outtype)(f),'data')
            data.(outtype)(f).data(goodchan_idx,:) = data.(outtype)(f).data;
            data.(outtype)(f).data(badchan_idx,:)  = nan;
        end
        if isfield(data.(outtype)(f),'stdev')
            data.(outtype)(f).stdev(goodchan_idx,:) = data.(outtype)(f).stdev;
            data.(outtype)(f).stdev(badchan_idx,:) = 0;
        end
    elseif ndims(data.(outtype)(f).(datfield)) == 3
        if isfield(data.(outtype)(f),'data')
          data.(outtype)(f).data(goodchan_idx,:,:) = data.(outtype)(f).data;
          data.(outtype)(f).data(badchan_idx,:,:)  = nan;
        end
        if isfield(data.(outtype)(f),'power')
            data.(outtype)(f).power(goodchan_idx,:,:) = data.(outtype)(f).power;
            data.(outtype)(f).power(badchan_idx,:,:)  = nan;
        end
        if isfield(data.(outtype)(f),'cmplx')
            data.(outtype)(f).cmplx(goodchan_idx,:,:) = data.(outtype)(f).cmplx;
            data.(outtype)(f).cmplx(badchan_idx,:,:)  = nan;
        end
    elseif ndims(data.(outtype)(f).(datfield)) == 4
        if isfield(data.(outtype)(f),'data')      
          data.(outtype)(f).data(goodchan_idx,:,:,:) = data.(outtype)(f).data;
          data.(outtype)(f).data(badchan_idx,:,:,:)  = nan;
        end
        if isfield(data.(outtype)(f),'power')
            data.(outtype)(f).power(goodchan_idx,:,:,:) = data.(outtype)(f).power;
            data.(outtype)(f).power(badchan_idx,:,:,:)  = nan;
        end
        if isfield(data.(outtype)(f),'cmplx')
            data.(outtype)(f).cmplx(goodchan_idx,:,:,:) = data.(outtype)(f).cmplx;
            data.(outtype)(f).cmplx(badchan_idx,:,:,:)  = nan;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validate the result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errors = ts_checkdata(data);

if (length(errors)~=0)
    mmil_error('Errors in converting object: %s.', sprintf('\t%s\n', errors{:}));
end;

  
