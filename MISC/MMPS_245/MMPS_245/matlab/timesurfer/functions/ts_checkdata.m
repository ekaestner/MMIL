function errors = ts_checkdata(data, objtype, errtype)
%function errors = ts_checkdata(data, [type]);
%
% Purpose:
%   Verifies data satifies the requirements to be a valid TimeSurfer epoch
%   or average
%   data structure.
%
% Input Args:
%   type:    what data type to verify
% Output Args:
%
% Created:  08/29/07 by Rajan Patel
% Last Mod: 03/23/11 by Don Hagler
%

% Revision: 10/05/2008 by Jason Sherfey
%   Commented out code to validate # of fields

if (~exist('objtype','var'))
  objtype = 'default';
end;
if (~exist('errtype','var'))
  errtype = {'warning','error'};
end;

% 
errors = {};
fields = fieldnames(data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These input types are known object fields/types
switch (objtype)

  % generic to avg, epoch, timefreq
  case 'default'

    % Check known fields to make sure sub-fields are
    % set as appropriate
    for i=1:length(fields)

      cur_field=fields{i};

      switch cur_field

        case 'num_sensors'
          %% To do: Check if num_sensors is correct or agrees with
          %% size of one of the epoch.data dimensions

        case 'sfreq'
          %% To do: Check if sfreq is an appropriate rate.
        case {'obj','cfg','parms'}
          %% 

        % sub-objects get evaluated separately
        case {'epochs','averages','timefreq',...
              'coor_trans', 'noise'}

          errors{end+1} = ts_checkdata(data.(cur_field), cur_field,  errtype);

      end;
    end;

  %%%%%%%%%%%%%%%%%%%%

  case {'averages', 'timefreq', 'epochs'}

%      switch(objtype)
%        case 'averages'
%          if (length(fields)~=6)
%            errors{end+1} = 'Inappropriate number of fields for averages.';
%          end;
%        case 'epochs'
%          if (length(fields)~=5)
%            errors{end+1} = 'Inappropriate number of fields for epochs.';
%          end;
%        case 'timefreq'
%          if (length(fields)<5)
%            errors{end+1} = 'Inappropriate number of fields for timefreq.';
%          end;
%      end;

    for j=1:length(fields)
      cur_field=fields{j};

      switch cur_field
        case 'event_code'

        case 'num_trials'
          %% To do:  Check if matches a dimension of data.

        case 'time'
          %%  To do: Check if time is appropriate

        case 'num_rejects'
          errors{end+1} = ts_checkdata([data.(cur_field)], cur_field, errtype);

        case 'data'
          %%  To do: Check if size of data is appropriate

        case 'stdev'
          if (~strcmp(objtype,'averages'))
            errors{end+1} = 'stdev can only be defined on average data';
          end;

        % Modified by Jason (added 'cmplx')  
        case {'frequencies', 'coeff', 'power','cmplx'}
          if (~strcmp(objtype,'timefreq') && ~strcmp(objtype,'averages'))
            errors{end+1} = 'stdev can only be defined on timefreq data';
          end;
        otherwise
        % Modified by Jason (comment otherwise errors)    
          %errors{end+1} = sprintf('%s contains unknown field: %s',objtype, cur_field);
      end
    end

  %%%%%%%%%%%%%%%%%%%%

  case 'epochs'

    if (length(fields)~=5)
      errors{end+1} = 'Inappropriate number of fields for epochs.';
    end;

    for j=1:length(fields)
      cur_field=fields{j};

      switch cur_field
        case 'event_code'

        case 'num_trials'
          %% To do:  Check if matches a dimension of data.

        case 'time'
          %%  To do: Check if time is appropriate

        case 'data'
          %%  To do: Check if size of data is appropriate

        case 'num_rejects'
          errors{end+1} = ts_checkdata(data.(cur_field), cur_field,  errtype);

        otherwise
        % Modified by Jason (comment otherwise errors)    
          %errors{end+1} = sprintf('%s contains unknown field: %s',objtype,cur_field);
      end
    end

  %%%%%%%%%%%%%%%%%%%%

  case 'noise'
    if length(fields)~=3
      errors{end+1} = 'Inappropriate number of fields for noise.';
    end;

    for j=1:length(fields)
      cur_field=fields{j};
      switch cur_field
        case {'num_trials','num_samples','covar'}
        otherwise
        % Modified by Jason (comment otherwise errors)    
          %errors{end+1} = sprintf('Unknown field for noise: %s', cur_field);
      end
    end

  %%%%%%%%%%%%%%%%%%%%

  case 'sensor_info'

    if length(fields)==7
      errors{end+1} = 'Inappropriate number of fields for sensor_info.';
    end;

    for j=1:length(fields)
      cur_field=fields{j};
      switch cur_field
        case 'type'
        case 'label'
          %%  To do:  Check if number of labels and
          %%  num_sensors agree?
        case 'loc'
          %%  To do:  Check size of matrix
        case 'badchan'
          %%  To do:  Make sure either 0 or 1
        case 'typestring'  
          for k=1:length(data.sensor_info)
            cur_ts = data.sensor_info(k).typestring;
            switch cur_ts
              case {'eeg','eog','grad','grad1','grad2','mag','sti'}
              otherwise
                errors{end+1} = sprintf('Sensor %d is an unknown sensor type: %s',k,cur_ts);
            end
          end
        case 'kind'
        case 'lognum'
        otherwise
          errors{end+1} = sprintf('sensor_info contains unknown field: %s',cur_field);
      end
    end

  %%%%%%%%%%%%%%%%%%%%

  case 'coor_trans'

    if length(fields)~=2
      errors{end+1} = 'Inappropriate number of fields for coor_trans.';
    end;

    for j=1:length(fields)
      cur_field=fields{j};
      switch cur_field
        case 'device2head'
          %% To do:  Check matrix.
        case 'mri2head'
          %% To do:  Check matrix.
        otherwise
          errors{end+1} = sprintf('coor_trans contains unknown field: %s',cur_field);
        end
    end

  %%%%%%%%%%%%%%%%%%%%

  case 'num_rejects'

    if length(fields)~=6
      errors{end+1} = 'Inappropriate number of fields for num_rejects.'; 
    end;

    for k=1:length(data)

      for l=1:length(fields)
        cur_field=fields{l};
        switch cur_field
          case {'mag','grad','eog','eeg','manual','skip'}
            %if (data.(k).num_rejects.(cur_rejfield) >= data.(k).num_trials)
            %    display('Warning: Number of rejects is equal to or exceeds num_trials!');
            %end
          otherwise
            errors{end+1} = sprintf('Unknown field for num_rejects: %s.', cur_field);
        end
      end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finally, remove all empty 'errors' by flattening the cell array
tmp = {};
for i=1:length(errors)
    if (~iscell(errors{i}))
        tmp = {tmp{:} errors{i}};
    else
        tmp = {tmp{:} errors{i}{:}};
    end;
end;
errors = tmp;

for i=length(errors):-1:1
  if (mmil_isempty(errors{i},true))
    errors = errors(1:i-1, i+1:end);
  end;
end;

