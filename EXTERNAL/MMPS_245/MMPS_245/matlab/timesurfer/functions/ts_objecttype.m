function objecttype = ts_objecttype(data)
%
%
% Purpose:
%   Determines which type of timesurfer object the input object is,
%     validating it in the process.
%
% Input Args:
%   data - a timesurfer data structure
%
% Output Args:
%   objecttype: a string denoting the timesurfer object type:
%               'averages'  - average_data structure
%               'epochs'    - epoch_data structure
%               'timefreq'  - timefreq_data structure
%               'stats'     - stat_data structure
%               ''         - an unknown data structure
%

  % Attempt to validate the object
  errors = ts_checkdata(data);
  
  % Determine object type based on data field
  if (isfield(data,'epochs'))
    objecttype = 'epochs';
  elseif (isfield(data,'epoch'))
		objecttype = 'epoch';
  elseif (isfield(data,'averages'))
    objecttype = 'averages';
  elseif (isfield(data,'average'))
		objecttype = 'average';
  elseif (isfield(data,'timefreq'))
    objecttype = 'timefreq';
  elseif (isfield(data,'stat'))
		objecttype = 'stat';
  elseif (isfield(data,'stats'))
    objecttype = 'stats';
  else
    error('%s: Unknown input datatype'); %objecttype = '';
  end;
  
  % Error for invalid KNOWN objects
  if (~isempty(objecttype) && ~isempty(errors))
    error('%s: Corrupt known %s object:\n%s', mfilename, objecttype, ...
          sprintf('\t%s\n', errors{:}));
  end;
    
  % At this point we either know the object from top to bottom, or it's
  % 'unknown'
