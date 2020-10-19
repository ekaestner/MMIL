function mmil_error( parms, varargin )
%function mmil_error( parms, varargin )
%
% Purpose: 
%   Abstract how to generate and handle errors in MMIL code.
%
%   Currently, this is done via the 'error' function.  We
%   also take care to log the error to any log file.
%
% Parameters:
%  parms: mmil optional argument 'parms' object
%
%   
%  Created By:       Ben Cipollini on 08/20/2007
%  Last Modified By: Ben Cipollini on 08/20/2007

  mmil_check_nargs( nargin, 1 );
  
  % Log to any external log stream
  if ((isfield(parms, 'logfile') && ~isempty(parms.logfile)) ...
    || (isfield(parms, 'logfid') && parms.logfid ~= 1))
    mmil_logstr( parms, varargin{:} );
  end;
  
  all_callers = dbstack;
  caller      = all_callers(2).name;
  
  dbstack;
  error( ['%s: ' varargin{1}], caller, varargin{2:end} );
