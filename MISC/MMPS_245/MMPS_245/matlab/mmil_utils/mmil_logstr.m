function mmil_logstr(varargin)
%function mmil_logstr(varargin)
%
% Purpose:
%   string logging with consistent formatting
%     always output mfilename and newline
%
% Input Parameters:
%   parms    : (optional) input parameter object with above fields defined
%     logfile: file where you want to log to 
%       OR       
%     logfid: file ID where to log to (default: 1, or stdout)
%     verbose' [0|1]: whether to output or not
%   str      : string to "log"
%   varargin : arguments to "print" into the format string str
%
% Example parameter parsing:
%
%  parms   = mmil_args2parms( varargin, ...
%                            { 'verbose',        true,     [false true], ...
%                              'logfile',        [],       [],...
%                              'logfid',         1,        [] 
%                            });
%
% See also: fprintf
%
% Created:  07/15/2007 by Ben Cipollini
% Prev Mod: 08/15/2007 by Ben Cipollini
% Last Mod: 09/08/2016 by Don Hagler
%

if ~mmil_check_nargs(nargin, 1), return; end;

logfid = 1;

if ~isstruct(varargin{1})
  % no parms specified
  logfid = 1;
else
  % parms object tells us whether to print or not
  mmil_check_nargs(nargin, 2, 'Not enough input args; must specify a string to log.');
  parms = varargin{1};
  % should not log anything
  if (isfield(parms,'verbose') && ~parms.verbose)
    return;
  end;
  if (isfield(parms,'logfile') && ~isempty(parms.logfile)) &&...
     (~isfield(parms,'logfid') || parms.logfid==1)
    logfid = fopen(parms.logfile, 'a');
    we_opened_logfid = true;
  elseif isfield(parms,'logfid')
    logfid = parms.logfid;
  else
    logfid = 1;
  end;
  % shift args so output string is 1
  varargin = varargin(2:end);
end;

% if here, should print
M = dbstack;
  
if (length(M)<2)
  caller = '';
else
  caller = M(2).name;
  if (strcmp(caller, 'logstr')==1 ||strcmp(caller,'mmil_error')==1)
    if (length(M) < 3)
      caller = '';
    else
      caller = M(3).name;
    end;
  end;
end;

if (strcmp(caller, '')==1)
  prefix = '';
else
  prefix = sprintf('%s: ', caller);
end;

% prep the args
formatstr = varargin{1};
if (length(varargin) > 1)
  args = {varargin{2:end}};
else
  args = {};
end;

% output the info
outstr = sprintf(formatstr, args{:});
fprintf(logfid, [prefix outstr '\n']);

% if necessary, close the logfid (e.g. if we opened it)
if (exist('we_opened_logfid', 'var') && we_opened_logfid)
  fclose(logfid);
end;

