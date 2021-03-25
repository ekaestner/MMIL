function [status,result]=mmil_unix(cmd,batch_limit)
%function [status,result]=mmil_unix(cmd,[batch_limit])
%
% Purpose: execute a unix command
%   moves MATLAB paths to bottom of LD_LIBRARY_PATH stack
%   to simulate non-MATLAB shell environment.
%
% Required Input:
%   cmd : UNIX command-line to execute
%
% Optional Input:
%   batch_limit: number of lines to execute at a time
%     {default = 250}
%
% Output:
%   status: 0 if successful, 1 if error
%   result: stdout/stderr output from running cmd
%
% See also: unix
%
% Created:  08/01/07 by Ben Cipollini (as smartunix)
% Last Mod: 08/05/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('batch_limit','var'), batch_limit = 250; end;
if isempty(cmd), error('cmd is empty'); end;

% remove MATLAB directories from LD_LIBRARY_PATH
path = getenv('LD_LIBRARY_PATH');
paths = regexp(path,':','split');
newpaths = paths(cellfun(@isempty,regexp(paths,'matlab')));
if isempty(newpaths)
  newpath = '';
else
  newpath = [sprintf('%s:',newpaths{1:end-1}) newpaths{end}];
end;
setenv('LD_LIBRARY_PATH',newpath);

% split command into multiple lines
cmd_list = regexp(cmd,'\n','split');
ncmds = length(cmd_list);
if ncmds<=1
  [status,result] = unix(cmd);
else
  status = 0;
  result = [];
  nbatches = ceil(ncmds/batch_limit);
  for b=1:nbatches
    j = 1 + batch_limit*(b-1);
    k = min(j + batch_limit - 1,ncmds);
    tmp_cmds = cmd_list(j:k);
    tmp_cmd = sprintf('%s\n',tmp_cmds{:});
    [status,tmp_result] = unix(tmp_cmd);
    result = [result tmp_result];
    if status, return; end;
  end;
end;

% restore LD_LIBRARY_PATH
setenv('LD_LIBRARY_PATH',path);

