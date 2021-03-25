function mmil_write_script(fname,cmd,args,tags,parms,quitflag,appendflag,...
  header)
%function mmil_write_script(fname,cmd,args,tags,parms,quitflag,appendflag,...
%  header)
%
% Purpose: Create matlab script for batch processing
%   Handles multiple data types, including vectors and cell arrays
%
% Required Input:
%   fname: output file name
%   cmd: matlab function name
%   args: cell array of arguments to come before ('key',value) pairs
%   tags: cell array of names of fields in parms (the 'key's)
%   parms: struct containing parameter values (the values)
%
% Optional Input:
%   quitflag: [0|1] whether to add 'exit' statement at end of script
%     {default = 1}
%   appendflag: [0|1] whether to append to existing file
%     {default = 0}
%   header: string containing line or lines to appear at beginning of script
%     {default = []}
%
% Created:  05/17/09 by Don Hagler
% Last Mod: 03/29/11 by Don Hagler
%

if (~mmil_check_nargs(nargin,5)) return; end;
if ~exist('quitflag','var') | isempty(quitflag), quitflag = 1; end;
if ~exist('appendflag','var') | isempty(appendflag), appendflag = 0; end;
if ~exist('header','var'), header = []; end;

if appendflag
  fid = fopen(fname,'a');
else
  fid = fopen(fname,'w');
end;
if fid==-1
  error('failed to open file %s for writing',fname);
end;
fprintf(fid,'%s(...\n',cmd);

if ~isempty(header)
  fprintf(fid,'%s\n',header);
end;

if ~iscell(args), args = {args}; end;
for a=1:length(args)
  mmil_write_value(fid,args{a});
  if a<length(args) || length(tags)~=0
    fprintf(fid,',');
  end;
  fprintf(fid,'...\n');
end;

mmil_write_tags(fid,tags,parms);

fprintf(fid,');\n');
if quitflag, fprintf(fid,'exit\n'); end;
fclose(fid);

