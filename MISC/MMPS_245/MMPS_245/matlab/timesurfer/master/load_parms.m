function allparms = load_parms(varargin)
% Purpose: read csv parameter file and set up processing operations
% Output: type 'help load_parameters'

[defaultfile  ext] = find_masterfile('defaults' );
[controlfile  ext] = find_masterfile('funspecs' );
[filenamefile ext] = find_masterfile('filenames');

if nargin==1 && ischar(varargin{1})
  % csv file
  parmfile = varargin{1};
  if ~iscell(parmfile), parmfile = {parmfile}; end
  [parms cflags] = ts_read_csv_parms(parmfile{1});
elseif nargin==1 && isstruct(varargin{1})
  % parm structure 
  parms = varargin{1};          % transform into proper structure below
elseif nargin==2 && ischar(varargin{1}) && strcmpi(varargin{1},'gui_study')
  [parms cflags] = read_study_parms(varargin{2});
else
  % key/value pairs (not working)
%   parms = args2parms(varargin); % transform into proper structure below
end

% set the wrapper control files
if issubfield(parms,'global.defaultfile')
  defaultfile = parms.global.defaultfile;
else
%   parms.global.defaultfile = defaultfile;
end
if issubfield(parms,'global.funspecfile')
  controlfile = parms.global.funspecfile;
else
%   parms.global.funspecfile = controlfile;
end
if ~issubfield(parms,'global.filenamefile')
%   parms.global.filenamefile = filenamefile;
end
if ~isfield(parms,'global')
  parms.global = [];
end

% get function names
fun = fieldnames(parms);
fun = fun(~strcmp(fun,'global'));

% set defaults
defaults = ts_read_csv_parms(defaultfile);
dfields = fieldnames(defaults);
for f = 1:length(dfields)
  tmpfld = dfields{f};
  if ~isfield(parms,tmpfld),continue; end
  dparms = defaults.(tmpfld);
  fnames = fieldnames(dparms);
  for i = 1:length(fnames)
    if ~isfield(parms.(tmpfld),fnames{i})
      parms.(tmpfld).(fnames{i}) = dparms.(fnames{i});
    end
  end
end

% propagate global
if isstruct(parms.global)
  gfld = fieldnames(parms.global);
  for i = 1:length(fun)
    fld = fieldnames(parms.(fun{i}));
    % empty params
    pnull = fld(structfun(@(x) isempty(x),parms.(fun{i})));
    pnull = intersect(pnull,gfld);
    % in global but not function
    pglob = setdiff(gfld,fld);
    % parameters to propagate
    pprop = unique({pglob{:} pnull{:}});
    for j = 1:length(pprop)
      parms.(fun{i}).(pprop{j}) = parms.global.(pprop{j});
    end
  end
end
%     fnames = fieldnames(parms.global);
%     for i = 1:length(fun)
%         for j = 1:length(fnames)
%           if ~isfield(parms.(fun{i}),fnames{j}) || isempty(parms.(fun{i}).(fnames{j}))
%             parms.(fun{i}).(fnames{j}) = parms.global.(fnames{j});
%           end
%         end
%     end
% end
if issubfield(parms,'global.script_flag') && parms.global.script_flag
  autoflag = 1;
else
  autoflag = 0;
end

% make the output structure
allparms.global = parms.global;
for i = 1:length(fun)
  % function name
  allparms.function(i).name = fun{i};
  if isfield(parms.(fun{i}),'script') && ischar(parms.(fun{i}).script)
    [flags,types,cmd] = read_funspecs(parms.(fun{i}).script,parms.(fun{i}));
  else
    [flags,types,cmd] = read_funspecs(fun{i},parms.(fun{i}));
  end
  % function call
  allparms.function(i).funcall   = cmd;
  % function definition and flow control flags
  allparms.function(i).spec_flags = flags;
  if isfield(parms.(fun{i}),'save_flag') && ~isempty(parms.(fun{i}).save_flag)
    allparms.function(i).spec_flags.save_flag = parms.(fun{i}).save_flag;
  end
  % allowable data types
  allparms.function(i).itype     = types.itype;
  allparms.function(i).otype     = types.otype;
  % process flags
  if isfield(parms.(fun{i}),'script_flag') && ~isempty(parms.(fun{i}).script_flag)
    allparms.function(i).process_flags.script_flag = parms.(fun{i}).script_flag;
  else
    allparms.function(i).process_flags.script_flag = autoflag;
  end
  if isfield(parms.(fun{i}),'itype') && ischar(parms.(fun{i}).itype)
    allparms.function(i).itype   = parms.(fun{i}).itype;
    allparms.function(i).funcall = strrep(allparms.function(i).funcall,['(' types.itype],['(' parms.(fun{i}).itype]);
  end
  if isfield(parms.(fun{i}),'otype') && ischar(parms.(fun{i}).otype)
    allparms.function(i).otype = parms.(fun{i}).otype;
    allparms.function(i).funcall = strrep(allparms.function(i).funcall,[types.otype ' ='],[parms.(fun{i}).otype ' =']);
    allparms.function(i).funcall = strrep(allparms.function(i).funcall,[types.otype '='],[parms.(fun{i}).otype '=']);
  end  
    runflag = [cflags{i,2}]; 
    allparms.function(i).process_flags.run_flag     = runflag;
    allparms.function(i).process_flags.cluster_flag = 0;
  if runflag == 2
    allparms.function(i).process_flags.run_flag = 0;
    allparms.function(i).process_flags.cluster_flag = 1;
  end    
  if runflag == 3
    allparms.function(i).process_flags.run_flag     = 0;
    allparms.function(i).process_flags.cluster_flag = 0;
    allparms.function(i).process_flags.script_flag  = 1;
  end
  % functions-specific parameters
  allparms.function(i).parms = parms.(fun{i});
    if size(cflags,2) >= 3
      funid = [cflags{i,3}]; 
      if ~isempty(funid) && isnumeric(funid)
        allparms.function(i).parms.function_id = funid;
      end
    end
    if size(cflags,2) >= 4
      inpid = [cflags{i,4}];
      if ~isempty(inpid) && isnumeric(inpid)
        allparms.function(i).parms.input_id = inpid;
      end    
    end
    if isfield(allparms.function(i).parms,'funcall') && ischar(allparms.function(i).parms.funcall)
      allparms.function(i).funcall = allparms.function(i).parms.funcall;
      allparms.function(i).parms   = rmfield(allparms.function(i).parms,'funcall');
    end
    clear funid inpid
end

%% Subfunctions
function [parms cflags] = read_study_parms(STUDY)
% parse STUDY from GUI into parms for ts_session
funs  = STUDY.funlist;
cnt = 0;
for f = 1:length(funs)
  % parameters
  name = funs{f};
  if ismember('.',name)
    [jnk name] = fileparts(name);
  end
  func = STUDY.protocol{f};
  if isfield(STUDY,'controls')
    cont = STUDY.controls{f};
  else
    cont{f}.id = nan;
  end
  parm = fieldnames(func);
  for j = 1:length(parm)
    tmp = func.(parm{j});
    if ischar(tmp) && strcmpi(tmp,'nan')
      val = [];
    elseif ~iscell(tmp) && all(isnan(tmp))
      val = [];
    else
      val = tmp;
    end
    parms.(name).(parm{j}) = val;
  end
  if strcmpi(name,'global')
    continue;
  else
    cnt = cnt + 1;
  end
  % control flags
  if isfield(cont{f},'action')
    switch cont{f}.action
      case 'execute'
        val = 1;
      case 'cluster'
        val = 2;
      case 'mscript'
        val = 3;
      otherwise
        val = 1;
    end
  else
    val = 1;
  end
  cflags(cnt,1:3) = {name,val,cont{f}.id};
  if isfield(cont{f},'input')
    val = cont{f}.input;
  else
    val = nan;
  end
  cflags(cnt,4) = {val};
  % override if given explicitly
  if isfield(func,'function_id') && ~isempty(func.function_id)
    parms.(name).function_id = func.function_id;
    cflags{cnt,3} = func.function_id;
  end
  if isfield(func,'input_id') && ~isempty(func.input_id)
    parms.(name).input_id = func.input_id;
    cflags{cnt,4} = func.input_id;
  end  
  if isfield(func,'cluster_flag') && isnumeric(func.cluster_flag) && ~isempty(func.cluster_flag) && func.cluster_flag==1
    cflags{cnt,2} = 2;
  elseif isfield(func,'script_flag') && isequal(func.cluster_flag,1)
    cflags{cnt,2} = 3;
  end
end 

%% TODO: 
% 1. Allow key/value pairs
% 2. Allow multiple parameter files
%
% Starting Point:
% 
% function parms = args2parms(varargin)
% if mod(length(varargin),2)
%   error('List of arguments must be even (must have name/value pair)');
% end
% for i = 1:length(varargin)
%   if mod(i,2) && (~ischar(varargin{i}) || isempty(varargin{i}))
%     error('names must be specified as strings');
%   end
% 	if mod(i,2)
%     name = varargin{i};
%   else
%     parms.(name) = varargin{i};
%   end
% end
% 
% loop over parmfiles
% for f = 1:length(parmfile)
%   fname = parmfile{f};
%   if f==1
%     parms = ts_read_csv_parms(fname);
%   else
%     tmp = ts_read_csv_parms(fname);
%     for i = 1:numfunfields
%       if ~isfield(parms,thisfun), parms.(thisfun) = tmp.(thisfun); continue; end
%       for j = 1:numparmfields
%         if ~isfield(parms.(thisfun),thisparm)
%           parms.(thisfun).(thisparm) = tmp.(thisfun).(thisparm);
%         end
%       end
%     end
%   end
% end
