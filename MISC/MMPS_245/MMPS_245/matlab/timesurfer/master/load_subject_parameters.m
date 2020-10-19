function subj = load_subject_parameters(opt,varargin)
% Purpose: creates an options structure for each subject within which
% subject-specific parameters are contained in a global field.
% Notes:
% - subject-specific parameters are "derived" from other parameters
% - derived parameters = f( parmfile parameters, constant parameters )
% - opt
%     global    -- can differ b/w subjects
%     function  -- the same for all subjects
% Example:
%   derived parameter: rootoutdir = rootdir/subjectid{s}/inpath
%   constant parameters:   rootdir, subjectid (defined in spreadsheet)
%   parmfile parameter:   inpath (the value assigned to 'inpath' in the input opt structure)
%
% Created by Jason Sherfey on 18-Dec-2008

parms = mmil_args2parms(varargin,{...
        'subjectfile',opt.global.subjectfile,[],...
        },false);

if ~exist(parms.subjectfile,'file')
  subj = opt;
  return;
end

if strfind(parms.subjectfile,'.xls')
  subj = load_subject_parameters_xls(opt,varargin{:});
  return;
elseif strfind(parms.subjectfile,'.csv')
%   dat = ts_csv2funparms(parms.subjectfile);
else
  error('unrecognized subject file format.\n');
end
% 
% delimiter = ',';
% comment   = '';
% quotes    = '';
% options   = 'empty2NaN';
% [data, result] = mmil_readtext(fname, delimiter, comment, quotes, options);
% data = parsecell(data);

% read subjectfile
dat = ts_csv2funparms(parms.subjectfile);

% backward compatibility
if isfield(dat,'base') && ~isfield(dat,'constant')
  dat.constant = dat.base;
  dat = rmfield(dat,'base');
end
if isfield(dat,'value') && ~isfield(dat,'parmfile')
  dat.parmfile = dat.value;
  dat = rmfield(dat,'value');
elseif isfield(dat,'parms') && ~isfield(dat,'parmfile')
  dat.parmfile = dat.parms;
  dat = rmfield(dat,'parms');
end

% set up constant, derived, and parmfile parameters
try
    deriv = dat.derived;       fd = fieldnames(deriv); 
catch
    deriv = [];                fd = [];
end
try
    derivset=dat.derived_set;  fds = fieldnames(derivset); 
catch
    derivset=[];               fds = [];
end
try
    parmfile = dat.parmfile;   fv = fieldnames(parmfile);
catch
    parmfile = [];             fv = [];
end
try
    constant = dat.constant;
    if isfield(constant,'fun_delimiter')
      delimiter = constant.fun_delimiter;
      constant = rmfield(constant,'fun_delimiter');
    elseif isfield(constant,'delimiter')
      delimiter = constant.delimiter;
      constant = rmfield(constant,'delimiter');
    else
      delimiter = '__'; % to specify params for individual functions
    end
    fb = fieldnames(constant);
catch
    constant = []; fb =[]; delimiter = '__';
end

% evaluate parmfile fields
for k = 1:length(fv)
  try
    eval(sprintf('parmfile.(fv{k})=opt.global.%s;',parmfile.(fv{k})));
    if iscell(parmfile.(fv{k})), parmfile.(fv{k}) = []; end
  catch
    parmfile = rmfield(parmfile,fv{k});
  end
end

% obtain key name (key = cell array of strings over which to loop)
% note: this function was written to loop over subjects; therefore the
% default key is "subjectid"
toggle = [];
if isfield(dat,'key')  % key.(keyid) => cell array = # of subjects
  key = dat.key;
  if isfield(key,'toggle')
    toggle = key.toggle;
    key    = rmfield(key,'toggle');
  end
  tmp = fieldnames(key);
  keyparm = tmp{1};
elseif isfield(constant,'subjectid')
  key.subjectid = constant.subjectid;
  keyparm = 'subjectid';
  constant = rmfield(constant,'subjectid');
end

% add keyparm to "constant" fieldnames
if   iscell(fb), fb{end+1} = keyparm;
else fb = keyparm;
end
keyid = key.(keyparm); if ~iscell(keyid), keyid = {keyid}; end
nsub  = length(keyid);

% parse cell array of subject-specific values
if ~isempty(derivset)
  derivset = parse_derivset(derivset,fds,delimiter,nsub);
  keep = ones(1,length(fds));
  for k = 1:length(fds)
    if length(derivset.(fds{k})) ~= nsub
      keep(k)  = 0;
    end
  end
  derivset = rmfield(derivset,fds(keep==0));
  fds      = fds(keep==1);
end

% toggle subjects if appropriate
if ~isempty(toggle) && length(toggle) == nsub
  for k = 1:length(fds)
    derivset.(fds{k}) = derivset.(fds{k})(toggle==1);
  end
  keyid = keyid(toggle==1);    
  nsub  = length(keyid);
end

% initialize structure array of parameters for each subject
subj(1:nsub) = opt;
funs         = {subj(1).function.name};

% loop over subjects
for s = 1:nsub
  constant.(keyparm) = keyid{s};
  subj(s).(keyparm) = keyid{s};
  der = deriv;
  % add derived sets to deriv and fd
  for i = 1:length(fds)
    if length(derivset.(fds{i})) == nsub
      der.(fds{i}) = derivset.(fds{i}){s};
    else
      der.(fds{i}) = derivset.(fds{i});
    end
    if ~ismember(fds{i},fd)
      fd = {fds{i} fd{:}}';
%       fd{end+1} = fds{i};
    end
  end
  for i = 1:length(fd)    % derived fields
    for j = 1:length(fb)  % constant fields
      % insert constant fields into derived field
      if ~ischar(constant.(fb{j})), continue; end
      if ischar(der.(fd{i}))
        der.(fd{i}) = strrep(der.(fd{i}),fb{j},constant.(fb{j}));
      elseif iscell(der.(fd{i})) && ischar(der.(fd{i}){1})
        dtmp = der.(fd{i}); % this derived cell array
        for k = 1:length(der.(fd{i}))
          dtmp{k} = strrep(dtmp{k},fb{j},constant.(fb{j}));
        end
        der.(fd{i}) = dtmp;
      end
    end
    for j = 1:length(fv)  % parmfile fields
      % insert parmfile fields into derived field
      if ~isfield(parmfile,fv{j}) || ~ischar(parmfile.(fv{j})), continue; end
      if ischar(der.(fd{i}))
        der.(fd{i}) = strrep(der.(fd{i}),fv{j},parmfile.(fv{j}));
      elseif iscell(der.(fd{i})) && ischar(der.(fd{i}){1})
        dtmp = der.(fd{i}); % this derived cell array
        for k = 1:length(der.(fd{i}))
          dtmp{k} = strrep(dtmp{k},fv{j},parmfile.(fv{j}));
        end
        der.(fd{i}) = dtmp;
      end
    end
    % insert derived field values into this derived field
    der.(fd{i}) = check_derivations(der.(fd{i}),der);
    % insert constants into updated derived sets
      for j = 1:length(fb)  % constant fields
        % insert constant fields into derived field
        if ~ischar(constant.(fb{j})), continue; end
        if ischar(der.(fd{i}))
          der.(fd{i}) = strrep(der.(fd{i}),fb{j},constant.(fb{j}));
        elseif iscell(der.(fd{i})) && ischar(der.(fd{i}){1})
          dtmp = der.(fd{i}); % this derived cell array
          for k = 1:length(der.(fd{i}))
            dtmp{k} = strrep(dtmp{k},fb{j},constant.(fb{j}));
          end
          der.(fd{i}) = dtmp;
        end
      end
    
    funidx = 0; 
    for j = 1:length(funs)
      % look for function name in this derived field
      if strmatch(funs{j},fd{i})
        funidx = j;
        funfld = fd{i}(findstr(fd{i},delimiter)+length(delimiter):end);
      end
    end
    if any(funidx)
      % this value is for a particular function
      subj(s).function(funidx).parms.(funfld) = der.(fd{i});
    else
      % this value is global
      subj(s).global.(fd{i}) = der.(fd{i});
      subj(s).function = update_function_parms(subj(s).function,fd{i},der.(fd{i}));
    end
  end
%   subj(s).global.subjectid = subjectid{s};
end

function parmarray = update_function_parms(parmarray,fld,val)
for i = 1:length(parmarray)
  parmarray(i).parms.(fld) = val;
end

function str = check_derivations(str, derived)
f = fieldnames(derived);
if ischar(str)
  for i = 1:length(f)
    if ~ischar(derived.(f{i})), continue; end
    if iscell(derived.(f{i})),  continue; end
    expr = strcat('derived.',f{i});
    if ~ischar(expr),           continue; end
    str = strrep(str,expr,derived.(f{i}));
  end
elseif iscell(str) && ischar(str{1})
  tmp = str;
  for n = 1:length(str)
    str = tmp{n};
    if ~ischar(str), continue; end
    for i = 1:length(f)
      if ~ischar(derived.(f{i})), continue; end
      if iscell(derived.(f{i})),  continue; end
      expr = strcat('derived.',f{i});
      if ~ischar(expr),           continue; end
      str = strrep(str,expr,derived.(f{i}));
    end
    tmp{n} = str;
    clear str;
  end
  str = tmp;
end

function derivset = parse_derivset(derivset,fds,delimiter,nsub)
% create cell array of cells if delimiter is found
for k = 1:length(fds)
  this = derivset.(fds{k});
  if length(this) == nsub
    continue; 
  end
  idx = find(ismember(this,delimiter));
  if length(idx)+1 ~= nsub
    continue;
  end
  n         = 1;
  for j = 1:length(idx)
    temp{j} = this(n:idx(j)-1);
    n       = idx(j)+1;
  end
  temp{j+1} = this(n:end);
  % convert to numeric array if appropriate
  for i = 1:length(temp)
    idx = regexp(temp{i},'\[*\]');
    if iscell(idx), idx=[idx{:}]; end
    if ~isempty(idx)
      n = 1;
      for j = 1:length(idx)
        ix = idx(j);
        temptemp{j} = str2num(temp{i}{1}(n:ix));
        n = ix+1;
      end
      temp{i} = temptemp;
      clear temptemp
    end
  end
  derivset.(fds{k}) = temp;
  clear temp
end

