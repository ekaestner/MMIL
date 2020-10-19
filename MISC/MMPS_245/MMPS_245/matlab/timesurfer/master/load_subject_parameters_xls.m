function subj = load_subject_parameters(opt,varargin)
% Purpose: creates an options structure for each subject within which
% subject-specific parameters are contained in a global field.
% Notes:
% - subject-specific parameters are "derived" from other parameters
% - derived parameters = f( value parameters, base parameters )
% - opt
%     global    -- can differ b/w subjects
%     function  -- the same for all subjects
% Example:
%   derived parameter: rootoutdir = rootdir/subjectid{s}/inpath
%   base parameters:   rootdir, subjectid (defined in spreadsheet)
%   value parameter:   inpath (the value assigned to 'inpath' in the input opt structure)
%
% Created by Jason Sherfey on 18-Dec-2008


 
parms = mmil_args2parms(varargin,{...
        'subjectfile',opt.global.subjectfile,[],...
        },false);

if strfind(parms.subjectfile,'.xls')
  dat = ts_read_excel_setup(parms.subjectfile);
elseif strfind(parms.subjectfile,'.csv')
  dat = ts_csv2funparms(parms.subjectfile);
else
  error('unrecognized subject file format.\n');
end

% base, derived, and value parameters
try
    deriv = dat.derived;  fd = fieldnames(deriv); 
catch
    deriv = [];           fd = [];
end
try
    value = dat.value;    fv = fieldnames(value);
catch
    value = [];           fv = [];
end
try
    base = dat.base;
    if isfield(base,'fun_delimiter')
      delimiter = base.fun_delimiter;
      base = rmfield(base,'fun_delimiter');
    else
      delimiter = '__'; % to specify params for individual functions
    end
    fb = fieldnames(base);
catch
    base = []; fb =[]; delimiter = '__';
end

% evaluate value fields
for k = 1:length(fv)
  try
    eval(sprintf('value.(fv{k})=%s;',value.(fv{k})));
    if iscell(value.(fv{k})), value.(fv{k}) = []; end
  catch
    value = rmfield(value,fv{k});
  end
end

% subjects
subjectid    = base.subjectid; 
nsub         = length(subjectid);
subj(1:nsub) = opt;
funs         = {subj(1).function.name};

if ~iscell(subjectid), subjectid = {subjectid}; end

for s = 1:nsub
  base.subjectid = subjectid{s};
  subj(s).subjectid = subjectid{s};
  der = deriv;
  for i = 1:length(fd)    % derived fields
    for j = 1:length(fb)  % base fields
      % insert base fields into derived field
      der.(fd{i}) = strrep(der.(fd{i}),fb{j},base.(fb{j}));
    end
    for j = 1:length(fv)  % value fields
      % insert value fields into derived field
      der.(fd{i}) = strrep(der.(fd{i}),fv{j},value.(fv{j}));
    end
    % insert derived field values into this derived field
    der.(fd{i}) = check_derivations(der.(fd{i}),der);
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
      subj(s).function(funidx).opt.(funfld) = der.(fd{i});
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
  parmarray(i).opt.(fld) = val;
end

function str = check_derivations(str, derived)
f = fieldnames(derived);
for i = 1:length(f)
  if iscell(derived.(f{i})), continue; end
  expr = strcat('derived.',f{i});
  str = strrep(str,expr,derived.(f{i}));
end