function opt = ts_make_output_filename(opt,varargin)
% field (optional): cell array of strings listing option structures to use
% name (default: 'all')
%
% Last Mod: 09/15/12 by Don Hagler
%

% do nothing if filename already exists
if isfield(opt,'filename') && ~isempty(opt.filename), return; end
if isfield(opt,'outfile') && ~isempty(opt.outfile), 
    opt.filename = opt.outfile; 
end

% set up prefix
if ~isfield(opt,'prefix') || isempty(opt.prefix) 
  if issubfield(opt,'previous.filename')
    if iscell(opt.previous.filename), fname = opt.previous.filename{1};
    else fname = opt.previous.filename; end
    [pathstr,name,ext] = fileparts(fname);
    [jnk,name] = fileparts(fullfile(pathstr,name));  % remove id preceding extension
    opt.prefix = [name];
  else
    opt.prefix = 'proc';
  end
end
if ~isfield(opt,'filenamefile') || ~exist(opt.filenamefile,'file')
  [opt.filenamefile ext] = find_masterfile('filenames');
end
if isfield(opt,'filenamefile') && exist(opt.filenamefile,'file')
  [raw, result] = mmil_readtext(opt.filenamefile,',','','','empty2NaN');
else
  % read filename mappings from master control file
  if ~isfield(opt,'funspecfile') || ~exist(opt.funspecfile,'file')
      [filename,pathname] = uigetfile('.xls','Master control file not found.  Pick an excel spreadsheet');
      spreadsheet = fullfile(pathname,filename);
  else
      spreadsheet = opt.funspecfile;
  end
  if ~exist(spreadsheet,'file')
      error('%s: Master control file with output_filenames worksheet not found.',mfilename);
  end
  warning('off','MATLAB:xlsread:Mode')
  [num,txt,raw]=xlsread(spreadsheet,'output_filenames');
  warning('on','MATLAB:xlsread:Mode')
end

% keep rows for this function
row  = strmatch(opt.function,raw(:,1),'exact');
col1 = raw(row,1);
col2 = raw(row,2);
col3 = raw(row,3);
col4 = raw(row,4);

% keep rows with yes flags
row  = find(strcmp(col3,'yes'));
fld = col2(row);
col4 = col4(row);

% keep rows that are in the opt structure
[com row jnk] = intersect(fld,fieldnames(opt));
row = sort(row);
fld = fld(row);
str = col4(row);

% replace 'value' with parameter value
id = strmatch('value',str);
for i=1:length(id)
  % check for prefix string after 'value'
  tmp = str{id(i)};
  if length(tmp) > 5  % more than value
    strval = tmp(6:end);
  else
    strval = '';
  end
  this = fld{id(i)};
  val = opt.(this);
  if isempty(val), strval = '';
  elseif isstr(val), strval = strcat(strval,val);
  elseif isnumeric(val)
    if length(val)==1, strval = strcat(strval,num2str(val));
    else 
%      for j = 1:length(val)
        strval = sprintf('%s%g-%g',strval,val(1),val(end));
%      end
    end
  else
    strval = strcat(strval,'');
  end
  str{id(i)} = strval;
end
str = str(~strcmp('',str));
str = strcat(repmat('.',size(str,1),1),str);
str = [str{:}];

% savepath
if ~isfield(opt,'savepath') || isempty(opt.savepath)
  if isfield(opt,'rootoutdir')
    opt.savepath = opt.rootoutdir;
  elseif issubfield(opt,'previous.rootoutdir')
    opt.savepath = opt.previous.rootoutdir;
  else
    opt.savepath = pwd;
  end
  try opt.savepath = fullfile(opt.savepath,opt.outpath); end
end
opt.filename = fullfile(opt.savepath,[opt.prefix str '.mat']);
