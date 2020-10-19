function filename = ts_create_filename(fun,varargin)
% {condition} {value}
% eval(condition) must be logical
% eval(value) must be a string to be appended to the filename
%
% other notes:
% use the word "this" to refer to the current parameter
% use parms.PARAMETER to refer to the parameter "PARAMETER"
% use .. instead of commas
% use " instead of '
%
% example:
% to append the cutoff freq (lpfreq) of a LPF if LPF is set to 'yes',
% enter: {strcmp(this.."yes")} {parms.lpfreq}
%
% to append "toi" and time limits,
% enter: {~isempty(this)} {["toi" parms.toi(1) "-" parms.toi(end)]}

inparms = mmil_args2parms(varargin,...
						{'verbose',0,{0,1},...
             'filename',[],[],...
             'filenamefile',[],[],...
             'prefix','proc',[],...
             'rootoutdir',pwd,[],...
             'outpath','',[],...
						},false);

% return if filename already exists
if ~isempty(inparms.filename)
  filename = inparms.filename;
  return;
end

% get file with filename specifications
if isempty(inparms.filenamefile) || ~exist(inparms.filenamefile,'file')
  inparms.filenamefile = find_masterfile('filenames');
end

% load filename specifications
delimiter = ',';
comment   = '';
quotes    = '';
options   = 'empty2NaN';
[data, result] = mmil_readtext(inparms.filenamefile, delimiter, comment, quotes, options);

nanMask = result.emptyMask;
numMask = result.numberMask;
strMask = result.stringMask;

keep = strMask(:,1) & strMask(:,2) & ~nanMask(:,3);
data = data(keep,1:3);
data = data(strmatch(fun,data(:,1)),:);

% replace "" with "
for r = 1:size(data,1)
  try
%     data{r,3} = strrep(data{r,3},'""','''');
  end
end

data = parsecell(data);

% string replacements
dlim = {'..',delimiter};  % replace "__" with delimiter
this = 'parms.';          % replace "this" with this

keep = [];
% convert statements into valid matlab expressions
for r = 1:size(data,1)
  % eliminate parms with spaces in col 1 or 2
  if any(regexp(data{r,1},'\W')) || any(regexp(data{r,2},'\W'))
    continue;
  end
  cond = data{r,3};
  if ~iscell(cond), cond = {cond}; end
  % loop over elements of the cell in col 3
  for k = 1:length(cond);
    % replace "" with "
    if any(strfind(cond{k},'""'))
      cond{k} = strrep(cond{k},'""','''');
    end
    % replace __ with ,
    if any(strfind(cond{k},dlim{1}))
      cond{k} = strrep(cond{k},dlim{1},dlim{2});
    end
    % replace "this" with parms.PARAMETER
    if any(strfind(cond{k},'this'))
      cond{k} = strrep(cond{k},'this',['parms.' data{r,2}]);
    end
  end
  data{r,3} = cond;
  clear cond;
  keep = [keep r];
end
data = data(keep,:);

% parse varargin
pname = {};
parms = [];
for i = 1:length(varargin)
  if ischar(varargin{i}) && mod(i,2)
    pname = {pname{:} varargin{i}};
    parms.(pname{end}) = varargin{i+1};
  end
end

filename = inparms.prefix;
if isstruct(parms)
  parmlist = fieldnames(parms);
else
  parmlist = {};
end
% % loop over conditional statements and construct filename
for r = 1:size(data,1)
  if length(data{r,3}) == 1
    cond = {};
    vals = data{r,3}{1};
  elseif length(data{r,3}) == 2
    cond = data{r,3}{1};
    vals = data{r,3}{2};    
  else
    continue;
  end
  if ~isempty(cond)
    try
      if ~eval(cond)
        continue;
      end
    catch
      continue;
    end
  end
  try
    str = [eval(vals)];
    if ~isempty(str)
      filename = [filename '_' str];
    end
  end
end
filename = [filename '.mat'];
filename = fullfile(inparms.rootoutdir,inparms.outpath,filename);
