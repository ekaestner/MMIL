function [parms cflag,rawdata] = ts_read_csv_parms(fname,varargin)
% Purpose: read function parameters from csv file
% Output: parms.(functions).(parameters) = values

% located in second column to right of function name in first column
% 0 - skip function
% 1 - run on local machine
% 2 - run on cluster
% 3 - auto script only

% load csv file
delimiter = ',';
comment   = '';
quotes    = '';
options   = 'empty2NaN';
[data, result] = mmil_readtext(fname, delimiter, comment, quotes, options);

if size(data,2) < 3
  data(:,3)               = deal({nan});
  result.emptyMask(:,3)   = ones(size(data,1),1);
  result.numberMask(:,3)  = zeros(size(data,1),1);
  result.stringMask(:,3)  = zeros(size(data,1),1);
end

nanMask = result.emptyMask;
numMask = result.numberMask;
strMask = result.stringMask;

% split cflags from parms
ckeep = strMask(:,1) & numMask(:,2);
pkeep = strMask(:,1) & strMask(:,2) & ~nanMask(:,3);

cflag = data(ckeep,:);
pdata = data(pkeep,1:3);

% get functions with cflags set to 0
ckeep = [cflag{:,2}] ~= 0;
rmfun = cflag(~ckeep,1);
cflag = cflag(ckeep,:);

% eliminate parms with spaces in col 1 or 2
pkeep = [];
for r = 1:size(pdata,1)
  if ~any(regexp(pdata{r,1},'\W')) && ~any(regexp(pdata{r,2},'\W'))
    pkeep = [pkeep r];
  end
end
pdata = pdata(pkeep,:);

% eliminate functions with cflags set to 0
for f = 1:length(rmfun)
  pdata = pdata(~ismember(pdata(:,1),rmfun{f}),:);
end

% convert strings to cells, arrays, etc.
pdata = parsecell(pdata,'all');  

% make parms structure
for r = 1:size(pdata,1)
  parms.(pdata{r,1}).(pdata{r,2}) = pdata{r,3};
end

cfuns = cflag(:,1);
if exist('parms','var')
  pfuns = fieldnames(parms);
else
  pfuns = {};
end

% note: setdiff(X,Y) => returns values in X that are not in Y

% given cflag but no parms => set parms.function = function name
temps = setdiff(cfuns,pfuns);
if ~isempty(temps)
  for f = 1:length(temps)
    parms.(temps{f}).function = temps{f};
  end
end
clear temps
  
% given parms but no cflag
temps = setdiff(pfuns,{cfuns{:} 'global'});
if ~isempty(temps)
  for f = 1:length(temps)
    cflag{end+1,1} = temps{f};
    cflag{end,2} = 1;
    cflag(end,3:end) = deal({nan});
  end
end
clear temps

% order parms according to cflags
if isfield(parms,'global')
  parms = orderfields(parms,{'global' cflag{:,1}});
else
  parms = orderfields(parms,cflag(:,1));
end
