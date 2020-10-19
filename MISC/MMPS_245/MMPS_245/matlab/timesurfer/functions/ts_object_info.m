function [datatype,datafield,dataparam] = ts_object_info(data,varargin)
% [object]                       = ts_object_info(data)
% [object,datafield]             = ts_object_info(data)
% [object,datafield,dataparam]   = ts_object_info(data)
%
% Outputs:
% object    - string
% datafield - string
% dataparam - cell array of strings

parms = mmil_args2parms(varargin,{...
  'opt',[],[],...
  'tsinfotype',[],[],...
  'funspecs',[],[],...
  'verbose',1,{0,1},...
  'logfile',      [],[],...
  'logfid',       [1],[], ...     
  },false);

if ~isstruct(data)
  mmil_logstr(parms,'%s: Input must be a timesurfer data structure.\n',mfilename);
%   if parms.verbose, fprintf('%s: Input must be a timesurfer data structure.\n',mfilename); end
  datatype = []; datafield = []; dataparam = []; return;
end
if isempty(parms.funspecs)
  [parms.funspecs ext] = find_masterfile('funspecs');
end
flds = fieldnames(data); % input data structure field names
dat0 = parsecell(mmil_readtext(parms.funspecs,',','','','empty2nan'));
dat  = {};
for r = 1:size(dat0,1)
  elm1 = dat0{r,1}; % potential datatype string
  elm2 = dat0{r,2}; % potential cell array of acceptable field names
  elm3 = dat0{r,3}; % potential cell array of acceptable param names
  % look for datafield
  if ~any(isnan(elm1)) && ~any(isnumeric(elm2))
    [id jk] = match_str(elm2,flds);
    if length(id)>1
      if parms.verbose, warning('found >1 datafield match; defaulting to first match.'); end
      id = id(1);
    end
    if ~isempty(id) % found a field name match in col2 of this row
      if ~iscell(elm2), elm2 = {elm2}; end
      elm2 = elm2{id};
      % look for dataparam 
      [id jk] = match_str(elm3,fieldnames(data.(elm2)));
      if ~isempty(id)
        % found a dataparam match in col3 of this row
        dat(end+1,:) = dat0(r,:); % keep this row since it matches datafield & dataparam
        if ~iscell(dat{end,3}), dat{end,3} = {dat{end,3}}; end
        dat{end,3}   = dat{end,3}(id); % remove params not in data
      end
    end
  end
end
clear dat0

if isempty(dat)
  mmil_logstr(parms,'Could not find object metadata in: %s\n',parms.funspecs);
%   if parms.verbose, fprintf('Could not find object metadata in: %s\n',parms.funspecs); end
  datatype = []; datafield = []; dataparam = []; return;
end

% find the unique datatype
[datatypes i j] = unique(dat(:,1));
if length(datatypes)>1
  mmil_logstr(parms,'More than one object has matching metadata in: %s\n',parms.funspecs);
%   if parms.verbose, fprintf('More than one object has matching metadata in: %s\n',parms.funspecs); end
  datatype = []; datafield = []; dataparam = []; return;  
end
datatype  = dat{1,1};
datafield = dat{1,2};
dataparam = dat{1,3};
if ~iscell(dataparam), dataparam = {dataparam}; end
