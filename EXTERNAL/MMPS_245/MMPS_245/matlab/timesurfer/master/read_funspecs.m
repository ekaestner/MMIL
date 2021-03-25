function [flags,types,cmd] = read_funspecs(fun,parms)

if ~exist('parms','var')
  parms = [];
end

funspecfile = find_masterfile('funspecs');
funtypefile = find_masterfile('funtypes');

[funspecs,specres] = mmil_readtext(funspecfile, ',','','','empty2zero');
[funtypes,typeres] = mmil_readtext(funtypefile, ',','','','empty2zero');

funspecs = funspecs(specres.stringMask(:,1),1:5);
funtypes = funtypes(typeres.numberMask(:,1),1:8);

specs = funspecs(ismember(funspecs(:,1),fun),:);

% funspecs.csv: iflag & oflag
try iflag = specs{4}; catch iflag = 0; end
try oflag = specs{5}; catch oflag = 0; end

% funtypes.csv: loadflag & saveflag
typid = [funtypes{:,2}]==iflag & [funtypes{:,3}]==oflag;
try lflag = funtypes{typid,4}; catch lflag = 0; end
try sflag = funtypes{typid,5}; catch sflag = 0; end

% funspecs.csv & parms: itype & otype
try itype = specs{2}; catch itype = 0; end
try otype = specs{3}; catch otype = 0; end

if ischar(itype)
  itype = typeselect(parms,itype,'input',  1);
  itype = typeselect(parms,itype,'itype',  0);
end
if ischar(otype)
  otype = typeselect(parms,otype,'output', 1);
  otype = typeselect(parms,otype,'otype',  0);
end

% funtypes.csv: function call
istr = funtypes{typid,6};
ostr = funtypes{typid,7};

istr(regexp(istr,'\w\s\w')+1)=',';
if strcmp(ostr,'NaN') || isempty(ostr) || ~ischar(ostr), ostr = ''; end

if length(regexp(istr,';'))>1
  ind = regexp(istr,';');
  ind = ind(end-1);
  cmd = [istr(1:ind+1) ostr istr(ind+2:end)];
else
  cmd = [ostr ' ' istr];
end

cmd = strrep(cmd,'fun',fun);
if ischar(itype), cmd = strrep(cmd,'itype',itype); end
if ischar(otype), cmd = strrep(cmd,'otype',otype); end

flags.iflag     = iflag;
flags.oflag     = oflag;
flags.load_flag = lflag;
flags.save_flag = sflag;
types.itype     = itype;
types.otype     = otype;

function type = typeselect(parms,types,typefield,datcheck)
if ismember('{',types) && ismember('}',types)
  type  = [];
  types = eval(types);
  if isfield(parms,typefield) && ischar(parms.(typefield)) && ismember(parms.(typefield),types)
    type = parms.(typefield);
  elseif datcheck && isfield(parms,'hdr') && exist(parms.hdr,'file')
    load(parms.hdr);
    if issubfield(hdr,'parms.filename')
      if ~iscell(hdr.parms.filename), hdr.parms.filename = {hdr.parms.filename}; end
      if exist(hdr.parms.filename{1},'file')
        S = who('-file',hdr.parms.filename{1}); % cell array of variable names
        S = S{1}; % pick the first name
        if ismember(S,types)
          type = S;
        end
      end
    end
  elseif datcheck && isfield(parms,'datafile') && ~isempty(parms.datafile)
    if ~iscell(parms.datafile), parms.datafile = {parms.datafile}; end
    if exist(parms.datafile{1},'file')
      S = who('-file',parms.datafile{1}); % cell array of variable names
      S = S{1};
      if ismember(S,types)
        type = S;
      end    
    end
  end
  if isempty(type)
    type = types{1};
  end
else
  type = types;
end