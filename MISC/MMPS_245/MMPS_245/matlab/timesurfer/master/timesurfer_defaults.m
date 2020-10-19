% default_params()
DEFAULT_FILE    = '/space/monkeys/1/home/jsherfey/svn/dev/timesurfer/guidefaults.csv';
DEFAULT_UISTYLE = 'edit';
DEFAULT_DETAILS = 'advanced';
DEFAULT_HELP = 'no documentation';
% DEFAULT_MODE = 'study';
CURRENTSET   = 0;

% load defaults once at startup & keep in a global variable
% display only functions/parameters selected by the user
if exist(DEFAULT_FILE,'file')
  defaultfile  = DEFAULT_FILE;
else
  defaultfile  = which('guidefaults.csv');
end
[defaults,res] = mmil_readtext(defaultfile, ',','','','empty2NaN');
defaults       = defaults(res.stringMask(:,1) & res.stringMask(:,2),:);
defaultstr     = cellfun(@(x) num2str(x),defaults,'uniformoutput',false);
% defaultval     = parsecell(defaultstr);

clear str cfg
[funs,i,j] = unique(defaultstr(:,1),'first');
[jnk,j]    = sort(i);
funs       = funs(j);
labs       = funs;
[tmp,i]  = unique(defaultstr(:,9),'first');
ix       = ~strcmp(tmp,'NaN') & cellfun(@(x) length(x)<30,defaultstr(i,9));
tmp      = defaultstr(i(ix),1);
[s1 s2]  = match_str(funs,tmp);
tmp      = defaultstr(i(ix),9);
labs(s1) = tmp(s2);

PARAMS = [];
for f = 1:length(funs)
  fun = funs{f};
  row = ismember(defaultstr(:,1),fun);
  par = deblank(defaultstr(row,2));
  def = deblank(defaultstr(row,3));
  typ = deblank(defaultstr(row,4));
  adv = deblank(defaultstr(row,5));
  lab = deblank(defaultstr(row,6));
  txt = deblank(defaultstr(row,7));
  aux = deblank(defaultstr(row,8));
%   answer = inputdlg(par,fun,1,def);
  answer = def;
  clear uitype advanc tmp
  for p = 1:length(par)
    str{f}.(par{p}) = answer{p}; % deblank
    tmp             = parsecell(answer(p));
    cfg{f}.(par{p}) = tmp{1};
    if strcmp(typ{p},'NaN'), typ{p} = DEFAULT_UISTYLE; end
    if strcmp(adv{p},'NaN'), adv{p} = DEFAULT_DETAILS; end
    if strcmp(txt{p},'NaN'), txt{p} = DEFAULT_HELP;    end
    uitype.(par{p}) = typ{p};
    advanc.(par{p}) = adv{p};    
    labels.(par{p}) = lab{p};
    txtdoc.(par{p}) = txt{p};
    auxcon.(par{p}) = aux{p};
%     uitype(f).(par{p}) = typ{p};
%     advanc(f).(par{p}) = adv{p};
  end
  if ~isfield(cfg{f},'function_id'), cfg{f}.function_id = 'NaN'; end
  cfg{f}.function    = fun;
  cfg{f}.funlabel    = labs{f};
  PARAMS.defaults{f} = cfg{f};
  PARAMS.controls{f}.uitype   = uitype;
  PARAMS.controls{f}.level    = advanc;
  PARAMS.controls{f}.labels   = labels;
  PARAMS.controls{f}.help     = txtdoc;
  PARAMS.controls{f}.uiaux    = auxcon;
  PARAMS.controls{f}.function = fun;
end

PARAMS.currfuns         = {};
PARAMS.currfunlist      = {};
PARAMS.currfunindex     = [];
PARAMS.history.defaults = PARAMS.defaults;
% PARAMS.program.mode     = DEFAULT_MODE;

STUDY          = [];
STUDY.protocol = {};
STUDY.funlist  = {};
STUDY.funindex = [];
STUDY.options  = [];

FULLSTUDY                   = [];
FULLSTUDY.ID                = {};
FULLSTUDY.subjects          = [];
FULLSTUDY.subjects.ID       = {};
FULLSTUDY.subjects.notes          = [];
FULLSTUDY.subjects.notes.SubjID   = '';
FULLSTUDY.subjects.notes.Filename = '';
FULLSTUDY.subjects.notes.Comments = '';
FULLSTUDY.subjects.notes.Created  = '';
FULLSTUDY.subjects.notes.Modified = '';
FULLSTUDY.subjects.batches    = [];
FULLSTUDY.options             = [];
% FULLSTUDY.rules.labels      = {};
FULLSTUDY.rules.definitions = {};
FULLSTUDY.rules.constants   = {};
FULLSTUDY.rules.derivatives = {};
FULLSTUDY.rules.parameters  = {};





