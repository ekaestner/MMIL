function rules = mmil_read_classify_rules(fname_rules)
%function rules = mmil_read_classify_rules(fname_rules)
%
% Purpose: read csv file containing dicom classification rules
%
% Required Input:
%   fname_rules: full path name of csv spreadsheet
%
%   Notes:
%     rows of fname_rules correspond to different series types to be classified
%     column headers name the rules for classification
%     must have one column called "SeriesType"
%     other column headers should correspond to field names in SeriesInfo
%       e.g. SeriesDescription, SequenceName, etc. (see mmil_classify_dicoms)
%     column headers may have a prefix string defining the operation
%       i.e. 'match', 'matchi', 'exact', 'expr', or 'set'
%     if no prefix is present, the 'set' operation is assumed
%
%   Operations:
%     'match'  : use regexp to match a part of the string
%     'matchi' : use regexpi to match string, ingnoring case
%     'exact' : use strcmp to match the entire string or == to match a value
%     'expr'  : use eval to evaluate a logical expresssion (e.g. >0, <100, >=5)
%               NOTE: value from SeriesInfo always on left side of expression
%     'set'   : set SeriesInfo field with this value
%
% Output:
%   rules: struct array containing rules for each series in fname_rules
%     organized by operation
%
%   e.g. rules(1).SeriesType = 'MPR'
%        rules(1).match.SeriesDescription = 'MPRAGE'
%        rules(1).exact.SequenceName = 'EFGRE3D'
%        rules(1).expr.FlipAngle = '<10'
%        rules(1).set.rfrxcoiltype = 'BODY'
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 11/16/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
rules = [];

operations = {'match','matchi','exact','expr','set'};
nops = length(operations);

if ~exist(fname_rules,'file')
  error('file %s not found',fname_rules);
end;

tmp_rules = mmil_csv2struct(fname_rules);

if ~isfield(tmp_rules,'SeriesType')
  error('fname_rules %s does not contain ''SeriesType'' column',fname_rules);
end;
nseries = length(tmp_rules);

rule_strings = setdiff(fieldnames(tmp_rules),'SeriesType');
nrules = length(rule_strings);

rule_operations = cell(nrules,1);
rule_varnames = cell(nrules,1);
for r=1:nrules
  operation = 'set';
  varname = rule_strings{r};
  [tok,rem] = strtok(varname,'_');
  if ismember(tok,operations)
    if length(rem)<=1
      error('missing variable name in rule %s',varname);
    end;
    operation = tok;
    varname = rem(2:end);
  end;
  rule_operations{r} = operation;
  rule_varnames{r} = varname;
end;

for s=1:nseries
  if isempty(tmp_rules(s).SeriesType)
    rules(s).SeriesType = 'JUNK';
  else
    rules(s).SeriesType = tmp_rules(s).SeriesType;
  end;
  for j=1:nops
    rules(s).(operations{j}) = [];
  end;
  for r=1:nrules
    rulestr = rule_strings{r};
    operation = rule_operations{r};
    varname = rule_varnames{r};
    if ~isempty(tmp_rules(s).(rulestr))
      rules(s).(operation).(varname) = tmp_rules(s).(rulestr);
    end;
  end;
end;

