function SeriesInfo = mmil_classify_by_rules(SeriesInfo,rules)
%function SeriesInfo = mmil_classify_by_rules(SeriesInfo,rules)
%
% Purpose: classify series of dicom files based on set of rules
%
% Required Input:
%   SeriesInfo: struct array containing information about
%     multiple dicom image series returned from mmil_classify_dicoms
%   rules: struct array containing rules for each series in fname_rules
%     organized by operation
%
%   Operations:
%     'match'  : use regexp to match a part of the string
%     'matchi' : use regexpi to match string, ingnoring case
%     'exact'  : use strcmp to match the entire string or == to match a value
%     'expr'   : use eval to evaluate a logical expresssion (e.g. >0, <100, >=5)
%                NOTE: value from SeriesInfo always on left side of expression
%     'set'    : set SeriesInfo field with this value
%
%   e.g. rules(1).SeriesType = 'MPR'
%        rules(1).match.SeriesDescription = 'MPRAGE'
%        rules(1).exact.SequenceName = 'EFGRE3D'
%        rules(1).expr.FlipAngle = '<10'
%        rules(1).set.rfrxcoiltype = 'BODY'
%
% Output:
%   SeriesInfo: will contain updated SeriesType field plus
%     additional information set from rules
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 04/26/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

operations = {'match','matchi','exact','expr','set'};
nops = length(operations);

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore, continue; end;
  for r=1:length(rules)
    match_flag = 1;
    for j=1:nops
      operation = operations{j};
      tmp_rule = rules(r).(operation);
      if isempty(tmp_rule), continue; end;
      varnames = fieldnames(tmp_rule);
      nvars = length(varnames);
      for v=1:nvars
        varname = varnames{v};
        rule_value = rules(r).(operation).(varname);
        if ~strcmp(operation,'set')
          if isempty(regexp(varname,'_dot_'))
            series_value = mmil_getfield(SeriesInfo(s),varname);
          else
            series_value = get_dot_value(SeriesInfo,s,varname);
          end;
          if isempty(series_value), match_flag = 0; break; end;
          if isstr(series_value)
            series_value = mmil_rowvec(series_value);
          end;
        end;
        switch operation
          case 'match'
            if isempty(regexp(series_value,rule_value))
              match_flag = 0; break;
            end;
          case 'matchi'
            if isempty(regexpi(series_value,rule_value))
              match_flag = 0; break;
            end;
          case 'exact'
            if isnumeric(series_value) || islogical(series_value)
              if ischar(rule_value), rule_value = str2double(rule_value); end;
              if islogical(series_value), series_value = 1.0 * series_value; end;
              if series_value ~= rule_value, match_flag = 0; break; end;
            else
              if isnumeric(rule_value), rule_value = num2str(rule_value); end;
              if ~strcmp(series_value,rule_value), match_flag = 0; break; end;
            end;
          case 'expr'
            if isnumeric(series_value)
              series_value = num2str(series_value);
            end;
            try
              expr = [series_value ' ' rule_value];
            catch
              error('unable to construct expression for "%s":\n%s',...
                varname,lasterr);
            end;
            try
              res = eval(expr);
            catch
              error('unable to evaluate expression "%s":\n%s',...
                expr,lasterr);
            end;
            if res~=1, match_flag = 0; break; end;
          case 'set'
            if isempty(regexp(varname,'_dot_'))
              SeriesInfo(s).(varname) = rule_value;
            else
              SeriesInfo = set_dot_value(SeriesInfo,s,varname,rule_value);
            end;
          otherwise
            error('invalid rule operation');
        end;
      end;
      if ~match_flag, break; end;
    end;
    if ~match_flag, continue; end;
    SeriesInfo(s).SeriesType = rules(r).SeriesType;
    break;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = get_dot_value(SeriesInfo,s,varname)
  % replace '_dot_' in varname
  try
    value = eval(['SeriesInfo(s).(''' regexprep(varname,'_dot_',''').(''') ''')']);
  catch
    error('failed to get SeriesInfo value for %s',varname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SeriesInfo = set_dot_value(SeriesInfo,s,varname,value)
  % replace '_dot_' in varname
  lhs = ['SeriesInfo(s).(''' regexprep(varname,'_dot_',''').(''') ''')'];
  % add quotes around a string or convert number to string
  if ischar(value)
    rhs = ['''' value ''''];
  elseif isnumeric(value)
    rhs = num2str(value);
  else
    error('set value is neither char nor numeric');
  end;
  try
    eval([lhs ' = ' rhs ';']);
  catch
    error('failed to set SeriesInfo value for %s',varname);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

