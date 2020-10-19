function studyinfo = ts_load_studyinfo(fname)
%function studyinfo = ts_load_studyinfo(fname)
%
%  Purpose: load csv file containing study info for group analysis
%
%  Input:
%   fname: name of study info csv file
%
%  Output:
%   studyinfo: structure with fields containing study info in cell arrays
%      studyinfo.subjnames: subject names (optional)
%      studyinfo.groupnames: group names (optional)
%      studyinfo.datafiles: file names for each subject and condition
%        (cell matrix)
%      studyinfo.condnames: condition names
%        (1 for each column in studyinfo.datafiles)
%      studyinfo.hemis: hemisphere labels ('lh' or 'rh')
%        (1 for each column in studyinfo.datafiles)
%
%  Structure of study info csv (comma separated value) file:
%
%    Row 1 contains column headers.
%    Subsequent rows define subjects and their data files.
%    
%    1 column may contain FreeSurfer subject names,
%      and the column header should be "Subject"
%    1 column may contain group names,
%      and the column header should be "Group"
%    Multiple columns may contain the full path or relative file names
%      of data files containing surface timecourses.
%      The column header should begin with "Condition", followed
%      by the condition name, followed by cortical hemisphere ('lh' or 'rh')
%        e.g. for condition "A" and left hemisphere,
%        the column header should be "Condition A lh"
%
%      For ts_groupavg, data files should be mgh format and resampled to ico.
%
%    Row 1 may look like this (with two conditions): 
%     "Subject","Group","Condition A lh","Condition A rh",...
%                       "Condition B lh","Condition B rh"
%    Or this:
%     "Subject","Condition A lh","Condition A rh",...
%               "Condition B lh","Condition B rh"
%    Or this:
%     "Condition A lh","Condition A rh","Condition B lh","Condition B rh"
%
%    Specifying subject names is not required except for
%      batch preprocessing to sample surface stats to ico (not yet implemented)
%    Specifying groups is optional and only makes sense when there are
%      two or more groups that will be compared in subsequent analyses
%
%  Created By:        Don Hagler     01/12/2008
%  Last Modified By:  Don Hagler     01/13/2008
%

%% todo: document GLM regressors in help

verbose = 1;
studyinfo = [];

% get info from csv file
if verbose
  fprintf('%s: reading study info from %s...\n',mfilename,fname);
end;
rawinfo = mmil_readtext(fname);
colnames = lower(rawinfo(1,:));

% identify columns
subj_colnum = find(strncmp('subject',colnames,length('subject')));
group_colnum = find(strncmp('group',colnames,length('group')));
cond_colnum = find(strncmp('condition',colnames,length('condition')));
reg_colnum = find(strncmp('reg',colnames,length('reg')));

% get subject names (optional)
if length(subj_colnum)>1
  error('multiple subject columns in %s',fname);
elseif ~isempty(subj_colnum)
  if verbose
    fprintf('%s: column %d contains subject names\n',mfilename,subj_colnum);
  end;
  studyinfo.subjnames = rawinfo(2:end,subj_colnum);
else
  if verbose
    fprintf('%s: no columns with subject names\n',mfilename);
  end;
  studyinfo.subjnames = [];
end;

% get group names (optional)
if length(group_colnum)>1
  error('multiple group columns in %s',fname);
elseif ~isempty(group_colnum)
  if verbose
    fprintf('%s: column %d contains group names\n',mfilename,group_colnum);
  end;
  % remove leading or trailing spaces
  studyinfo.groupnames = strtrim(rawinfo(2:end,group_colnum));
  % replace any spaces with underscores
  studyinfo.groupnames = regexprep(studyinfo.groupnames,'\s+','_');
else
  if verbose
    fprintf('%s: no columns with group names\n',mfilename);
  end;
  studyinfo.groupnames = [];
end;

% get data files
if isempty(cond_colnum)
  error('no condition columns specified in studyinfo file %s',fname);
end;
nconds = length(cond_colnum);
if verbose
  tmp = num2str(cond_colnum);
  tmp = regexprep(tmp,'\s+',',');
  if nconds>1
    fprintf('%s: columns %s contain data files\n',mfilename,tmp);
  else
    fprintf('%s: column %s contains data files\n',mfilename,tmp);
  end;
end;
studyinfo.datafiles = rawinfo(2:end,cond_colnum);

% get condnames and hemis
conditions = colnames(cond_colnum);
studyinfo.condnames = {};
studyinfo.hemis = {};
for j=1:nconds
  n = regexp(lower(conditions{j}),'^condition (?<condname>.+) (?<hemi>[lr]h$)','names');
  if isempty(n)
    error('condition header %s does not match expected pattern e.g. "condition A lh"',...
      conditions{j});
  end;
  studyinfo.hemis{j} = n.hemi;
  % remove leading or trailing spaces
  studyinfo.condnames{j} = strtrim(n.condname);
end;
% replace any spaces with underscores
studyinfo.condnames = regexprep(studyinfo.condnames,'\s+','_');

% get regressors
if isempty(reg_colnum)
%  fprintf('%s: no regressor columns specified\n',mfilename);
%  studyinfo.regressors = [];
else
  nreg = length(reg_colnum);
  if verbose
    tmp = num2str(reg_colnum);
    tmp = regexprep(tmp,'\s+',',');
    if nregs>1
      fprintf('%s: columns %s contain regressors\n',mfilename,tmp);
    else
      fprintf('%s: column %s contains regressors\n',mfilename,tmp);
    end;
  end;
  studyinfo.regressors = cell2mat(str2double(rawinfo(2:end,reg_colnum)));

  % get regnames
  regcolname = colnames(reg_colnum);
  studyinfo.regnames = {};
  for j=1:nregs
    n = regexp(lower(regcolname{j}),'^reg (?<regname>.+)','names');
    if isempty(n)
      error('regressor header %s does not match expected pattern e.g. "reg X"',...
        regcolname{j});
    end;
    studyinfo.regnames{j} = n.regname;
  end;
end;
