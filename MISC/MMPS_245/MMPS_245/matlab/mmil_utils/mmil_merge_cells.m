function [vals,ind1,ind2] = mmil_merge_cells(vals1,vals2,merge_field,merge_flag)
%function [vals,ind1,ind2] = mmil_merge_cells(vals1,vals2,merge_field,merge_flag)
%
% purpose: merge two cell matrices
%
% required input:
%   vals1: cell matrix with first row containing column headers
%   vals2: cell matrix with first row containing column headers
%
% optional input:
%   merge_field: string specifiying name of column on which to merge
%     if not supplied, will use first column header in common
%     rows with empty values in this column will be removed
%     {default = []}
%   merge_flag: how to handle missing rows in each cell array
%     0: intersection of rows in vals1 and vals2
%     1: all rows in vals1, with blanks for missing vals2
%     2: all rows in vals2, with blanks for missing vals1
%     3: union of rows in vals1 and vals2
%     {default = 0}
%
% output:
%   vals: merged cell array
%   ind1: index to rows of vals1 used for vals
%   ind2: index to rows of vals2 used for vals
%
% Created:  12/06/12 by Don Hagler
% Prev Mod: 08/13/13 by Don Hagler
% Last Mod: 01/02/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('merge_field','var'), merge_field = []; end;
if ~exist('merge_flag','var') || isempty(merge_flag), merge_flag = 0; end;
vals = []; ind1 = []; ind2 = [];

% check for bad merge_flag values
if ~ismember(merge_flag,[0:3])
  error('invalid merge_flag: %d',merge_flag');
end;

% get column headers from each cell matrix
cols1 = vals1(1,:);
cols2 = vals2(1,:);

% determine merge field if not supplied
if isempty(merge_field)
  [cols,ind1,ind2] = intersect(cols1,cols2);
  if isempty(cols)
    error('no common variable names to merge datasets');
  end;
  ind1 = sort(ind1);
  merge_field = cols1{ind1(1)};
end;

% find column index of merge_field
ind_m1 = find(strcmp(cols1,merge_field));
ind_m2 = find(strcmp(cols2,merge_field));

% find rows in common
% find columns in common, keeping original order
cols = cat(2,cols1,cols2);
[cols,ind_c] = unique(cols,'first');
[ind_c,ind_sort] = sort(ind_c);
cols = cols(ind_sort);

% find indexes back to columns of vals1
[tmp,ind_c1,ind_c_to1] = intersect(cols1,cols);
[ind_c1,ind_sort1] = sort(ind_c1);
ind_c_to1 = ind_c_to1(ind_sort1);  

% find indexes back to columns of vals2
[tmp,ind_c2,ind_c_to2] = intersect(cols2,cols);
[ind_c2,ind_sort2] = sort(ind_c2);
ind_c_to2 = ind_c_to2(ind_sort2);

% exclude redundant columns in vals2
[ind_c_to2,ind_diff] = setdiff(ind_c_to2,ind_c_to1);
ind_c2 = ind_c2(ind_diff);

% exclude merge field from vals1
ind_c1 = setdiff(ind_c1,ind_m1);

rows1 = vals1(2:end,ind_m1);
rows2 = vals2(2:end,ind_m2);

% remove rows from vals1 with empty values in merge field
ind_empty = find(cellfun(@isempty,rows1));
if ~isempty(ind_empty)
  fprintf('%s: WARNING: removing %d rows with empty %s values from vals1\n',...
    mfilename,length(ind_empty),merge_field);
  ind_nonempty = find(~cellfun(@isempty,rows1));
  vals1 = vals1([1;ind_nonempty+1],:);
  rows1 = vals1(2:end,ind_m1);
end;

% remove rows from vals2 with empty values in merge field
ind_empty = find(cellfun(@isempty,rows2));
if ~isempty(ind_empty)
  fprintf('%s: WARNING: removing %d rows with empty %s values from vals2\n',...
    mfilename,length(ind_empty),merge_field);
  ind_nonempty = find(~cellfun(@isempty,rows2));
  vals2 = vals2([1;ind_nonempty+1],:);
  rows2 = vals2(2:end,ind_m2);
end;

% remove rows from vals1 with duplicate values in merge field
[uniq_rows1,ind_u1] = unique(rows1);
if length(uniq_rows1)~=length(rows1)
  fprintf('%s: WARNING: removing %d rows with duplicate %s values from vals1\n',...
    mfilename,length(rows1)-length(uniq_rows1),merge_field);
  vals1 = vals1([1;ind_u1+1],:);
  rows1 = vals1(2:end,ind_m1);
end;

% remove rows from vals2 with duplicate values in merge field
[uniq_rows2,ind_u2] = unique(rows2);
if length(uniq_rows2)~=length(rows2)
  fprintf('%s: WARNING: removing %d rows with duplicate %s values from vals2\n',...
    mfilename,length(rows2)-length(uniq_rows2),merge_field);
  vals2 = vals2([1;ind_u2+1],:);
  rows2 = vals2(2:end,ind_m2);
end;

if merge_flag<=2
  [rows,ind_r1,ind_r2] = intersect(rows1,rows2);
else
  rows = union(rows1,rows2);
end;

tmp_vals1 = vals1;
tmp_vals2 = vals2;

% get subset of datasets
switch merge_flag
  case 0
    ind1 = [1;1+ind_r1];
    ind2 = [1;1+ind_r2];
    vals1 = tmp_vals1(ind1,ind_c1);
    vals2 = tmp_vals2(ind2,ind_c2);
  case 1
    ind1 = [1:size(vals1,1)];
    vals1 = tmp_vals1(:,ind_c1);
    ind2 = zeros(size(ind1));
    ind2([1;1+ind_r1]) = [1;1+ind_r2];
    vals2 = cell(size(vals1,1),length(ind_c2));
    vals2([1;1+ind_r1],:) = tmp_vals2([1;1+ind_r2],ind_c2);
    rows = tmp_vals1(2:end,ind_m1);
  case 2
    ind1 = [1:size(vals2,1)];
    vals2 = tmp_vals2(:,ind_c2);
    ind1 = zeros(size(ind2));
    ind1([1;1+ind_r2]) = [1;1+ind_r1];
    vals1 = cell(size(vals2,1),length(ind_c1));
    vals1([1;1+ind_r2],:) = tmp_vals1([1;1+ind_r1],ind_c1);
    rows = tmp_vals2(2:end,ind_m2);
  case 3
    [tmp_rows1,ind_r1,ind_r1_union] = intersect(rows1,rows);
    [tmp_rows2,ind_r2,ind_r2_union] = intersect(rows2,rows);
    vals1 = cell(1+length(rows),length(ind_c1));
    vals2 = cell(1+length(rows),length(ind_c2));
    vals1([1;1+ind_r1_union],:) = tmp_vals1([1;1+ind_r1],ind_c1);
    vals2([1;1+ind_r2_union],:) = tmp_vals2([1;1+ind_r2],ind_c2);
    ind1 = zeros(1,length(rows));
    ind1([1;1+ind_r1_union]) = [1;1+ind_r1];
    ind2 = zeros(1,length(rows));
    ind2([1;1+ind_r2_union]) = [1;1+ind_r2];
end;

% concatenate subsets of both datasets
vals = cat(2,vals1,vals2);

% concatenate merge field values
merge_vals = cat(1,merge_field,rows);
vals = cat(2,merge_vals,vals);

