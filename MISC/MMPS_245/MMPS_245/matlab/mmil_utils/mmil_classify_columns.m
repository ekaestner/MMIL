function [vals,ind_num,ind_str] = classify_columns(vals,labels,verbose)
%function [vals,ind_num,ind_str] = classify_columns(vals,labels,verbose)
%
% Purpose: classify the columns of a cell array as numeric or string
%
% Required Input:
%   vals: cell matrix of values
%
% Optional Input:
%   labels: cell array of column labels
%     {default = []}
%   verbose: [0|1] display warning messages
%     {default = 1}
%
% Output:
%   vals: cell matrix of values with empty values replaced with NaNs
%           and values in each column forced to be all numeric or all string
%   ind_num: vector of indices of numeric columns
%   ind_str: vector of indices of string columns
%
% Created:  02/24/14 by Don Hagler
% Last Mod: 02/24/14 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('labels','var'), labels = []; end;
if ~exist('verbose','var') || isempty(verbose), verbose = 0; end;

ind_num = [];
ind_str = [];
ncols = size(vals,2);
for i=1:ncols
  if isempty(labels)
    col = num2str(i);
  else
    col = labels{i};
  end;
  tmp_vals = vals(:,i);
  % replace empty cells with NaNs
  ind_empty = find(cellfun(@isempty,tmp_vals));
  if ~isempty(ind_empty)
    if verbose
      fprintf('%s: WARNING: replacing empty values in column %s with NaNs...\n',...
        mfilename,col);
    end;
    for j=1:length(ind_empty)
      tmp_vals{ind_empty(j)} = NaN;
    end;
    vals(:,i) = tmp_vals;
  end;
  % replace string NaNs with numeric NaNs
  ind_char = find(cellfun(@ischar,tmp_vals));
  if ~isempty(ind_char)
    ind_char_nan = find(strcmpi(tmp_vals(ind_char),'nan') |...
                        strcmpi(tmp_vals(ind_char),'na'));
    if ~isempty(ind_char_nan)
      if verbose
        fprintf('%s: WARNING: replacing string NaNs in column %s with numeric NaNs...\n',...
          mfilename,col);
      end;
      for j=1:length(ind_char_nan)
        tmp_vals{ind_char(ind_char_nan(j))} = NaN;
      end;
    end;
    vals(:,i) = tmp_vals;
  end;
  if all(cellfun(@isnumeric,tmp_vals))
    ind_num = [ind_num,i];
  elseif all(cellfun(@ischar,tmp_vals))
    ind_str = [ind_str,i];
  else
    if verbose
      fprintf('%s: WARNING: replacing numeric values in column %s with strings...\n',...
        mfilename,col);
    end;
    % convert numeric to string
    tmp_vals = cellfun(@num2str,tmp_vals,'UniformOutput',false);
    ind_str = [ind_str,i];
    vals(:,i) = tmp_vals;
  end;
end;


