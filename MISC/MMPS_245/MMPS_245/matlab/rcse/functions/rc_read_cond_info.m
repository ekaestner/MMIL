function cond_info = rc_read_cond_info(fname_conds)
%function cond_info = rc_read_cond_info(fname_conds)
%
% Created:  02/02/09 by Don Hagler
% Last Mod: 02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist(fname_conds,'file')
  error('file %s not found',fname_conds);
end;

% load cond file
raw_cond_info = mmil_readtext(fname_conds);
% convert cond info to struct array
[nrows,ncols]=size(raw_cond_info);
collabels = {raw_cond_info{1,:}};
cond_info = [];

for i=2:nrows
  tmp = [];
  for j=1:ncols
    fieldname = regexprep(collabels{j},' ','_');
    tmp = setfield(tmp,fieldname,raw_cond_info{i,j});
  end;
  if isempty(cond_info)
    cond_info = tmp;
  else
    cond_info(i-1) = tmp;
  end;
end;

for i=1:length(cond_info)
  cond_info(i).cond_number = i;
end;


