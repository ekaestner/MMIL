function [filename ext] = find_masterfile(type)
% type: the kind of master file to find {'funspecs','defaults','filenames'}

if strcmp(type,'filename') && ~isempty(whos('global','filepartnames'))
  global filepartnames
  filename = filepartnames;
  ext      = 'csv';
else
  S = which(mfilename);
  S = S(1:find(S=='/',1,'last'));
  filename = [S type '.csv'];
  if exist(filename), ext='csv'; return; end
  filename = [S type '.mat'];
  if exist(filename), ext='mat'; return; end
  filename = ['/home/halgsvn/matlab/master/' type '.csv'];
  if exist(filename), ext='csv'; return; end
  filename = usehardcode(type);
  ext      = 'default';
end
function cellmat = usehardcode(type)
switch type
  case 'funspecs'
  case 'defaults'
  case 'filenames'
end