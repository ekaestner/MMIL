function cond_key = get_cond_key(conditionfile)

if exist(conditionfile,'file') && findstr(conditionfile,'.mat')
  load(conditionfile);
elseif exist(conditionfile,'file') && findstr(conditionfile,'.csv')
  ts_makecondkey(conditionfile);
  fpath = fileparts(conditionfile);
  load(fullfile(fpath,'cond_key.mat'));
else
  cond_key = [];
end


