function ts_write_condition_key(filename)
% skips the first line

[pathstr,fname,ext,vrsn] = fileparts(filename);
if strcmp(ext,'.xls')
  % excel spreadsheet
  warning('off','MATLAB:xlsread:Mode');
  [num,txt,raw]=xlsread(filename);
  warning('on','MATLAB:xlsread:Mode');
  ncond = size(raw,1)-1;
  [cond_key(1:ncond).cond] = deal(raw{2:end,1});
  [cond_key(1:ncond).event_code] = deal(raw{2:end,1});
  [cond_key(1:ncond).name] = deal(raw{2:end,3});
end
save(strrep(filename,ext,'.mat'),'cond_key');