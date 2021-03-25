function D = recursive_dir_info(d)

D = [];

files = dir(d);
if isempty(files)
  return
end

for i=1:length(files)
  dirname = files(i).name;
  if files(i).isdir 
    if ~strcmp(dirname,'.') & ~strcmp(dirname,'..')
      D = [D, recursive_dir_info(fullfile(d,dirname))]; % recursive calling of this function.
    end
  else
    tmp_D = [];
    tmp_D.name = fullfile(d,dirname);
    tmp_D.bytes = files(i).bytes;
    tmp_D.date = files(i).date;
    tmp_D.datenum = files(i).datenum;
    D = [D, tmp_D];
  end
end
