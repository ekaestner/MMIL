function dlist = recursive_dir_regexp(d,regexpstr)

%fprintf(1,'recursive_dir_regexp(%s,%s)\n',d,regexpstr);

dlist = [];

files = dir(d);
if isempty(files)
  return
end

for i=1:length(files)
  dirname = files(i).name;
  if  isempty(regexp(dirname,regexpstr))
    if  files(i).isdir & ~strcmp(dirname,'.') & ~strcmp(dirname,'..')
      dlist = [dlist, recursive_dir_regexp(fullfile(d,dirname),regexpstr)]; % recursive calling of this function.
    end
  else
    dlist = [dlist, {fullfile(d,dirname)}];
  end
end
