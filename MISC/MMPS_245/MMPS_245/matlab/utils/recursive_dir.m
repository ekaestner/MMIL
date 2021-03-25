function flist = recursive_dir(d)

flist = [];

files = dir(d);
if isempty(files)
  return
end

for i=1:length(files)
  dirname = files(i).name;
  if  files(i).isdir 
    if ~strcmp(dirname,'.') & ~strcmp(dirname,'..')
      flist = [flist, recursive_dir(fullfile(d,dirname))]; % recursive calling of this function.
    end
  else
    flist = [flist, {fullfile(d,dirname)}];
  end
end
