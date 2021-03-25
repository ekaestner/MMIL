function files = ts_find_timesurfer_files(opt,objecttype)
if isfield(opt,'rootindir')
  directory = opt.rootindir;
else
  directory = pwd;
end

if isfield(opt,'subdir')
  directory = fullfile(directory,opt.subdir);
end

res = what(directory);
res = res.mat;

files = {};
for i=1:length(res)
  if strcmp(who('-file',res{i}),objecttype)
    files{end+1} = res{i};
  end
end